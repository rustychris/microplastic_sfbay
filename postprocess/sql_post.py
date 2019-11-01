"""
Build database for sql-based postprocessing
"""
import os
import pandas as pd
import glob
import time
import six
from sqlite3 import dbapi2 as sql
import xarray as xr

import numpy as np
from stompy.grid import unstructured_grid
from stompy import utils
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from stompy.model.fish_ptm import ptm_tools

##

use_spatial=False

def add_grid_to_db(g,con):
    con.execute(f"""
    create table grid_cells (
      id INTEGER NOT NULL PRIMARY KEY,
      center_x DOUBLE,
      center_y DOUBLE
    );
    """)

    if use_spatial:
        con.execute(f"""
        SELECT AddGeometryColumn('grid_cells', 'poly',
                                 -1, 'POLYGON', 'XY');
        """)

    cur=con.cursor()
        
    # Insert cell data:
    ids=list(range(g.Ncells()))
    cc=g.cells_center()
    center_xs=cc[:,0]
    center_ys=cc[:,1]

    if use_spatial:
        wkbs=[g.cell_polygon(c).wkb for c in range(g.Ncells())]
        row_valuesb=zip(ids,
                        center_xs,
                        center_ys,
                        wkbs)
        cur.executemany(f"""INSERT INTO grid_cells
                             (id,center_x,center_y,poly)
                             VALUES (?,?,?,GeomFromWKB(?,-1));""",
                        row_valuesb)
    else:
        row_values=zip(ids,
                       center_xs,
                       center_ys)
        cur.executemany(f"""INSERT INTO grid_cells
                             (id,center_x,center_y)
                            VALUES (?,?,?);""",
                        row_values)

    con.commit()                



##

# So the more efficient approach is probably to do the point->cell sorting
# on the python side
# insert it all into sql
# and then go about the queries.

def init_ptm_tables(con):
    curs=con.cursor()
    curs.execute(f"""
    CREATE TABLE ptm_run (
     id integer primary key,
     run_dir text not null
    );""")

    curs.execute(f"""
    CREATE TABLE ptm_group (
     id integer primary key,
     name text,
     filename text,
     run_id integer not null,
     FOREIGN KEY(run_id) REFERENCES ptm_run(id)
    );""")

    # Record the data from release_log
    curs.execute(f"""
    CREATE TABLE ptm_release (
     id integer primary key,
     time integer not null,
     count integer not null,
     group_id integer not null,
     start_id integer not null, -- inclusive
     end_id   integer not null, -- inclusive
     -- water volume represented by a single particle from this release
     -- assuming the release is not duplicated by another release
     volume DOUBLE,
     FOREIGN KEY(group_id) REFERENCES ptm_group(id)
    );""")

    # each particle is seen potentially 1000x times.
    # is this table really needed, or could I go straight to
    # the release table?
    # going from particle_loc to release means that particle_loc
    # has to have both release_id and gid, but if I go through
    # a particle table, then only need unique particle id
    curs.execute(f"""
    CREATE TABLE particle (
     id integer primary key,
     ptm_id integer not null, -- unique only within groups.
     group_id integer not null, -- not really necessary, but helps with unique checking
     release_id integer not null,
     FOREIGN KEY(release_id) REFERENCES ptm_release(id),
     UNIQUE (ptm_id,group_id)
    );""")

    curs.execute("""
    CREATE TABLE particle_loc (
     particle_id integer not null,
     time integer not null, 
     cell integer not null,
     z_from_bed DOUBLE,
     z_from_surface DOUBLE,
     FOREIGN KEY(cell) REFERENCES grid_cells(id),
     FOREIGN KEY(particle_id) REFERENCES particle(id)
    );""")
    
    # so each particle observation in the PTM output uses 32 bytes
    # as written, each particle observation takes 3 ints, 2 doubles.
    # docs suggest that sqlite automatically chooses the width of
    # the integers.
    # so hopefully that means the above
    # 2+4+4+8+8= 26 bytes.
    # if z_from bed were instead an integer interpreted as cm,
    # could shave 2*4 bytes off that, to get 18 bytes
    # but rowid is hidden in there, for an extra 4 bytes.
    # updated approach: 3*4 + 2*8

def clean_ptm_tables(con):
    curs=con.cursor()
    curs.execute("DELETE FROM particle_loc;")
    curs.execute("DELETE FROM particle;")
    curs.execute("DELETE FROM ptm_release;")
    curs.execute("DELETE FROM ptm_group;")
    curs.execute("DELETE FROM ptm_run;")

def add_ptm_run_to_db(run,con,grid):
    """
    run: a postprocess_v00.PtmRun instance
    con: database connection
    grid: the corresponding grid (might be supplanted)
    """
    curs=con.cursor()

    curs.execute("INSERT into ptm_run (run_dir) VALUES (?)",
                 [run.run_dir])
    run_id=curs.lastrowid
    
    for group in run.groups():
        add_ptm_group_to_db(group,run,run_id,con,curs,grid)

def add_ptm_group_to_db(group,run,run_id,con,curs,grid):
    pbf=run.open_binfile(group)

    name=group
    curs.execute("""INSERT into ptm_group (name,filename,run_id)
                    VALUES (?,?,?)""",
                 [name,name,run_id])
    group_id=curs.lastrowid
    print(f"Run [{run_id}] {run.run_dir}   Group [{group_id}] {name}")

    # Record the data from release_log
    release_log=run.open_release_log(group)
    release=release_log.intervals
    release['group_id']=group_id
    release['epoch']=(release['time'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

    # HERE - add in volume information.
    Qfunc=run.get_Qfunc_for_group(group)
    # no negative flows, which can happen with SJ river
    Q=Qfunc(release['time'].values).clip(0,np.inf)
    
    grp_hour_per_rel=(release['time'].values[1] - release['time'].values[0])/np.timedelta64(3600,'s')
    # m3/s * s/hour * hour/release => m3/release
    volume=Q*3600 * grp_hour_per_rel
    
    rows=zip( release['epoch'],
              release['count'],
              release['group_id'],
              release['gid_min'],
              release['gid_max'],
              volume )
    
    curs.executemany(f"""
    INSERT into ptm_release (time,count,group_id,start_id,end_id,volume)
      VALUES (?,?,?,?,?,?)""",
                     rows)

    print("Populating unique particles based on releases")
    # Pull the releases back out to help populate particle_id
    group_gid_to_particle={}
    
    releases=curs.execute("""select id,start_id,end_id from ptm_release 
                              where group_id=?
                              order by start_id""",[group_id]).fetchall()
    # create the particles
    for rel in releases:
        for i in range(rel[1],rel[2]+1):
            # print(f"insert: {i:5d} {group_id:5d}")
            curs.execute("""INSERT INTO particle (ptm_id,group_id,release_id)
                             VALUES (?,?,?);""",
                         (i,group_id,rel[0]) )
            particle_id=curs.lastrowid
            # This mapping can then be used to fill in the remaining data below
            # when populating particle_loc
            group_gid_to_particle[ (group_id,i) ]=particle_id
                     
    # Assemble these one group at a time.
    Nsteps=pbf.count_timesteps()
    gids=[] # the group-unique ids
    part_ids=[] # the globally unique ids
    epochs=[]
    xys=[]

    # Fastest to compute cells for all of the points at once, then cache the
    # result for returning by read_timesteps.  Since this mapping is grid-dependent
    # require that the grid be provided, and the returned key can be efficiently
    # used to return the correct cache data later.
    cell_fn=pbf.precompute_cells(grid)
    cell_ds=xr.open_dataset(cell_fn)
    all_cell=[int(c) for c in cell_ds['cell'].values]
    
    for ts in utils.progress(range(Nsteps)):
        # caching of the mapping from point to grid is 
        dnum,parts=pbf.read_timestep(ts)
        epoch=utils.to_unix(dnum)
        if ts%100==0:
            print(f"{ts} {dnum} {len(parts)}")
        # xys.append( parts['x'][:,:2].copy() )
        gids.append( parts['id'].copy() )
        part_ids.append( [ group_gid_to_particle[(group_id,gid)]
                           for gid in parts['id'] ] )
        epochs.append( epoch*np.ones(len(parts['id']),np.int32) )

    # all_xy=np.concatenate(xys)
    # all_gid=[int(i) for i in np.concatenate(gids)]
    all_part_id=[pid for sublist in part_ids for pid in sublist] # flatten
    all_epoch=[int(t) for t in np.concatenate(epochs)]
    
    # # How does that compare to what I can do in straight python?
    # # 42 seconds.
    # t=time.time()
    # # be sure to insert these as regular int
    # all_cell=[int(c) for c in grid.points_to_cells(all_xy)]
    # elapsed=time.time() - t
    # print("Python mapping time: %.3fs"%elapsed)

    rows=zip(all_part_id,all_epoch,all_cell)
    curs.executemany("""
                     INSERT INTO particle_loc (particle_id,time,cell)
                     VALUES (?,?,?)""",
                     rows)
    
    # Still TODO: 
    # z_from_bed DOUBLE,
    # z_from_surface DOUBLE,
    try:
        con.commit()
    except sql.OperationalError as exc:
        print("Commit got an operational error")
        print(exc)
        print("Sleeping")
        time.sleep(5)
        print("Trying again")
        con.commit()
        
##

import postprocess_v00 as post

run=post.PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source_020/20170715/w0.0")

grid=run.grid()

fn=os.path.join(run.run_dir,"ptm_and_grid.db")

clean=True

if clean :
    os.path.exists(fn) and os.unlink(fn)

con = sql.connect(fn)
if use_spatial:
    con.enable_load_extension(True)
    con.execute('SELECT load_extension("mod_spatialite");')
    con.execute('SELECT InitSpatialMetadata()')

if clean:
    add_grid_to_db(grid,con)
    init_ptm_tables(con)

clean_ptm_tables(con)

add_ptm_run_to_db(run,con,grid)
for group in run.groups():
    # flow data for assigning particle mass
    Qfunc=run.get_Qfunc_for_group(group)

##
