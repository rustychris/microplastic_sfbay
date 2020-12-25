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

import postprocess_v00 as post

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
    # note that the unique constraint was added late in this code, so 021b databases and
    # before do not have this.
    curs.execute(f"""
    CREATE TABLE ptm_run (
     id integer primary key,
     run_dir text not null,
     UNIQUE(run_dir)
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

def add_ptm_run_to_db(run,con,grid,z_extractor=None,on_exists='skip',profile=False):
    """
    run: a postprocess_v00.PtmRun instance
    con: database connection
    grid: the corresponding grid (might be supplanted)
    z_extractor: instance of EtaExtractor.

    on_exists: 'skip': if the run is already present in the database, do nothing
      may add an option in the future to clean out the old run.  And may need
      other options to allow more flexible matching of run_dir, since it is currently
      stored as a full, absolute path.
      'error':  raise an error
      'continue': if a run already exists, move on to check its groups

    return: 1 if database was modified, 0 if not.
    """
    curs=con.cursor()
    # Concurrency: wrap these steps in a transaction. If two processes hit this section
    # at the same time, both will get a SHARED lock at the time of the first SELECT
    # below. Then nobody can get a RESERVED lock for the INSERT, we'll fail, and move on.
    # Or one process is already to the select, gets the RESERVED lock, and the late-comer
    # cannot get a SHARED lock.
    modified=0

    try:
        curs.execute("BEGIN TRANSACTION")
        try:
            run_id=None
            try:
                existing=curs.execute("select id from ptm_run where run_dir=?",[run.run_dir]).fetchall()
            except sql.OperationalError:
                print("While checking for ptm_run, got error. Assuming locked, and we should move on")
                return modified

            if len(existing):
                if on_exists=='skip':
                    print(f"Run {run.run_dir} already in database. Skipping")
                    return modified
                elif on_exists=='error':
                    raise Exception(f"Run {run.run_dir} already in database.")
                elif on_exists=='continue':
                    run_id=existing[0][0]
                    print(f"Run {run.run_dir} already in database. Will proceed with run_id={run_id}")
                else:
                    raise Exception(f"Bad value for on_exists='{on_exists}'")

            if run_id is None:
                curs.execute("INSERT into ptm_run (run_dir) VALUES (?)",
                             [run.run_dir])
                run_id=curs.lastrowid
                modified+=1
        finally:
            # To help coordinate between multiple workers, commit ptm_run quickly.
            # Is this necessary?
            curs.execute("END TRANSACTION")
            con.commit()
    except sql.OperationalError as exc:
        print("Operational error (%s) while inserting run.  Assume locked."%str(exc))
        return 0

    if z_extractor is None:
        z_extractor=post.EtaExtractor(ptm_runs=[run])
        
    for group in run.groups():
        if profile:
            print("PROFILING!")
            import cProfile
            from pstats import SortKey
            pr = cProfile.Profile()
            pr.enable()
        curs.execute("BEGIN TRANSACTION")
        try:
            modified+=add_ptm_group_to_db(group,run,run_id,con,curs,grid,z_extractor=z_extractor,
                                          on_exists='skip')
        finally:
            curs.execute("END TRANSACTION")
            con.commit()
            
        if profile:
            print("END PROFILING")
            pr.disable()
            pr.print_stats(SortKey.CUMULATIVE)
            raise Exception("Will bail out -- check profiling")
        
    return modified

def add_ptm_group_to_db(group,run,run_id,con,curs,grid,z_extractor=None,
                        on_exists='error'):
    """
    returns: 1 if database was modified. 0 otherwise.
    """
    modified=0
    t0=time.time()
    
    pbf=run.open_binfile(group)

    name=group
    existing=curs.execute("""select id from ptm_group
                              where filename=?
                                and run_id=? """,
                          [name,run_id]).fetchall()
    if len(existing):
        group_id=existing[0][0]
        print(f"Run [{run_id}] {run.run_dir}   Group (existing) [{group_id}] {name}")
        if on_exists=='error':
            raise Exception("Group already exists in database")
        elif on_exists=='skip':
            return modified
        else:
            raise Exception(f"on_exists='{on_exists}' not understood in add_ptm_group_to_db")
    else:
        curs.execute("""INSERT into ptm_group (name,filename,run_id)
                        VALUES (?,?,?)""",
                     [name,name,run_id])
        group_id=curs.lastrowid
        print(f"Run [{run_id}] {run.run_dir}   Group [{group_id}] {name}")
        modified+=1

    # Record the data from release_log
    release_log=run.open_release_log(group)
    release=release_log.intervals
    release['group_id']=group_id
    release['epoch']=(release['time'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

    # Ham-handed fix up of truncated release logs.
    typical_count=int( np.median( release['count'] ) )
    bad=release['count']!=typical_count
    # a little bit smart -- only step in when it's the last step that's
    # different.
    if (not np.any(bad[:-1])) and bad[-1]:
        print("Yikes - release_log might be missing some particles.")
        print("  Changing reported count of %d to %d"%(release['count'].values[-1],
                                                       typical_count))
        print("  And punting on id_max,gid_max")
        release.loc[bad,'count']=typical_count
        # careful of inclusive indexes
        release.loc[bad,'id_max']=release.loc[bad,'id_min']+typical_count-1
        release.loc[bad,'gid_max']=release.loc[bad,'gid_min']+typical_count-1

    # add in volume information.
    Qfunc=run.get_Qfunc_for_group(group)
    # no negative flows, which can happen with SJ river
    Q=Qfunc(release['time'].values).clip(0,np.inf)

    # Note that the volume is for the full release. Downstream code
    # must divide by count (i.e. 5) in order to get the per-particle volume.
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
    modified+=curs.rowcount
    
    print("add_ptm_group_to_db: Populating unique particles based on releases")
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
            modified+=1
            
    # Assemble these one group at a time.
    Nsteps=pbf.count_timesteps()
    gids=[] # the group-unique ids
    part_ids=[] # the globally unique ids
    epochs=[]
    xys=[]

    t_elapsed=time.time()-t0
    t0+=t_elapsed
    print(f"add_ptm_group_to_db: prep time {t_elapsed:.2f}s")
    
    # Fastest to compute cells for all of the points at once, then cache the
    # result for returning by read_timesteps.  Since this mapping is grid-dependent
    # require that the grid be provided, and the returned key can be efficiently
    # used to return the correct cache data later.
    cell_fn=pbf.precompute_cells(grid)
    cell_ds=xr.open_dataset(cell_fn)
    cell_offsets=cell_ds.offset.values
    cell_ends=cell_offsets + cell_ds['count'].values
    cells=cell_ds['cell'].values
    all_cell=[int(c) for c in cells]

    t_elapsed=time.time()-t0
    t0+=t_elapsed
    print(f"add_ptm_group_to_db: precompute_cells {t_elapsed:.2f}s")
    
    z_from_beds=[]
    z_from_surfs=[]
    
    for ts in utils.progress(range(Nsteps)):
        # caching of the mapping from point to grid is 
        dnum,parts=pbf.read_timestep(ts)
        epoch=utils.to_unix(dnum)
        if ts%100==0:
            print(f"{ts} {dnum} {len(parts)}")
        # xys.append( parts['x'][:,:2].copy() )
        part_z=parts['x'][:,2].copy()
        
        eta=z_extractor.eta_for_time(utils.to_dt64(dnum))
        c=cells[cell_offsets[ts]:cell_ends[ts]]
        z_eta=eta[c]
        z_bed=grid.cells['z_bed'][c]
        part_z_from_bed=part_z-z_bed # should be >=0
        part_z_from_surf=part_z-z_eta
        z_from_beds.append(part_z_from_bed)
        z_from_surfs.append(part_z_from_surf)
        
        gids.append( parts['id'].copy() )
        part_ids.append( [ group_gid_to_particle[(group_id,gid)] 
                           for gid in parts['id'] ] )
        epochs.append( epoch*np.ones(len(parts['id']),np.int32) )

    t_elapsed=time.time()-t0
    t0+=t_elapsed
    print(f"add_ptm_group_to_db: populate z dims, particle ids {t_elapsed:.2f}s")

    # all_xy=np.concatenate(xys)
    # all_gid=[int(i) for i in np.concatenate(gids)]
    all_part_id=[pid for sublist in part_ids for pid in sublist] # flatten
    all_epoch=[int(t) for t in np.concatenate(epochs)]
    all_z_from_bed=np.concatenate(z_from_beds)
    all_z_from_surf=np.concatenate(z_from_surfs)
    
    # # How does that compare to what I can do in straight python?
    # # 42 seconds.
    # t=time.time()
    # # be sure to insert these as regular int
    # all_cell=[int(c) for c in grid.points_to_cells(all_xy)]
    # elapsed=time.time() - t
    # print("Python mapping time: %.3fs"%elapsed)

    rows=zip(all_part_id,all_epoch,all_cell,all_z_from_bed,all_z_from_surf)
    curs.executemany("""
                     INSERT INTO particle_loc (particle_id,time,cell,z_from_bed,z_from_surface)
                     VALUES (?,?,?,?,?)""",
                     rows)
    modified+=curs.rowcount
    t_elapsed=time.time()-t0
    t0+=t_elapsed
    print(f"add_ptm_group_to_db: insert into particle_loc {t_elapsed:.2f}s")
    
    # Still TODO: 
    # z_from_bed DOUBLE,
    # z_from_surface DOUBLE,

    t_elapsed=time.time()-t0
    t0+=t_elapsed
    print(f"add_ptm_group_to_db: commit {t_elapsed:.2f}s")
    return modified
        
