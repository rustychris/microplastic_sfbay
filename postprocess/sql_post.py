import os
import pandas as pd
import glob
import time
import six
from sqlite3 import dbapi2 as sql

import numpy as np
from stompy.grid import unstructured_grid
from stompy import utils
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from stompy.model.fish_ptm import ptm_tools

## 
g=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/splice-merge-05-filled-edit70.nc")

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

# how does it stack up for 3M particles? 60s.
##

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
     weight DOUBLE,
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

def add_ptm_run_to_db(run_dir,con):
    curs=con.cursor()

    curs.execute("INSERT into ptm_run (run_dir) VALUES (?)",
                 [run_dir])
    run_id=curs.lastrowid
    
    groups=glob.glob(os.path.join(run_dir,"*_bin.out"))
    groups.sort()

    for group in groups:
        add_ptm_group_to_db(group,run_dir,run_id,con,curs)

def add_ptm_group_to_db(group,run_dir,run_id,con,curs):
    pbf=ptm_tools.PtmBin(group)
    name=pbf.release
    curs.execute("""INSERT into ptm_group (name,filename,run_id)
                    VALUES (?,?,?)""",
                 [name,group,run_id])
    group_id=curs.lastrowid
    print(f"Run [{run_id}] {run_dir}   Group [{group_id}] {name}")

    # Record the data from release_log
    release_log=ptm_tools.ReleaseLog(group.replace('_bin.out','.release_log'))
    release=release_log.intervals
    release['group_id']=group_id
    release['epoch']=(release['time'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')
    
    rows=zip( release['epoch'],
              release['count'],
              release['group_id'],
              release['gid_min'],
              release['gid_max'] )
    
    curs.executemany(f"""
    INSERT into ptm_release (time,count,group_id,start_id,end_id)
      VALUES (?,?,?,?,?)""",
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
            print(f"insert: {i:5d} {group_id:5d}")
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
    
    for ts in utils.progress(range(Nsteps)):
        dnum,parts=pbf.read_timestep(ts)
        epoch=utils.to_unix(dnum)
        if ts%100==0:
            print(f"{ts} {dnum} {len(parts)}")
        xys.append( parts['x'][:,:2].copy() )
        gids.append( parts['id'].copy() )
        part_ids.append( [ group_gid_to_particle[(group_id,gid)]
                           for gid in parts['id'] ] )
        epochs.append( epoch*np.ones(len(parts['id']),np.int32) )

    all_xy=np.concatenate(xys)
    # all_gid=[int(i) for i in np.concatenate(gids)]
    all_part_id=[pid for sublist in part_ids for pid in sublist] # flatten
    all_epoch=[int(t) for t in np.concatenate(epochs)]
    
    # How does that compare to what I can do in straight python?
    # 42 seconds.
    t=time.time()
    # be sure to insert these as regular int
    all_cell=[int(c) for c in g.points_to_cells(all_xy)]
    elapsed=time.time() - t
    print("Python mapping time: %.3fs"%elapsed)

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
fn="ptm_and_grid.db"
clean=True

six.moves.reload_module(ptm_tools)

if clean :
    os.path.exists(fn) and os.unlink(fn)

con = sql.connect(fn)
if use_spatial:
    con.enable_load_extension(True)
    con.execute('SELECT load_extension("mod_spatialite");')
    con.execute('SELECT InitSpatialMetadata()')

if clean:
    add_grid_to_db(g,con)
    init_ptm_tables(con)

clean_ptm_tables(con)

add_ptm_run_to_db("/home/rusty/src/sfb_ocean/ptm/all_sources/all_source_select_w_const",
                  con)

# That's running. see where it gets us.
# maybe need more commits?

##

# ok - so before going any further with the details like weights and z-coord,
# how does this do with queries?

con = sql.connect(fn)
curs=con.cursor()

## 
def db(sql,*a,row_limit=5):
    t=time.time()
    curs.execute(sql,*a)
    results=curs.fetchall()
    elapsed=time.time() -t
    print("Query time: %.2fs"%elapsed)
    print("Returned %d rows"%len(results))
    if row_limit:
        print(results[:row_limit])
    else:
        print(results)

##
# a few seconds
curs.execute("ANALYZE")
con.commit()

##

# index on time?
# takes 115s.
db("CREATE INDEX IF NOT EXISTS particle_time_idx on particle_loc (time);")
db("ANALYZE") # 13s

## 
# 187_585_005  particle locations.
# 118_877_472  on the second go-around. failed earlier.
# this query was a bit slow, like 10s.
db("select count(1) from particle_loc")

## 
time_start=int( utils.to_unix( np.datetime64("2017-07-15") ) )
time_stop =int( utils.to_unix( np.datetime64("2017-07-16") ) )

# 4_146_976
# query time: 13s w/o index.
# 0.12s with index.
db("""
   select count(1) from particle_loc
    where time>=? and time < ?""",
   [time_start,time_stop])

##

# closer to the real deal. 0.6s
db("""
   select min(loc.cell),count(1) 
    from particle_loc as loc
    where loc.time>=? and loc.time < ?
      and loc.cell>=0
    group by loc.cell""",
   [time_start,time_stop])

##

# also 0.6s
db("""
   select min(loc.cell),count(1)
    from particle_loc as loc, particle as p
    where loc.time>=? and loc.time < ?
      and loc.particle_id=p.id
      and loc.cell>=0
    group by loc.cell""",
   [time_start,time_stop])

##
# db("update ptm_release set weight=1.0")

# basically the real deal. 0.86s for 24 hour
# for 15 days, that becomes 12s.
db("""
   select min(loc.cell),sum(rel.weight)
    from particle_loc as loc, particle as p, ptm_release as rel
    where loc.time>=? and loc.time < ?
      and loc.particle_id=p.id
      and p.release_id=rel.id
      and loc.cell>=0
    group by loc.cell""",
   [time_start,time_start+15*86400])

## 
db("""
   explain query plan
   select min(loc.cell),sum(rel.weight)
    from particle_loc as loc, particle as p, ptm_release as rel
    where loc.time>=? and loc.time < ?
      and loc.particle_id=p.id
      and p.release_id=rel.id
      and loc.cell>=0
    group by loc.cell""",
   [time_start,time_start+15*86400],
   row_limit=None)

##

# This is a pretty close to a real query, and takes about 2s.
# the data load was incomplete, so it doesn't reflect the full
# size of the dataset.
# looks like 36/54 groups were loaded.
# total bin_out size is 5.6G.
# so the database currently holds about 3.73G worth of bin.out
# database size is 4 222 980 096. so we're slightly less efficient
# than the binary data.
curs.execute("""
   select min(loc.cell),sum(rel.weight)
    from particle_loc as loc, particle as p, ptm_release as rel
    where loc.time>=? and loc.time < ?
      and loc.particle_id=p.id
      and p.release_id=rel.id
      and loc.cell>=0
    group by loc.cell""",
   [time_start,time_stop])

data=np.array( curs.fetchall() )

cell_count=np.zeros(g.Ncells(),np.float64)
cell_count[data[:,0].astype(np.int32)]=data[:,1]

cell_2dconc=cell_count/g.cells_area()

##
plt.figure(1).clf()

ccoll=g.plot_cells(values=cell_2dconc.clip(1e-5,np.inf),
                   cmap='jet',norm=LogNorm())
ccoll.set_edgecolor('face')
plt.axis('equal')

##
run_dir="/home/rusty/src/sfb_ocean/ptm/all_sources/all_source_select_w_const"
groups=glob.glob(os.path.join(run_dir,"*_bin.out"))
