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
