"""
Post-processing for the newer all_source_021b PTM runs

adapted from sql_post_all_source_020b.py
"""
from sql_post import *

base_dir="/opt2/sfb_ocean/ptm/all_source_021b"

import os
import glob
import numpy as np

# Only include sources that have run all the way through
# 20180416
all_srcs=[os.path.dirname(f)
          for f in glob.glob(os.path.join(base_dir,"*/20170720"))]
complete_srcs=[os.path.dirname(f)
               for f in glob.glob(os.path.join(base_dir,"*/20180416"))]
all_runs=[f
          for src_dir in complete_srcs
          for f in glob.glob(os.path.join(src_dir,"201[78]????") ) ]

print(f"Processing {len(complete_srcs)} completed sources of {len(all_srcs)} total, for {len(all_runs)} runs")

# group by period
periods=np.unique( [ r.split('/')[-1] for r in all_runs] )
# and within each period enumerate sources
sources=np.unique( [ r.split('/')[-2] for r in all_runs] )

clean=False # True => existing databases are removed first
update=True # True => allow partial update of existing databases

for period in periods:
    # Try shoving all of each 10 day period into the same database
    run_dirs=[os.path.join(base_dir,source,period)
              for source in sources ]

    # target database
    fn=os.path.join(base_dir,period,"ptm_and_grid.db")

    os.makedirs(os.path.dirname(fn),exist_ok=1)

    if 1: # building the original databases
        if os.path.exists(fn):
            # for the moment, play it safe and skip anything that appears to
            # have been run already.
            # But warn if it looks stale:
            fn_mtime=os.stat(fn).st_mtime

            input_mtime=0
            newer=[]
            for rd in run_dirs:
                for bin_out in glob.glob(os.path.join(rd,'*_bin.out')):
                    if os.stat(bin_out).st_mtime > fn_mtime:
                        newer.append(rd)
                        print(f"Run {rd} appears newer than {fn}")
                        break # No need to look at more in this run
            if not update:
                print("Update is not set and database exists, so moving on")
                continue

        if clean:
            os.path.exists(fn) and os.unlink(fn)
            create=True
        else:
            create=not os.path.exists(fn)

        # Create the file quickly to keep subsequent invocations from
        # going here
        sql.connect(fn).close()
        assert os.path.exists(fn)

        for run_idx,run_dir in enumerate(run_dirs):
            run=post.PtmRun(run_dir=run_dir)
            grid=run.grid()

            con = sql.connect(fn)
            if run_idx==0:
                if create:
                    add_grid_to_db(grid,con)
                    init_ptm_tables(con)

            add_ptm_run_to_db(run,con,grid,on_exists='skip') 
            # open and close each time to catch IO errors sooner
            con.commit()
            con.close()

    if 1: # building index on particle time
        print(fn)
        if not os.path.exists(fn):
            print("Can't build index - %s doesn't exist"%fn)
            continue

        con = sql.connect(fn)
        curs=con.cursor()
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

        # index on time?
        db("CREATE INDEX IF NOT EXISTS particle_time_idx on particle_loc (time);")
        db("ANALYZE") #
                
        con.commit()
        con.close()
