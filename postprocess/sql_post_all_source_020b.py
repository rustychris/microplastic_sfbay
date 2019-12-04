"""
Post-processing for the newer all_source_020b PTM runs
"""
from sql_post import *


base_dir="/opt2/sfb_ocean/ptm/all_source_020b"

import os
import glob
import numpy as np

all_runs=glob.glob(os.path.join(base_dir,"*/201[78]????"))
# group by period
periods=np.unique( [ r.split('/')[-1] for r in all_runs] )
# and within each period enumerate sources
sources=np.unique( [ r.split('/')[-2] for r in all_runs] )

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
            continue

        clean=False

        if clean :
            os.path.exists(fn) and os.unlink(fn)
            create=True
        else:
            create=not os.path.exists(fn)

        for run_idx,run_dir in enumerate(run_dirs):
            run=post.PtmRun(run_dir=run_dir)
            grid=run.grid()

            con = sql.connect(fn)
            if use_spatial:
                con.enable_load_extension(True)
                con.execute('SELECT load_extension("mod_spatialite");')
                con.execute('SELECT InitSpatialMetadata()')

            if run_idx==0:
                if create:
                    add_grid_to_db(grid,con)
                    init_ptm_tables(con)
                elif clean:
                    # stale logic
                    clean_ptm_tables(con)

            add_ptm_run_to_db(run,con,grid)
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
