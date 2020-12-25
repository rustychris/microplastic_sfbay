"""
Post-processing for the newer all_source_021b PTM runs

adapted from sql_post_all_source_020b.py, then 21b, now
22a.

"""
import os, sys
import glob
import numpy as np

import time

from sql_post import *
import post_local

import argparse

parser=argparse.ArgumentParser(description="Transcribe PTM data to sqlite")
parser.add_argument("-u","--update",help="Allow updating existing database with additional runs",
                    action='store_true')
parser.add_argument("-U","--update-groups",help="Allow updating existing runs with additional groups",
                    action='store_true')
parser.add_argument("-c","--clean",help="Remove existing databases")
parser.add_argument("-p","--profile",help="Add a singe group and profile the process",
                    action='store_true')
parser.add_argument("-s","--status",help="Scan all runs and groups and print overall progress",
                    action='store_true')

args=parser.parse_args()
if args.update_groups:
    args.update=True

base_dir=post_local.all_source_base_dir

# Only include sources that have run all the way through
# 20180416
# Need to rewrite this for chunks. Want to process all runs that
# are complete.
# Will come back on the sql side to report on whether there are
# sources which do not have complete data.

all_srcs=[os.path.dirname(f)
          for f in glob.glob(os.path.join(base_dir,"*/20170720"))]
complete_srcs=[os.path.dirname(f)
               for f in glob.glob(os.path.join(base_dir,"*/20180416"))]
all_runs=[f
          for src_dir in complete_srcs
          for f in glob.glob(os.path.join(src_dir,"201[78]????") ) ]

# group by period
periods=np.unique( [ r.split('/')[-1] for r in all_runs] )
# and within each period enumerate sources
sources=np.unique( [ r.split('/')[-2] for r in all_runs] )

print(f"Processing {len(complete_srcs)} completed sources/chunks of {len(all_srcs)} total, for {len(all_runs)} runs")
print(f"{len(periods)} unique time periods")
print(f"Profiling: {args.profile}")

db_basename='ptm_and_grid.db'

if args.status:
    print("Scanning to gather status information")
    print(f"Periods ({len(periods)}):" )
    for period in periods:
        print(f"  {period}")
        
    print(f"Chunks ({len(sources)}): ")
    for source in sources:
        print(f"  {source}")
        
    periods_with_db=0
    periods_without_db=0
    periods_partial_runs=0
    periods_partial_groups=0
    periods_complete=0
    runs_stale=0
    runs_fresh=0
    runs_absent=0
    db_read_errors=0
    
    for period in periods:
        print(f"Scanning period: {period}")

        run_dirs=[os.path.join(base_dir,source,period)
                  for source in sources ]

        # target database
        fn=os.path.join(base_dir,period,db_basename)

        if not os.path.exists(fn):
            periods_without_db+=1
            runs_absent+=len(run_dirs)
            continue
        periods_with_db+=1

        fn_mtime=os.stat(fn).st_mtime

        input_mtime=0
        newer=[]
        group_count=0 # all groups 
        for rd in run_dirs:
            for bin_out in glob.glob(os.path.join(rd,'*_bin.out')):
                group_count+=1
                if os.stat(bin_out).st_mtime > fn_mtime:
                    runs_stale+=1
                else:
                    runs_fresh+=1

        # Actually query the database to see if it's complete-ish
        con = sql.connect(fn)
        try:
            try:
                curs=con.cursor()
                # Are all runs for this period accounted for?
                existing=curs.execute("select count(1) from ptm_run").fetchall()
                db_run_count=existing[0][0]
                if db_run_count!=len(run_dirs):
                    print(f"{fn}: Expecting {len(run_dirs)} runs, database has {db_run_count}")
                    periods_partial_runs+=1
                else:
                    # Are there the expected number of groups in total?
                    existing=curs.execute("select count(1) from ptm_group").fetchall()
                    db_group_count=existing[0][0]
                    if db_group_count!=group_count:
                        print(f"Expecting {group_count} groups, database has {existing[0][0]}")
                        periods_partial_groups+=1
                    else:
                        periods_complete+=1
                # something like 4 chunks * 13 sources/chunk * 7 groups/source
                # except that the last chunk doesn't get 13 sources.
            except sql.OperationalError:
                print(f"{fn} read error. probably write-locked")
                db_read_errors+=1
        finally:
            con.close()

    print()
    print(f"Periods with database: {periods_with_db}")
    print(f"  {periods_partial_runs} missing runs")
    print(f"  {periods_partial_groups} missing groups")
    print(f"  {periods_complete} complete")
    print(f"  {db_read_errors} read errors")
    print(f"Periods w/o database:  {periods_without_db}")
    print(f"Stale runs:  {runs_stale}")
    print(f"Fresh runs:  {runs_fresh}")
    print(f"Absent runs: {runs_absent}")
    
    sys.exit(0)

for period in periods:
    print(f"Processing period: {period}")
        
    # Try shoving all of each 10 day period into the same database
    run_dirs=[os.path.join(base_dir,source,period)
              for source in sources ]

    # target database
    fn=os.path.join(base_dir,period,db_basename)

    os.makedirs(os.path.dirname(fn),exist_ok=1)

    modified=0
    
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
            if not args.update:
                print("Update is not set and database exists, so moving on")
                continue
            else:
                print("Database exists and update is enabled. Proceeding with caution")

        if args.clean:
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
                    modified+=1

            # Trying out new on_exists='continue' which will add missing groups
            # could also use 'continue', which would get missing groups.
            if args.update_groups:
                on_exists='continue'
            else:
                on_exists='skip'
            modified+=add_ptm_run_to_db(run,con,grid,
                                        on_exists=on_exists,profile=args.profile)
            # open and close each time to catch IO errors sooner
            # Note that groups are committed along the way.
            con.commit()
            con.close()

    if modified: # building index on particle time
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

        
# Where is the time going?
# A single group is 35s.
#         1    0.102    0.102   36.524   36.524 sql_post.py:206(add_ptm_group_to_db)
#      1681    0.054    0.000   12.533    0.007 postprocess_v00.py:600(eta_for_time)
#         2   10.677    5.338   10.677    5.338 {method 'executemany' of 'sqlite3.Cursor' objects}
#      1681    7.996    0.005    9.993    0.006 sql_post.py:319(<listcomp>)
#      8428    0.045    0.000    7.457    0.001 common.py:220(__getattr__)
#     76376    0.081    0.000    7.285    0.000 dataset.py:1225(__getitem__)
#     76376    0.780    0.000    7.147    0.000 dataset.py:1140(_construct_dataarray)
#      8427    0.009    0.000    6.848    0.001 dataset.py:1174(_attr_sources)
#      8427    0.022    0.000    6.837    0.001 dataset.py:1180(_item_sources)
#      8427    0.028    0.000    6.723    0.001 dataset.py:1187(<dictcomp>)
# 44490/10309    0.011    0.000    4.236    0.000 _asarray.py:16(asarray)
# 46023/11775    0.088    0.000    4.236    0.000 {built-in method numpy.array}

# After some memoizing, the initial scan is much faster.
# Building the point index is slow.
# Python mapping time 49.8s
#  -- that seems a lot slower. What's up?

# Maybe because precompute_cells is getting called?
#  eta_for_time is now just 10s.
