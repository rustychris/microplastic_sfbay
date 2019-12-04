"""
Post-processing for the earlier all_source PTM runs
"""
from sql_post import *

#run=post.PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source_020/20170715/w0.0")
#fn=os.path.join(run.run_dir,"ptm_and_grid.db")
# getting ready for the older runs.

months=[
    #"/opt2/sfb_ocean/ptm/all_source/20170715",
    "/opt2/sfb_ocean/ptm/all_source/20170815",
    "/opt2/sfb_ocean/ptm/all_source/20170915",
    "/opt2/sfb_ocean/ptm/all_source/20171015",
    "/opt2/sfb_ocean/ptm/all_source/20171115",
    "/opt2/sfb_ocean/ptm/all_source/20171215",
    "/opt2/sfb_ocean/ptm/all_source/20180115",
    "/opt2/sfb_ocean/ptm/all_source/20180215",
    "/opt2/sfb_ocean/ptm/all_source/20180315",
    "/opt2/sfb_ocean/ptm/all_source/20180415",
    "/opt2/sfb_ocean/ptm/all_source/20180515"
]


for month in months:
    # Try shoving all of one month into the same database
    speeds=[
        "w0.0",
        "w-0.0005",
        "w0.0005",
        "w-0.005",
        "w0.005",
        "w-0.05",
        "w0.05"
    ]
    
    run_dirs=[os.path.join(month,speed)
              for speed in speeds ]

    fn=os.path.join(month,"ptm_and_grid.db")

    if 0: # building the original databases
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
