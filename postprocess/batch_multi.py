import six
import multiprocessing as mp
import numpy as np
import datetime
import os
import glob
import matplotlib.pyplot as plt
import logging as log
from stompy import utils
from stompy import memoize
import xarray as xr
from matplotlib import colors

from stompy.grid import unstructured_grid
import postprocess_v00 as post

grid_fn="/opt2/sfb_ocean/suntans/runs/merged_018_20171227/ptm_average.nc_0000.nc"  
grid=post.grid_from_ptm_hydro(grid_fn)

out_dir="processed"
os.path.exists(out_dir) or os.makedirs(out_dir)
    
def process_batch(ptm_runs,
                  time_range,
                  patterns,
                  z_ranges,
                  max_age_days=15,
                  spinup=None,
                  version='v05'):
    # v03: shift to particles/m3 units (in sync with change in postprocess_v00).
    #   to speed things up, simply adjust v02 output if that exists.
    # v04: should be same as v03, but is recomputed, not just post-hoc scaled.
    # v05: start scaling up stormwater, too.

    # spinup: timedelta such that when processing each run, only record particles
    #   beyond spinup into the run.
    time_str=(utils.to_datetime(time_range[0]).strftime('%Y%m%d')
              + '_'
              + utils.to_datetime(time_range[1]).strftime('%Y%m%d'))
    
    for group_name,group_patt in patterns:
        log.info(f"Processing {group_name}, pattern: {group_patt}")
        chunk_dir=os.path.join(out_dir,time_str,group_name)
        os.path.exists(chunk_dir) or os.makedirs(chunk_dir)
  
        # calculated on-demand in the loop below
        @memoize.memoize()
        def parts():
            log.info("Extracting particles")
            result=post.query_runs(ptm_runs,
                                   group_patt=group_patt,
                                   time_range=time_range,
                                   z_range=None, # not ready
                                   max_age=np.timedelta64(max_age_days,'D'),
                                   spinup=spinup,
                                   conc_func=post.conc_func,
                                   grid=grid)
            log.info("Adding vertical info")
            result=post.add_z_info(result,grid,ptm_runs)
            return result
        for z_name,z_range in z_ranges: 
            # Just the particles for the period, with mass, but not filtered
            # on elevation:
            conc_nc_fn=os.path.join(chunk_dir,f'particles-{z_name}-{max_age_days}days-{version}.nc')
            #old_conc_nc_fn=os.path.join(chunk_dir,f'particles-{z_name}-{max_age_days}days-v02.nc')
            #if not os.path.exists(conc_nc_fn) and os.path.exists(old_conc_nc_fn):
            #    ds=xr.open_dataset(old_conc_nc_fn)
            #    ds.conc.values *= 1000 # scale up to account for old code that used particles/L
            #    ds.conc.attrs['units']='particles m-2'
            #    ds.to_netcdf(conc_nc_fn)
            #    ds.close()
            #    log.info("Rescaled v02 output to v03")
            if not os.path.exists(conc_nc_fn):
                log.info(f"writing to {conc_nc_fn}")
                p=parts()
                p=post.filter_by_z_range(p,z_range,grid,ptm_runs)
                conc=post.particle_to_density(p,grid,normalize='area')
                # could also use the z_bed, z_surf values to turn particle mass into
                # a mass/z, such that normalize by area then gives a volume concentration.
                # unless it's full water column, have to do some truncating
                
                # this should preserve most of the metadata
                ds=p.copy()
                particle_vars=[v for v in ds.variables if 'particle' in ds[v].dims]
                for v in particle_vars:
                    del ds[v]
                
                grid.write_to_xarray(ds=ds)
                ds['conc']=('face',),conc
                ds['conc'].attrs['units']='particles m-2'
                
                ds.to_netcdf(conc_nc_fn)
                log.info("done writing")

def process_batch_onearg( full ):
    a,k=full
    return process_batch(*a,**k)
        
if __name__=="__main__":                
    patterns=[
        ('-0.05','.*_down50000'),
        ('-0.005','.*_down5000'),
        ('-0.0005','.*_down500'),
        ('0.0','.*_none'),
        ('0.0005','.*_up500'),
        ('0.005','.*_up5000'),
        ('0.05','.*_up50000')
    ]

    z_ranges=[
        ('bed',[0,0.5]),
        ('surf',[-0.5,0]),
        ('avg',[0,0])
    ]

    
    # run_date: YYYYMMDD string for start of runs
    run_dates=["20170615",
               "20170715",
               "20170815",
               "20170915",
               "20171015",
               "20171115",
               "20171215",
               "20180115",
               "20180215",
               "20180315",
               "20180415",
               "20180515"
    ]
    
    log.info("Gathering list of runs")
    calls=[] # (*a,**kw) for calls to process_batch

    if 1:
        # the basic, 15 day setup
        for run_date in run_dates:
            ptm_runs=[post.PtmRun(run_dir=d) 
                      for d in glob.glob(f"/opt2/sfb_ocean/ptm/all_source/{run_date}/w*") ]
            assert len(ptm_runs)==7

            # just the time period with a full field for max_age=15D
            start=np.timedelta64(15,'D') + utils.to_dt64(datetime.datetime.strptime(run_date,'%Y%m%d'))

            time_range=[start,start+np.timedelta64(15,'D')]
            # process_batch(ptm_runs,time_range,patterns,z_ranges=z_ranges)
            for pattern in patterns:
                # one pattern at a time
                calls.append( ([ptm_runs,time_range,[pattern]],dict(z_ranges=z_ranges)) )
    if 1:
        # quirky 45-60 day setup.
        # any particular run starts on the 15th of the month.
        # runs for 60 days.
        # e.g. 2017-09-15 00:00 to 2017-11-14 00:00
        # exactly 60 days.
        #

        # in order to have [almost] a full field, to average over a spring-neap,
        # and to maximize max age, the idea is to combine two runs at a time,
        # and average over days 44-58 of the first run.
        for runA,runB in zip(run_dates[:-1],run_dates[1:]):
            log.info(f"considering a joint analysis of {runA}, {runB}")

            max_age_days=44 # because of differences in the lengths of months, this can't be
            # 45.
            runA_start=utils.to_dt64(datetime.datetime.strptime(runA,'%Y%m%d'))
            start=runA_start + np.timedelta64(44,'D')
            stop =runA_start + np.timedelta64(58,'D')
            time_range=[start,stop]

            for ver,run in [ ('v05A',runA),
                             ('v05B',runB) ]:
                # need the 'w' prefix to skip bak directories
                run_dirs=glob.glob(f"/opt2/sfb_ocean/ptm/all_source/{run}/w*")
                if len(run_dirs)!=7:
                    log.error(f"Expected 7 runs, got")
                    for d in run_dirs:
                        log.error(" " + d)
                    raise Exception("Unexpected number of runs")
                ptm_runs=[post.PtmRun(run_dir=d) 
                          for d in run_dirs ]

                for pattern in patterns:
                    # one pattern at a time
                    # will these get named okay? output goes into a folder with the
                    # start and end dates, but have to include v05A / v05B to distinguish
                    # the two parts, because the outputs do not include the run_dir.
                    # spinup has to be set to zero, otherwise we'd drop all of runB when
                    # it defaults to max age.
                    calls.append( ([ptm_runs,time_range,[pattern]],dict(max_age_days=max_age_days,
                                                                        spinup=np.timedelta64(0,'D'),
                                                                        version=ver,
                                                                        z_ranges=z_ranges)) )
    
    log.info("%d invocations"%len(calls))
    
    with mp.Pool(16) as pool:
        pool.map(process_batch_onearg, calls)
