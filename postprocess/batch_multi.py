import six
import multiprocessing as mp
import numpy as np
import datetime
import os
import glob
import matplotlib.pyplot as plt
import logging as log
from stompy import utils, memoize, xr_utils
import xarray as xr
from matplotlib import colors
import re

from stompy.grid import unstructured_grid
import postprocess_v00 as post

out_dir="processed"
os.path.exists(out_dir) or os.makedirs(out_dir)

def process_batch_incr(ptm_runs,
                       time_range,
                       patterns,
                       z_ranges,
                       conc_fn,
                       max_age_days=15,
                       spinup=None,
                       version='v08'):
    """
    A memory-conscious version of process_batch. Moves more logic here, and
    avoids collecting such large datasets.

    v08: switch to new runs of ptm with full stormwater hydro inputs
         also update the loops below to be more memory-conscious
    """
    time_str=(utils.to_datetime(time_range[0]).strftime('%Y%m%d')
              + '_'
              + utils.to_datetime(time_range[1]).strftime('%Y%m%d'))

    max_age=np.timedelta64(max_age_days,'D')
    if spinup is None:
        spinup=max_age

    conc_func=post.conc_func_full(conc_fn)
        
    for group_name,group_patt in patterns:
        log.info(f"Processing {group_name}, pattern: {group_patt}")
        chunk_dir=os.path.join(out_dir,time_str,group_name)
        os.path.exists(chunk_dir) or os.makedirs(chunk_dir)

        # A dictionary to accumulate the final, gridded data.
        # limit to z_names that are missing
        conc_per_z_range={} # datasets
        conc_nc_fn_per_z_range={} # filename to save to
        
        for z_name,z_range in z_ranges:
            # Just the particles for the period, with mass, but not filtered
            # on elevation:
            conc_nc_fn=os.path.join(chunk_dir,f'particles-{z_name}-{max_age_days}days-{version}.nc')
            if os.path.exists(conc_nc_fn):
                log.info(f"netcdf {conc_nc_fn} exists")
            else:
                # existence of the key signals that it needs to be processed. the actual
                # dataset is created using the first group as a template, below.
                conc_per_z_range[z_name]=None
                conc_nc_fn_per_z_range[z_name]=conc_nc_fn
                
        if len(conc_per_z_range)==0:
            log.info(f"   no work for {group_name}, {group_patt}")
            continue
        
        #------
        # code that was in query_runs:

        assert group_patt.endswith('$')
        assert group_patt.startswith('^')

        # compile an array of particles that are
        #  (i) observed within the time_range, limited to output steps
        #      at least spinup into that specific run
        #  (ii) not older than max_age
        # and report each particle's
        #  (a) mass
        #  (b) age

        # tracers are assumed to be specific to a behavior, but
        # aggregated across runs and sources.
        # mass is normalized to unit concentration from each source

        # note that to normalize mass, we also need to know how many
        # particles were released in this group, including particles not
        # observed at this time.

        # fist step is to scan each matching group of each run, and
        # populate everything but mass
        all_part_obs=[]
        groups=[]    
        for run_idx,run in enumerate(ptm_runs):
            for group in run.groups():
                if re.match(group_patt,group) is None:
                    continue
                log.info(f"{run.run_dir:30s}: {group}")

                src_name,behavior_name=run.group_to_src_behavior(group)
                if src_name in post.skip_source:
                    log.info(f"Will skip source {src_name} -- it's in skip_source")
                    continue
                # here convert conc from particles/L to particles/m3.
                # this corresponds to the switch to v03 outputs.
                conc=1000*conc_func(group,src_name,behavior_name)
                if conc==0.0:
                    log.info(f"Will skip source {src_name}, behavior {behavior_name}, its concentration is 0")
                    continue
                groups.append(group)

                part_obs=post.scan_group(run,group,time_range=time_range,
                                         weight=conc,# drop run weights altogether
                                         max_age=max_age,spinup=spinup,
                                         grid=grid, extra_fields=[('run_idx',np.int32)])
                part_obs['run_idx']=run_idx

                # Convert to xarray Dataset for easier processing down the line.
                result=xr_utils.structure_to_dataset(part_obs,'particle',{'x':('xyz',)})

                # and include some metadata
                result['group_pattern']=(),group_patt
                result['time_start']=(),time_range[0]
                result['time_end']  =(),time_range[1]
                result['max_age']=(),max_age
                result['ptm_runs']=('ptm_run',), [p.run_dir for p in ptm_runs]
                result['ptm_groups']=('ptm_groups',), groups

                log.info("Adding vertical info")
                result=post.add_z_info(result,grid,ptm_runs)
                
                # who sets grp_rel_per_hour?
                assert np.isnan(result['grp_rel_per_hour']).sum()==0

                # these steps used to be in an outside loop
                # awaiting refactor.
                for z_name,z_range in z_ranges:
                    if z_name not in conc_per_z_range:
                        continue # already have this data, no need to do it again.
                    
                    # Just the particles for the period, with mass, but not filtered
                    # on elevation:
                    log.info(f"filtering by z")
                    result_z=post.filter_by_z_range(result,z_range,grid,ptm_runs)
                    conc=post.particle_to_density(result_z,grid,normalize='area')

                    if conc_per_z_range[z_name] is None:
                        # create the dataset
                        # this should preserve most of the metadata
                        ds=result_z.copy()
                        particle_vars=[v for v in ds.variables if 'particle' in ds[v].dims]
                        for v in particle_vars:
                            del ds[v]
                        grid.write_to_xarray(ds=ds)
                        ds['conc']=('face',),conc
                        ds['conc'].attrs['units']='particles m-2'
                        conc_per_z_range[z_name]=ds
                    else:
                        ds=conc_per_z_range[z_name]
                        # just increment
                        ds['conc'].values += conc
                            
        # after loop
        for z_name,z_range in z_ranges:
            if z_name not in conc_per_z_range:
                continue # already have this data, no need to do it again.
            conc_nc_fn=conc_nc_fn_per_z_range[z_name]
            ds=conc_per_z_range[z_name]
            ds.to_netcdf(conc_nc_fn)
            log.info(f"done writing netcdf to {conc_nc_fn}")

def process_batch(ptm_runs,
                  time_range,
                  patterns,
                  z_ranges,
                  conc_fn,
                  max_age_days=15,
                  spinup=None,
                  version='v08'):
    # v03: shift to particles/m3 units (in sync with change in postprocess_v00).
    #   to speed things up, simply adjust v02 output if that exists.
    # v04: should be same as v03, but is recomputed, not just post-hoc scaled.
    # v05: start scaling up stormwater, too.
    # v06: include new concentration inputs that derate stormwater and wastewater
    #       per-category for blank contamination
    # v07: I think this is what was used for the report and symposium. have to check
    #      git logs to see the change from v06. probably something related to blanks.
    # v08: switch to new runs of ptm with full stormwater hydro inputs
    #      also update the loops below to be more memory-conscious

    # spinup: timedelta such that when processing each run, only record particles
    #   beyond spinup into the run.

    # Note that grid is global, and used in here.  that's because these calls
    # are made via multiprocessing, and I don't want to deal with the potential
    # headache of passing a grid object through multiprocessing.
    time_str=(utils.to_datetime(time_range[0]).strftime('%Y%m%d')
              + '_'
              + utils.to_datetime(time_range[1]).strftime('%Y%m%d'))
    
    for group_name,group_patt in patterns:
        log.info(f"Processing {group_name}, pattern: {group_patt}")
        chunk_dir=os.path.join(out_dir,time_str,group_name)
        os.path.exists(chunk_dir) or os.makedirs(chunk_dir)
            
        if 1: # old way, takes too much memory:
            # calculated on-demand in the loop below
            @memoize.memoize()
            def parts():
                log.info("Extracting particles")
                result=post.query_runs(ptm_runs,
                                       group_patt=group_patt,
                                       time_range=time_range,
                                       max_age=np.timedelta64(max_age_days,'D'),
                                       spinup=spinup,
                                       conc_func=post.conc_func_full(conc_fn),
                                       grid=grid)
                if result is not None:
                    log.info("Adding vertical info")
                    result=post.add_z_info(result,grid,ptm_runs)
                else:
                    log.info("No data for this set of runs")
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
                if os.path.exists(conc_nc_fn):
                    log.info(f"netcdf {conc_nc_fn} exists")
                else:
                    log.info(f"preparing data for {conc_nc_fn}")
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
                    log.info(f"done writing netcdf to {conc_nc_fn}")

def process_batch_onearg( full ):
    a,k=full
    return process_batch_incr(*a,**k)


# load this globally so each process gets a copy, rather than trying to
# pass it through the multiprocessing call
grid_fn="/opt2/sfb_ocean/suntans/runs/merged_020_20170610/ptm_average.nc_0000.nc"  
grid=post.grid_from_ptm_hydro(grid_fn)

if __name__=="__main__":                
    patterns=[
        ('-0.05','^.*_down50000$'),
        ('-0.005','^.*_down5000$'),
        ('-0.0005','^.*_down500$'),
        ('0.0','^.*_none$'),
        ('0.0005','^.*_up500$'),
        ('0.005','^.*_up5000$'),
        ('0.05','^.*_up50000$')
    ]

    # previously these went to 0.5, consistent with the behaviors.
    # new runs use 0.2
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

    conc_fns=[
        ('std',"../loads/plastic_loads-7classes-v03.nc"),
        ('nofiber',"../loads/plastic_loads-7classes-v03-nofiber.nc")
    ]

    # previous, as used in report and symposium ppt
    # ptm_base_dir="/opt2/sfb_ocean/ptm/all_source"
    # updated with new hydro, passive regions.
    ptm_base_dir="/opt2/sfb_ocean/ptm/all_source_020"
    
    log.info("Gathering list of runs")
    calls=[] # (*a,**kw) for calls to process_batch

    # v07: target for the report, symposium.
    # v08: updated with local hydro w/ full stormwater, and slimmer passive regions.
    base_version='v08'
    
    for conc_version,conc_fn in conc_fns:
        if 1:
            version=base_version+conc_version # e.g. v08nofiber
            
            # the basic, 15 day setup
            for run_date in run_dates:
                ptm_runs=[post.PtmRun(run_dir=d) 
                          for d in glob.glob(os.path.join(ptm_base_dir,f"{run_date}/w*")) ]
                if len(ptm_runs)!=7:
                    log.warning("Incomplete set of runs: ")
                    for r in ptm_runs:
                        log.warning(f"  {r}")
                    continue
                # PTM is taking a long time, so be kind to incomplete runs.
                #assert len(ptm_runs)==7

                # just the time period with a full field for max_age=15D
                start=np.timedelta64(15,'D') + utils.to_dt64(datetime.datetime.strptime(run_date,'%Y%m%d'))

                time_range=[start,start+np.timedelta64(15,'D')]

                for pattern in patterns:
                    # one pattern at a time
                    calls.append( ([ptm_runs,time_range,[pattern]],
                                   dict(z_ranges=z_ranges,
                                        conc_fn=conc_fn,
                                        version=version)) )
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

                for AB_version,run in [ ('A',runA),
                                        ('B',runB) ]:
                    version=base_version+AB_version+conc_version # e.g. v07Anofiber
                    # need the 'w' prefix to skip bak directories
                    run_dirs=glob.glob(os.path.join(ptm_base_dir,f"{run}/w*"))
                    if len(run_dirs)!=7:
                        log.warning(f"Expected 7 runs, got")
                        for d in run_dirs:
                            log.warning(" " + d)
                        continue
                        # raise Exception("Unexpected number of runs")
                    ptm_runs=[post.PtmRun(run_dir=d) 
                              for d in run_dirs ]

                    for pattern in patterns:
                        # one pattern at a time
                        # will these get named okay? output goes into a folder with the
                        # start and end dates, but have to include v06A / v06B to distinguish
                        # the two parts, because the outputs do not include the run_dir.
                        # spinup has to be set to zero, otherwise we'd drop all of runB when
                        # it defaults to max age.
                        calls.append( ([ptm_runs,time_range,[pattern]],dict(max_age_days=max_age_days,
                                                                            spinup=np.timedelta64(0,'D'),
                                                                            version=version,
                                                                            conc_fn=conc_fn,
                                                                            z_ranges=z_ranges)) )
    
    log.info("%d invocations"%len(calls))
    
    with mp.Pool(8) as pool:
        count=0
        for ret in pool.imap_unordered(process_batch_onearg, calls):
            count+=1
            log.info(f"--------{count} / {len(calls)} --------" )
