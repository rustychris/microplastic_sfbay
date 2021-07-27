"""
Postprocessing code to query PTM runs via multiprocessing
"""

import logging
import os, glob
import pandas as pd
import xarray as xr
import re
import numpy as np

#import dask
#import dask.dataframe as dd
import multiprocessing as mp

from stompy.model.fish_ptm import ptm_tools
from stompy.model.suntans import sun_driver
from stompy import utils
from stompy.grid import unstructured_grid

log=logging.getLogger('postproc')

from scipy import sparse

from stompy.spatial import proj_utils
ll2utm=proj_utils.mapper('WGS84','EPSG:26910')

# Is it any faster to just use basic multiprocessing to run queries?
# dask is proving to be tricky to debug and not particularly fast.
import postproc_dask as post

from postproc_dask import (config_malloc,get_load_data,
                           w_s_map,source_ptm_to_load,
                           sun_path_to_hydro, hydro_ds,
                           godwin,godwin_offset_h,
                           load_hydro_timestamps, extract_bc_ds,
                           bc_ds, source_Qdata, run_to_group_paths,
                           parse_group_path, criteria_to_groups,
                           load_conc,volume_per_particle,
                           query_group_particles, t_to_idx,
                           get_particle_attrs, filter_particles_post_attrs,
                           particles_to_conc, rec_to_cell_weights, 
                           group_weights )



from cfg_v01 import cfg

# Load the grid into... grid
hydro_path=cfg['sun_paths'][0]
ptm_ds=xr.open_dataset(os.path.join(hydro_path,"ptm_average.nc_0000.nc"))
grid=unstructured_grid.UnstructuredGrid.read_ugrid(ptm_ds,dialect='fishptm')
ptm_ds.close()

cfg['bc_ds_d']=bc_ds(cfg)
cfg['load_data_d']=get_load_data()

def helper(args):
    grp_fn,criteria=args
    particles=query_group_particles(grp_fn,criteria,cfg['load_data_d'],cfg['bc_ds_d'],cfg['info_version'])

    assert 'cell' in particles.columns
    cell=particles['cell'].values
    particles['z_bed']=grid.cells['z_bed'][cell].astype(np.float32)
    assert 'z_surface' in particles.columns
    age_s=(particles['time']-particles['rel_time'])/np.timedelta64(1,'s')
    particles['age_s']=age_s.astype(np.float32)

    # Filter at this level
    part_filtered=filter_particles_post_attrs(particles,criteria)
    
    return part_filtered

def query_particles(criteria,cfg,compute_kw={}):
    """
    Takes a criteria dictionary and returns an
    uncomputed dask dataframe with particles satisfying
    the criteria
    Includes particle counts scaled from loads, relative
    z coordinates.
    """
    # Warm start more like 2s. 
    groups=criteria_to_groups(criteria,cfg=cfg)
    #groups=groups[:50] # DBG
    
    if 'pool' in compute_kw:
        pool=compute_kw['pool']
        group_iter=pool.imap_unordered(helper,
                                       [ [grp_fn,criteria] for grp_fn in groups])
    else:
        group_iter=(helper(grp_fn,criteria) for grp_fn in groups)

    mem_total=0
    group_data=[]
    for particles in group_iter:
        mem_total+=particles.memory_usage().sum()
        group_data.append(particles)
        
    print("Total memory before concat: ",mem_total)
    parts=pd.concat(group_data)
    
    return parts

# This doesn't need to be specialized for mp, but we need to call the mp query_particles.
def query_particle_concentration(criteria,cfg,grid,decay=None,compute_kw={}):
    """
    criteria: dictionary with selection criteria
    cfg: dictionary with experiment configuration data
    decay: np.timedelta64 giving a decay e-folding timescale

    Map particles, with optional decay function, onto grid 
    as concentrations. Returns xr.Dataset suitable for BayConcFigure.
    """
    # First, take a sequential approach, time it, then come back
    # to push more postprocessing out to nodes.
    
    # For long-ish time periods this runs into serious memory issues.
    print("Query particles")
    particles=query_particles(criteria=criteria,cfg=cfg)
    print(f"Compute done. {len(particles)} particles")
 
    age_s=particles['age_s'].values
    age_max=criteria.get('age_max', np.timedelta64(60,'D'))
    age_max_s=age_max/np.timedelta64(1,'s')
    weights=np.where( age_s<age_max_s,1.0,0.0)
    
    if decay is not None:
        decay_s=decay/np.timedelta64(1,'s')
        weights*=np.exp(-age_s/tau_s)
        
    particles['weighted_count']=weights*particles['mp_per_particle'].values
    
    output_steps=1+int( (criteria['t_max'] - criteria['t_min'])/self.cfg['ptm_output_interval'])
    scale=1./output_steps
    
    ds=particles_to_conc(particles,grid,smooth=0,scale=scale,
                         count_field='weighted_count')

    ds['conc'].attrs['grp_filter']=criteria.get('behavior','')
    ds['criteria']=(),criteria
    
    # Manual determination of name of z_filter
    z_filters=[]
    if 'z_below_surface_max' in criteria:
        z_filters.append("surface %.3f"%criteria['z_below_surface_max'])
    if 'z_above_bed_max' in criteria:
        z_filters.append("bed %.3f"%criteria['z_above_bed_max'])
    
    ds['conc'].attrs['z_filter']=" and ".join(z_filters)
    ds['conc'].attrs['t_start']=criteria['t_min']
    ds['conc'].attrs['t_stop']=criteria['t_max']
    ds['conc'].attrs['max_age']=age_max
    ds['conc'].attrs['scale']=scale
    ds['conc'].attrs['output_steps']=output_steps
    
    return ds

def particles_for_date(rec_DATE,cfg,criteria={},include_godin=True,cache=True,
                       compute_kw={}):
    """
    Query particles centered around the given date, as a string.
    Assumes the date has no time information

    The time window is large enough to allow a godin tidal filter
    
    cfg['manta_out_dir'] controls where results are cached.

    criteria will be added to defaults below
    """
    t_sample=np.datetime64(rec_DATE)

    if cache and len(criteria):
        raise Exception("Caching only available for default criteria")
    
    if cache:
        out_dir=cfg['manta_out_dir']
        fn=os.path.join(out_dir,f"v01-{rec_DATE[:10]}.nc")
    
    if (not cache) or (not os.path.exists(fn)):
        # pull a generous buffer of particles here, and narrow
        # the time zones are annoying but I double-checked and
        # this does give enough of a buffer.
        
        # want to be able to, after the fact, query a full 25h tidal cycle
        # centered on the actual time of a sample that could fall anywhere
        # in this day.
        t_start=t_sample+np.timedelta64(8,'h') - np.timedelta64(12,'h')
        # and go for a tidal day
        t_stop =t_sample+np.timedelta64(8+24,'h') + np.timedelta64(13,'h')
        query_n_steps=25 # how many hours are included in the query.

        # 2020305 changes:
        t_start-=np.timedelta64(24,'h')
        t_stop+=np.timedelta64(24,'h')
        # This is just enough of a time window to calculate a godin
        # filter.
        query_n_steps += 48

        # Something like this:
        defaults=dict(t_min=np.datetime64(t_start), 
                      t_max=np.datetime64(t_stop), 
                      category='nonfiber',
                      z_below_surface_max=0.095,
                      age_max=np.timedelta64(60,'D'))
        defaults.update(criteria)
        df=query_particles(criteria=defaults,cfg=cfg,compute_kw=compute_kw)
        if cache:
            df.to_parquet(fn)
    else:
        df=pd.read_parquet(fn)

    if include_godin:
        # Add godin weights:
        # Noon, local, day of sampling
        t_center = t_sample+np.timedelta64(8,'h') + np.timedelta64(12,'h')
        delta_hours=((df['time']-t_center)/np.timedelta64(1,'h')).astype(np.int32)
        df['weight_time']=godwin[delta_hours+godwin_offset_h]
        
    return df

import resource

if __name__=="__main__":
    if not os.path.exists(cfg['manta_out_dir']):
        os.makedirs(cfg['manta_out_dir'])

    ru=resource.getrusage(resource.RUSAGE_SELF)
    print("RSS: ",ru.ru_maxrss)
        
    pool=mp.Pool(24)
    import time
    t=time.time()
    # import cProfile as prof
    df=particles_for_date("2017-10-21 00:00:00",cfg=cfg,cache=False,
                          compute_kw={'pool':pool})
    print("Rows returned: ",len(df))
    print("Elapsed seconds:",time.time() - t)
    print("Memory for final result:")
    print(df.memory_usage().sum())
    
    ru=resource.getrusage(resource.RUSAGE_SELF)
    print("RSS: ",ru.ru_maxrss)

# So this runs, maybe it takes 15 minutes to run the queries.
# but whether here or in dask, need to filter on manta sample
# earlier.

# with 50 groups, the pre-concat memory usage is 673_179_784
# So the troubling part is that the memory usage before concatenation
# is so large.  over 10x larger than it should be.
# and the post-concat memory usage 53_513_068
# max_rss: 1_778_096k, 1.8GB
# Changing to applying the postprocessing step on the fly on the iterator
# return values. No change in memory. runs the same.
# Change the data types to be more efficient. Drops 15-20% off the final result size.
# Does copying the individual dataframes slim up the aggregate memory usage?
# nope.
# Okay - it was because extra filtering was being saved for after the
# 

# So the full query now runs. Result size:
# 4_069_170_308
# Good. Took 1046s, using 24 processes.
# Main process 8_481_684k max RSS, so double the data size.  Seems reasonable.
# 
