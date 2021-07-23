"""
Running totals of number of bed and shoreline interactions
"""
# Map all particles to bed hits, shoreline hits.
# [timestep,id] => (n_bed,n_shore)
#

# HERE

from stompy.grid import unstructured_grid
import postproc_dask as post

import xarray as xr
import numpy as np
import os, glob

cfg=dict(
    ptm_base_dir="/opt2/sfb_ocean/ptm/all_source_022b",
    sun_base_dir="/opt2/sfb_ocean/suntans/runs",
    ptm_output_interval=np.timedelta64(1,'h')
)

cfg['ptm_run_patt']=os.path.join(cfg['ptm_base_dir'],"chunk??","20??????")
cfg['sun_patt']=os.path.join(cfg['sun_base_dir'],"merged_022_20??????")

ptm_run_paths=glob.glob(cfg['ptm_run_patt'])
ptm_run_paths.sort()
cfg['ptm_run_paths']=ptm_run_paths

sun_paths=glob.glob(cfg['sun_patt'])
sun_paths.sort()
cfg['sun_paths']=sun_paths

cfg['hydro_timestamps']=post.load_hydro_timestamps(sun_paths)

# Load the grid into... grid
hydro_path=sun_paths[0]
ptm_ds=xr.open_dataset(os.path.join(hydro_path,"ptm_average.nc_0000.nc"))
grid=unstructured_grid.UnstructuredGrid.read_ugrid(ptm_ds,dialect='fishptm')

ptm_ds.close()

tri,tsrcs=grid.mpl_triangulation(return_sources=True)
tf=tri.get_trifinder()

## 
from stompy.model.fish_ptm import ptm_tools

# Best to just specify
Nids=2400

def process_group(group):
    nc_fn=group+"-v00.nc"
    if os.path.exists(nc_fn):
        return
    print("Processing %s"%group)
    
    pbf=ptm_tools.PtmBin(group)
    times=pbf.time.astype('<M8[us]')

    # should be on the order of a 30M nc
    ds=xr.Dataset()
    ds['time']=('time',),times

    _,particles=pbf.read_timestep(len(times)-1)

    # id will be stored 0-based, even though particles['id'] is 1-based
    # Nids=particles['id'].max() # -1+1
    ds['id']=('id',),np.arange(Nids)

    tp_cell=np.zeros( (ds.dims['time'],ds.dims['id']), np.int32)
    tp_cell[...]=-1

    for ti,t in enumerate(times):
        datenum,particles=pbf.read_timestep(ti)
        points=particles['x'][:,:2]
        tcells=tf(points[:,0],points[:,1])
        cells=np.where(tcells>=0,
                       tsrcs[tcells],
                       -1)
        tp_cell[ti,particles['id']-1]=cells

    ds['cell']=('time','id'),tp_cell

    all_t,all_cell=xr.broadcast(ds['time'],ds['cell'])

    # This is a bit slower....  maybe 5s?
    z_surf=post.get_z_surface(all_t.values.ravel(),
                              all_cell.values.ravel(),
                              cfg)

    ds['z_surf']=('time','id'),z_surf.reshape([ds.dims['time'],ds.dims['id']]).astype(np.float32)

    # 32M of data, compresses to 14M
    encoding=dict(cell=dict(zlib=True,complevel=5),z_surf=dict(zlib=True,complevel=5))
    ds.to_netcdf(nc_fn,encoding=encoding)
    
# There are 152 runs
# Each run has about 90 groups
# So that's something like 13k groups?
# so 7h single core.
# with cell and z_surface, this is 5.0s (with hot data)
# 7s with cold data

if __name__ == '__main__':
    groups=[g for r in cfg['ptm_run_paths'] for g in post.run_to_group_paths(r)]
    print("%d groups"%len(groups))

    import multiprocessing

    with multiprocessing.Pool(16) as pool:
        for res in pool.imap_unordered(process_group,groups,chunksize=5):
            print(".",end="")

