# Map all particles to grid cells, saving, per group, a mapping from
# [timestep,id] => (cell,z_bed,z_surf)
#

from stompy.grid import unstructured_grid
import postproc_dask as post

import xarray as xr
import numpy as np
import os, glob

import cfg_v01 
cfg=dict(cfg_v01.cfg)
cfg['hydro_timestamps']=post.load_hydro_timestamps(cfg['sun_paths'])

# Load the grid into... grid
hydro_path=cfg['sun_paths'][0]
ptm_ds=xr.open_dataset(os.path.join(hydro_path,"ptm_average.nc_0000.nc"))
grid=unstructured_grid.UnstructuredGrid.read_ugrid(ptm_ds,dialect='fishptm')

ptm_ds.close()

tri,tsrcs=grid.mpl_triangulation(return_sources=True)
tf=tri.get_trifinder()

CS=post.CellStatus(cfg=cfg,grid=grid)

## 
from stompy.model.fish_ptm import ptm_tools

# Best to just specify
Nids=2400

def process_group(group):
    nc_fn=group+"-v01.nc"
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

    # bed and shoreline encounters:
    # Count up the encounters
    d_thresh=0.25 # arbitrary, but at least this makes it "consistent" with dry cells.
    tp_hit_bed=np.zeros( (ds.dims['time'],ds.dims['id']), np.int16)
    tp_hit_shore=np.zeros( (ds.dims['time'],ds.dims['id']), np.int16)
                           
    for ti,t in enumerate(times):
        datenum,particles=pbf.read_timestep(ti)
        points=particles['x'][:,:2]
        tcells=tf(points[:,0],points[:,1])
        valid=tcells>=0
        cells=np.where(valid, tsrcs[tcells], -1)
        pid=particles['id']-1
        tp_cell[ti,pid]=cells

        # Bed hits:
        z_part=particles['x'][:,2]
        # 999.0: if no cell was found assume it's dry.
        z_bed=np.where(valid, grid.cells['z_bed'][cells], 999.0 )
        hit_bed= (z_part-z_bed < d_thresh)
        tp_hit_bed[ti,pid]=1*hit_bed
        
        # And check whether we're adjacent to a boundary or 
        # dry cell
        # Pretty slow. 0.2 seconds per call, most of which is in cell_to_cells.
        # Got that down to 0.01s per call.
        # This could maybe be faster by integrating with get_z_surface, but in get_z_surface
        # we query specific cells, but for CS have to check all cells for adjacent dry cells
        cell_stats=CS(t) # 1: dry 2: dry or adjacent to dry, 4: boundary
        hit_shore=np.where(valid, cell_stats[cells]>0, True)
        tp_hit_shore[ti,pid]=1*hit_shore

    ds['bed_hits']=('time','id'), np.cumsum( tp_hit_bed, axis=0,dtype=np.int16)
    ds['shore_hits']=('time','id'), np.cumsum( tp_hit_shore, axis=0,dtype=np.int16)
    ds['cell']=('time','id'),tp_cell

    all_t,all_cell=xr.broadcast(ds['time'],ds['cell'])

    # This is a bit slower....  maybe 5s?
    z_surf=post.get_z_surface(all_t.values.ravel(),
                              all_cell.values.ravel(),
                              cfg)

    ds['z_surf']=('time','id'),z_surf.reshape([ds.dims['time'],ds.dims['id']]).astype(np.float32)
    
    # 32M of data, compresses to 14M
    comp=dict(zlib=True,complevel=5)
    encoding=dict(cell=comp,z_surf=comp,bed_hits=comp,shore_hits=comp)
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

    if 1: 
        import multiprocessing

        with multiprocessing.Pool(16) as pool:
            for res in pool.imap_unordered(process_group,groups,chunksize=5):
                print(".",end="")
    else:
        process_group(groups[0])
