import numpy as np
import os
import xarray as xr
import glob
from stompy import memoize
from stompy.grid import unstructured_grid
import postproc_dask

# Ideally this would be split into a ConfigBase with the logic
# and a Config subclass with specifics

class Config:
    ptm_base_dir="/opt2/sfb_ocean/ptm/all_source_022b"
    sun_base_dir="/opt2/sfb_ocean/suntans/runs"
    ptm_output_interval=np.timedelta64(1,'h')
    manta_out_dir="manta_sets_20211001"
    info_version='v01'
    load_version='v06_355'

    def __getstate__(self):
        s=dict(self.__dict__)
        s.pop('_memocache',-1)
        return s
                
    # Derived properties
    @property
    @memoize.imemoize()
    def ptm_run_patt(self):
        return os.path.join(self.ptm_base_dir,"chunk??","20??????")
    @property
    @memoize.imemoize()
    def sun_patt(self):
        return os.path.join(self.sun_base_dir,"merged_022_20??????")
    @property
    @memoize.imemoize()
    def ptm_run_paths(self):
        run_paths=glob.glob(self.ptm_run_patt)
        run_paths.sort()
        return run_paths
    @property
    @memoize.imemoize()
    def sun_paths(self):
        sp=glob.glob(self.sun_patt)
        sp.sort()
        return sp

    @property
    @memoize.imemoize()
    def bc_ds_d(self):
        return postproc_dask.bc_ds(self)

    @property
    @memoize.imemoize()
    def load_data_d(self):
        return postproc_dask.get_load_data(version=self['load_version'])

    # look like a dict, too.
    def __getitem__(self,k):
        return getattr(self,k)

    # Load the grid into... grid
    @property
    @memoize.imemoize()
    def grid(self):
        hydro_path=self.sun_paths[0]
        ptm_ds=xr.open_dataset(os.path.join(hydro_path,"ptm_average.nc_0000.nc"))
        grid=unstructured_grid.UnstructuredGrid.read_ugrid(ptm_ds,dialect='fishptm')
        ptm_ds.close()
        return grid
    @property
    def grid_d(self):
        return self.grid

    @property
    @memoize.imemoize()
    def hydro_timestamps(self):
        return postproc_dask.load_hydro_timestamps(self.sun_paths)    
    
cfg=Config()
