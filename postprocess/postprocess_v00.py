"""
Postprocess PTM output with constant particles/release, normalize particle weight
by flow.

The goal is to generate concentration fields on the computational grid reflecting
each source x settling velocity combination.

"""
import glob
import numpy as np
import xarray as xr
from stompy import utils, memoize, xr_utils
from stompy.model.fish_ptm import ptm_tools
from stompy.model.suntans import sun_driver
from stompy.grid import unstructured_grid

import os
import logging as log
log.basicConfig(level=log.INFO)

log.root.setLevel(log.INFO)

import re
import matplotlib.pyplot as plt

##

# Extract the relevant flow data from the BC files.
class PtmRun(object):
    run_dir=None
    # in case runs are copied from one place to another
    run_dir_mapping=[('/shared2/src/sfb_ocean/','/opt2/sfb_ocean/')]
    
    def __init__(self,**kw):
        utils.set_keywords(self,kw)

    def hydrodynamics_inp(self):
        return self.parse_sql(os.path.join(self.run_dir,'FISH_PTM_hydrodynamics.inp'))

    def ptm_hydro_files(self):
        """ List of the ptm_average*.nc files that are used by this run
        """
        path=None
        files=[]
        for tok in self.hydrodynamics_inp():
            if tok[0]=='HYDRO_FILE_PATH':
                path=self.remap_path(tok[1])
            elif tok[0]=='FILENAME':
                files.append( os.path.join(path,tok[1]) )
                path=None
        return files

    def remap_path(self,path):
        for src,tgt in self.run_dir_mapping:
            # might get fancier if the paths were weirder, but this
            # should be good enough here.
            if path.startswith(src):
                log.info(f"Remapping {src} => {tgt} in {path}")
                path=path.replace(src,tgt)
        return path
    
    @memoize.imemoize()
    def hydro_models(self):
        hydro=self.hydrodynamics_inp()
        all_paths=[tok[1] for tok in hydro if tok[0]=='HYDRO_FILE_PATH']
        paths=[]
        for p in all_paths:
            if paths and paths[-1]==p: continue
            paths.append(p)
        # Actually it's useful to load the model files to get the true model duration
        # as the BC files have extra fake data.
        models=[]
        for p in paths:
            p=self.remap_path(p)
            # loading the grid can slow it down, and really only need the
            # bc dataset.
            model=sun_driver.SuntansModel.load(p,load_grid=False)
            if model is None:
                log.warning("Could not load model from path %s"%p)
            else:
                models.append(model)
        return models

    def bc_ds(self):
        """ Extract the relevant parts of the BC data, return as a single dataset
        """
        compiled_fn=os.path.join(self.run_dir,'bc_extracted_v3.nc')
        
        if not os.path.exists(compiled_fn):
            dss=[]

            for model in self.hydro_models():
                # only care about point sources, and river inflows (i.e. ignore
                # ocean flux BCs, and any freesurface BCs
                model.load_bc_ds()
                ds=model.bc_ds.copy()
                ti_start,ti_stop=np.searchsorted(ds.time.values,[model.run_start,model.run_stop])
                ds=ds.isel(Nt=slice(ti_start,ti_stop+1))

                for extra in ['T','S','h','boundary_h','boundary_w','boundary_T',
                              'boundary_u','boundary_v','z',
                              'point_S','point_T','cellp','xv','yv','uc','vc','wc']:
                    if extra in ds: del ds[extra]

                type2_sel= ds['boundary_S'].isel(Nk=0,Nt=0)==0.0
                ds=ds.isel(Ntype2=type2_sel)
                del ds['boundary_S']
                dss.append(ds)

            trim_dss=[]
            for ds in dss:
                if not trim_dss:
                    trim_dss.append(ds)
                    continue
                else:
                    t_sel=ds.time.values>trim_dss[-1].time.values[-1]
                    if t_sel.sum():
                        trim_dss.append(ds.isel(Nt=t_sel))
                    else:
                        log.warning("BC dataset had no useful times?")

            ds=xr.concat(trim_dss,dim='Nt',data_vars='different')
            # somehow there is some 1e-9 difference between xe 
            for v in ['xe','ye']:
                if 'Nt' not in ds[v].dims: continue
                # make sure it's not a terrible issue
                assert ds[v].std(dim='Nt').max()<1.0
                ds[v]=ds[v].isel(Nt=0)
                
            ds.to_netcdf(compiled_fn)
            for model in self.hydro_models():
                model.bc_ds.close()
                model.bc_ds=None

        ds=xr.open_dataset(compiled_fn)
        return ds

    def parse_sql(self,fn):
        """
        semi-parse a text file that has -- comments, 
        key = value lines
        and other lines are returned as is (stripping whitespace), in
        a list of length 1
        """
        with open(fn, 'rt') as fp:
            def tok():
                for line in fp:
                    line=line.split('--')[0] # drop comments
                    if '=' not in line:
                        yield [line.strip()]
                    else:
                        k,v=line.split('=')
                        k=k.strip()
                        v=v.strip()
                        if v[0]==v[-1]=="'":
                            v=v[1:-1]
                        yield k,v
            tokens=list(tok())
        return tokens
    def open_binfile(self,group):
        return ptm_tools.PtmBin(os.path.join(self.run_dir,group+'_bin.out'))

    def groups(self):
        """
        list of all the group with bin output for this run
        """
        all_bins=glob.glob(os.path.join(self.run_dir,'*_bin.out'))
        return [os.path.basename(b).replace('_bin.out','') for b in all_bins]
    
    @classmethod
    def group_to_src_behavior(cls,group):
        m=re.match('(.*)_([^_]*)',group)
        return m.group(1),m.group(2)
    
    def get_Qdata_for_group(self,group):
        """
        group: name of PTM group
        returns a time series DataArray for the
        respective source's flow rate.
        """
        src_name,behavior_name=self.group_to_src_behavior(group)

        bc_ds=self.bc_ds()
        
        # get a flow time series for this specific group
        try:
            seg_i=list(bc_ds.seg_name.values).index(src_name)
            Q_time_series=bc_ds.set_coords('time')['boundary_Q'].isel(Nseg=seg_i)
        except ValueError:
            # point sources are labeled as srcNNN
            pnt_i=int(src_name.replace('src',''))
            Q_time_series=bc_ds.set_coords('time')['point_Q'].isel(Npoint=pnt_i)
        return Q_time_series
    def get_Qfunc_for_group(self,group):
        """
        thin wrapper to handle time interpolation, 
        """
        Q_time_series=self.get_Qdata_for_group(group)
        
        def Q_for_t(t,Q_time_series=Q_time_series):
            ti=np.searchsorted(Q_time_series.time.values,t)
            return Q_time_series.values[ti]
        return Q_for_t
    
# information that will be filled in by scan_group()
base_ret_dtype=[ ('x',np.float64,3), # location
                 ('group',object), # string name of group,
                 ('source',object),  # string name of source
                 ('behavior',object), # string name of behavior
                 ('part_id',np.int32), # id from fish-ptm within group
                 ('rel_time','M8[s]'), # time of release
                 ('obs_time','M8[s]'), # time of observation
                 # when this particle was released, how many others
                 # were released within the same group, per hour.
                 ('grp_rel_per_hour',np.float64),
                 ('mass',np.float64)]

def scan_group(self,group,time_range,z_range=None,grid=None,
               max_age=np.timedelta64(30,'D'),spinup=None,
               weight=1.0,
               extra_fields=[]):
    """
    group: name of the group to scan
    time_range: a pair of datetime64 denoting the range of time steps to include
    mass is returned per particle, based on the number of particles released per hour,
      the flow associated with the particular source, and the number of time steps 
      integrated over.

    z_range: a pair of z coordinates denoting what part of the water column to include.
     - this requires a grid which has per-cell z_bed with the same datum as
       the ptm run.
     - values are taken as positive-up, with bed=0.

     currently experimenting with having the z_range part in a post processing step.
    """
    if spinup is None: spinup=max_age

    ret_dtype=base_ret_dtype+extra_fields
    
    src_name,behavior_name=self.group_to_src_behavior(group)

    bf=self.open_binfile(group)

    nsteps=bf.count_timesteps()
    t0,_=bf.read_timestep(0)
    t0=utils.to_dt64(t0)

    # Array mapping particle index to a mass.
    # here mass is for a unit concentration

    # data to store for each particle -- this is an intermediate
    # data structure, not what gets returned.  it is indexed by
    # particle id, while the return value is indexed by
    # particle-sample, potentially including each particle multiple times.

    calc_dtype=[ ('rel_time','M8[s]'),
                 ('mass',np.float64),
                 ('grp_rel_per_hour',np.float64)]
    
    # mass=np.nan*np.ones(1000) # will be expanded as needed
    particles=np.zeros(1000,calc_dtype)

    # could calculate the number of particles released in an interval,
    # or just know that it's 5.
    count_per_release=5
    release_interval_s=3600 # and how often are particles released

    # accumulate per-time step observations, to be concatenated at the
    # end and returned
    ret_particles=[]

    Qfunc=self.get_Qfunc_for_group(group)

    # how many time steps of output fell within time_range (and
    # thus the denominator for averaging mass after the loop)
    n_steps_included=0

    for ti in range(nsteps):
        t,parts=bf.read_timestep(ti)
        t=utils.to_dt64(t)
        if t>=time_range[1]:
            log.debug("Read beyond the time range. Done with this group")
            break
        
        max_part_id=parts['id'].max()
        while max_part_id+1>len(particles):
            # double the size
            new=np.zeros(len(particles),calc_dtype)
            new['grp_rel_per_hour']=np.nan # mark uninitialized
            particles=np.concatenate([particles,new])

        # any particles with nan grp_rel_per_hour are assumed new
        # missing now has indices into parts['id']
        missing=np.nonzero(np.isnan(particles['grp_rel_per_hour'][parts['id']]))[0]
        # 1g/m3 * m3/s / (particles/release) * (seconds/release)
        # 1g/particle
        # mass[parts['id'][missing]]=Q / count_per_release * release_interval_s

        new_ids=parts['id'][missing] # index into particles for new particles this step
        # this could be inferred.  gets trickier with SJ, Sac and DDSD. FIX.
        grp_rel_per_hour=5.0
        particles['grp_rel_per_hour'][new_ids]=grp_rel_per_hour
        particles['rel_time'][new_ids]=t

        # Sac, SJ can have negative Q...  will have to return to that
        # later.  FIX.
        Q=max(0,Qfunc(t))
        
        # g/m3 * m3/s * s/hr / (particles/hour)
        # => g/particle
        # and added /nsteps to get an average
        particles['mass'][new_ids]=weight*Q*3600/grp_rel_per_hour
        
        if (t-t0>=spinup) and (t>=time_range[0]):
            n_steps_included += 1
            # at this point only filter on age.
            age=(t-particles['rel_time'][parts['id']])
            sel=age < max_age
            # this is where we could further filter on where the particles
            # are, further narrowing the sel array
            # but doing the z_range filter after the fact will make it easier
            # to streamline, and probably more efficient than doing it here.
            assert z_range is None,"Have not implemented z-range yet"
            
            ret=np.zeros(len(parts[sel]),ret_dtype)
            ret['x']=parts['x'][sel]
            ret['group']=group
            ret['source']=src_name
            ret['behavior']=behavior_name
            ret['part_id']=parts['id'][sel]
            ret['rel_time']=particles['rel_time'][parts['id'][sel]]
            ret['obs_time']=t
            ret['grp_rel_per_hour']=particles['grp_rel_per_hour'][parts['id'][sel]]
            ret['mass']=particles['mass'][parts['id'][sel]]
            ret_particles.append(ret)
        else:
            # already done the bookkeeping, but too early to actually output
            # particles
            pass

    msg=f"Run {self.run_dir}, group {group} time range {time_range} found no particles, with {n_steps_included} steps"
    if len(ret_particles)==0:
        log.warning(msg)
        raise Exception(msg)
        return None
    
    result=np.concatenate(ret_particles)
    # averaging in time:
    result['mass'] /= n_steps_included
    return result


# POTWs that I'm not worrying about, but are still left in
# the PTM runs.
# for the moment, I'm ignoring the Delta, too.
skip_source=['petaluma','sonoma_valley','ddsd','lg']

def query_runs(ptm_runs,group_patt,time_range,z_range=None,grid=None,
               max_age=np.timedelta64(30,'D'),
               spinup=None,
               conc_func=lambda group,source,behavior: 1.0,
               run_weights=None):
    """
    ptm_runs: List of PtmRun instancse
    groups: regular expression matching the group names of interest
    time_range: pair of datetime64s defining time period to include
    z_range: come back to this, but it will be a way to filter on 
     vertical position.
    max_age: ignore particles older than this
    spinup: don't return particles within this interval of the
     start of the run, defaults to max_age.

    conc_func: maps group, source behavior to concentration in particles/L.

    run_weights: list of factors, probably should sum to 1.0,
      determines how to scale each of the runs.  defaults to average.  that's
      not great if one run had a much larger particles_per_release.
      then again, currently particles_per_release is assumed!
    """
    if spinup is None:
        spinup=max_age

    if run_weights is None:
        # straight sum.  assumes that there is no duplication across runs
        run_weights=np.ones(len(ptm_runs))/float(len(ptm_runs))

    if not group_patt.endswith('$'):
        group_patt+="$"
    if not group_patt.startswith('^'):
        group_patt='^' + group_patt
        
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
            if src_name in skip_source:
                log.info(f"Will skip source {src_name} -- it's in skip_source")
                continue
            # here convert conc from particles/L to particles/m3.
            # this corresponds to the switch to v03 outputs.
            conc=1000*conc_func(group,src_name,behavior_name)
            if conc==0.0:
                log.info(f"Will skip source {src_name}, behavior {behavior_name}, its concentration is 0")
                continue
            groups.append(group)
            
            part_obs=scan_group(run,group,time_range=time_range,z_range=z_range,
                                weight=conc*run_weights[run_idx],
                                max_age=max_age,spinup=spinup,
                                grid=grid, extra_fields=[('run_idx',np.int32)])
            part_obs['run_idx']=run_idx
            
            all_part_obs.append(part_obs)
    result=np.concatenate(all_part_obs)
    assert np.isnan(result['grp_rel_per_hour']).sum()==0

    # Convert to xarray Dataset for easier processing down the line.
    ds=xr_utils.structure_to_dataset(result,'particle',{'x':('xyz',)})

    # and include some metadata
    ds['group_pattern']=(),group_patt
    ds['time_start']=(),time_range[0]
    ds['time_end']  =(),time_range[1]
    ds['max_age']=(),max_age
    ds['ptm_runs']=('ptm_run',), [p.run_dir for p in ptm_runs]
    ds['ptm_groups']=('ptm_groups',), groups
    return ds

##

# as of -v02, there is no scaling for 100%/70% of POTWs in this
# file, and instead that scaling is applied in conc_func
conc_ds=xr.open_dataset("../loads/plastic_loads-7classes-v02.nc")
# Group: UALAMEDA_up50000 UALAMEDA up50000

@memoize.memoize()
def conc_func(group,src,behavior):
    """
    For a release with group name group, parsed into src and
    behavior, return the estimated concentration in particles/l.
    """
    if behavior=='none':
        w_s=0.0
    else:
        w_s=float(behavior.replace('up','-').replace('down',''))/1e6

    # for laptop
    w_s_i=utils.nearest(conc_ds.w_s.values,w_s)
    
    if src in ['SacRiver','SJRiver']:
        v=0.0
        log.info(f"Got {src} -- returning {v}")
        return v
    else:
        source_map={'NAPA':'stormwater',
                    'COYOTE':'stormwater',
                    'SCLARAVCc':'stormwater',
                    'UALAMEDA':'stormwater',
                    'cccsd':'CCCSD',
                    'src000':'EBDA',
                    'src001':'EBMUD',
                    'sunnyvale':'SUNN',
                    'fs':'FSSD',
                    'palo_alto':'PA',
                    'src002':'SFPUC',
                    'san_jose':'SJ'}
        source=source_map[src]

    c=conc_ds.conc.isel(w_s=w_s_i).sel(source=source).item()
    # 
    if source=='stormwater':
        scale=1./0.33
    else:
        scale=1./0.70 # for wastewater
    return scale*c

def add_cells(particles,grid,overwrite=False):
    """
    add a 'cell' field to a structure ndarray filling with the cell containing
    each particle, or -1 when cell doesn't exist.
    particles: xr.Dataset with an 'x' field giving position.
    """
    if 'cell' in particles and not overwrite:
        return particles
    particles['cell']=('particle',),grid.points_to_cells(particles['x'].values[:,:2])
    return particles

def add_z_bed(particles,grid):
    particles=add_cells(particles,grid)

    # assume that whoever loaded/supplied the grid has created a
    # cells['z_bed'] field which is positive-up, and includes any
    # bathy offset

    part_z_bed=grid.cells['z_bed'][particles['cell'].values]
    valid=particles['cell'].values>=0
    part_z_bed[~valid]=np.nan
    # these are easy to calculate and just make the files bigger.
    #hab=particles['x'].values[:,2] - part_z_bed
    #particles['hab']=('particle',), hab
    particles['z_bed']=('particle',), part_z_bed
    
    return particles

class EtaExtractor(object):
    def __init__(self,ptm_runs):
        self.hydro_fns=ptm_runs[0].ptm_hydro_files()
        self.fn_index=-1
        self.ds=None
        assert len(self.hydro_fns),"Got zero hydro_fns"
        self.open_index(0)
    def open_next(self):
        return self.open_index(self.fn_index+1)
    def open_prev(self):
        return self.open_index(self.fn_index-1)
    def open_index(self,index):
        if (index<0) or (index>=len(self.hydro_fns)): return False
        if self.ds: 
            self.ds.close()
            self.ds=None
        self.ds=xr.open_dataset(self.hydro_fns[index])
        self.fn_index=index
        # more standard name for time
        self.ds['time']=self.ds['Mesh2_data_time']  
        return True
    def eta_for_time(self,t):
        while t>self.ds.time.values[-1]:
            if not self.open_next():
                return None
        while t<self.ds.time.values[0]:
            if not self.open_prev():
                return None
        ti=np.searchsorted(self.ds.time.values,t)
        if self.ds.time.values[ti]!=t:
            print(f"Time mismatch: {self.ds.time.values[ti]} (ds) != {t} (requested)")
        else:
            pass # print("Time matches")
        return self.ds.Mesh2_sea_surface_elevation.isel(nMesh2_data_time=ti).values

def add_z_surf(particles,grid,ptm_runs): 
    # process these in chronological order
    add_cells(particles,grid)
    z_surf=np.nan*np.ones(particles.dims['particle'])
    extractor=EtaExtractor(ptm_runs)
    for t,idxs in utils.enumerate_groups(particles['obs_time'].values):
        eta=extractor.eta_for_time(t)
        z_surf[idxs]=eta[particles['cell'].values[idxs]]
    z_surf[ particles['cell'].values<0 ] = np.nan
    particles['z_surf']=('particle',),z_surf
    return particles

def add_z_info(particles,grid,ptm_runs=None):
    particles=add_z_bed(particles,grid)
    if ptm_runs is not None:
        particles=add_z_surf(particles,grid,ptm_runs)
    return particles

def filter_by_z_range(particles,z_range,grid,ptm_runs=None):
    if ptm_runs is None: 
        assert (z_range[0]>=0) and (z_range[1]>0.0), "Surface referenced requires ptm_runs"
    particles=add_z_info(particles,grid,ptm_runs=ptm_runs)
    z=particles['x'].values[:,2]
    z_to_bed=z-particles.z_bed.values # should be >=0
    if 'z_surf' in particles:
        z_to_surf=z-particles.z_surf.values # should <=0

    with np.errstate(invalid='ignore'):
        # lower bound:
        sel=np.isfinite(z_to_bed)
        # choice of < vs. <= is quirky, but trying
        # to ensure that a lower bound of 0 includes the bed,
        # and an upper bound of 0 includes the surface, while
        # but a value used in the middle will not double count
        # particles.  not critical with any of the existing
        # analysis.
        if z_range[0]>=0:
            sel=sel&(z_to_bed>=z_range[0])
        else:
            sel=sel&(z_to_surf>z_range[0])
        if z_range[1]>0:
            sel=sel&(z_to_bed<z_range[1])
        else:
            sel=sel&(z_to_surf<=z_range[1])
    result=particles.isel(particle=sel)
    result['z_range']=('two',),z_range
    return result

def particle_to_density(particles,grid,normalize='area'):
    """
    particles: dataset with 'x' and 'mass', optionally with 'cell'
    normalize: 'area' normalize by cell areas to get a mass/area
               'volume': Not implemented yet.
    """
    mass=np.zeros(grid.Ncells(),np.float64)

    particles=add_cells(particles,grid)

    cells=particles['cell'].values
    masses=particles['mass'].values
    for i,cell in utils.progress(enumerate(cells)):
        if cell<0: continue
        mass[cell] += masses[i]

    if normalize=='area':
        mass/=grid.cells_area()
    elif normalize=='mass':
        pass
    else:
        raise Exception(f"Not ready for normalize={normalize}")
    return mass

def grid_from_ptm_hydro(grid_fn):
    ds=xr.open_dataset(grid_fn)
    ds['Mesh2']=(),1
    ds.Mesh2.attrs.update(dict(cf_role='mesh_topology',
                               node_coordinates='Mesh2_node_x Mesh2_node_y',
                               face_node_connectivity='Mesh2_face_nodes',
                               edge_node_connectivity='Mesh2_edge_nodes',
                               node_dimension='nMesh2_node',
                               edge_dimension='nMesh2_edge',
                               face_dimension='nMesh2_face',
                               face_coordinates='Mesh2_face_x Mesh2_face_y',
                               edge_coordinates='Mesh2_edge_x Mesh2_edge_y'))
    grid=unstructured_grid.UnstructuredGrid.from_ugrid(ds)
    # This is straight from the output, so no need to add bathy offset
    grid.add_cell_field('z_bed',-grid.cells['Mesh2_face_depth'])

    return grid
    

if __name__=='__main__':
    import socket
    if socket.hostname() == 'soling': # laptop
        ptm_runs=[
            PtmRun(run_dir="/home/rusty/src/sfb_ocean/ptm/all_sources/all_source_select_w_const")
        ]

        # grid=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/spliced-bathy.nc")
        # ptm_runs=[ PtmRun(run_dir="../../sfb_ocean/ptm/all_source/20180115/w-0.0005") ]

        if 0:
            # May not be the right grid -- would be better to copy in a grid from
            # one of the original ptm hydro paths
            grid=unstructured_grid.UnstructuredGrid.from_ugrid("../../sfb_ocean/suntans/grid-merged/spliced_grids_01_bathy.nc")
            # include a z_bed field which matches the hydro that was used to run the PTM.
            # in this case, I know that cells['depth'] in the above file is positive up, no shift, no clip.
            # also know that the minimum depth is always 0.1, and bathy is clipped to that before being written.
            # would be better to read in the real thing..
            grid.add_cell_field('z_bed',(grid.cells['depth'] - 5.0).clip(-np.inf,-0.1))
        else:
            grid_fn="/home/rusty/src/sfb_ocean/suntans/runs/merged_018_20171227/ptm_average.nc_0000.nc"
            ds=xr.open_dataset(grid_fn)
            ds['Mesh2']=(),1
            ds.Mesh2.attrs.update(dict(cf_role='mesh_topology',
                                       node_coordinates='Mesh2_node_x Mesh2_node_y',
                                       face_node_connectivity='Mesh2_face_nodes',
                                       edge_node_connectivity='Mesh2_edge_nodes',
                                       node_dimension='nMesh2_node',
                                       edge_dimension='nMesh2_edge',
                                       face_dimension='nMesh2_face',
                                       face_coordinates='Mesh2_face_x Mesh2_face_y',
                                       edge_coordinates='Mesh2_edge_x Mesh2_edge_y'))
            grid=unstructured_grid.UnstructuredGrid.from_ugrid(ds)
            # This is straight from the output, so no need to add bathy offset
            grid.add_cell_field('z_bed',-grid.cells['Mesh2_face_depth'])

        group_patt='.*_down2000'

    elif socket.hostname() == 'cws-linuxmodeling':
        # Ultimately the interface is probably something along the lines of
        ptm_runs=[
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w-0.05"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w-0.005"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w-0.0005"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w0.0"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w0.0005"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w0.005"),
            PtmRun(run_dir="/opt2/sfb_ocean/ptm/all_source/20170615/w0.05"),
        ]
        # May not be the right grid -- would be better to copy in a grid from
        # one of the original ptm hydro paths
        grid=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/spliced-bathy.nc")
        group_patt='.*_down500'

    else:
        raise Exception(f"unknown host {socket.gethostname()}")

    # A 1 hour window gives 27k particles
    part_obs=query_runs(ptm_runs,
                        group_patt=group_patt,
                        time_range=[np.datetime64("2017-07-30 00:00"),
                                    np.datetime64("2017-07-30 03:00")],
                        z_range=None, # not ready
                        max_age=np.timedelta64(50,'D'),
                        conc_func=conc_func,
                        grid=grid)


    # Filter by z_range in a separate step, so there's the possibility
    # of more speedup.
    z_range=[0,0.5]

    ## 
    plt.figure(1).clf() ; plt.hist(part_hab,200)
    # plt.figure(1).clf() ; plt.hist(part_z_bed,200)# looks fine
    # plt.figure(1).clf() ; plt.hist(part_obs['x'][:,2],200)

    ## That's showing a lot of cells that are below the bed.  wtf?

    # maybe those are cells that diffused horizontally too much?
    bad=np.isfinite(part_hab) & (part_hab<-1)
    good=np.isfinite(part_hab) & (part_hab>=0)

    ##

    # Look at the history of a particle that's below the bed
    bad_i=np.nonzero(bad)[0][0]
    print(f"part idx {bad_i}: z={part_obs['x'][bad_i,2]} z_bed={part_z_bed[bad_i]}  hab={part_hab[bad_i]}")
    bad_part=part_obs['part_id'][bad_i]
    for name in part_obs.dtype.names:
        print(f"  {name}: {part_obs[name][bad_i]}")

    ##

    ptm_run=ptm_runs[part_obs['run_idx'][bad_i]]

    bf=ptm_run.open_binfile(part_obs['group'][bad_i])

    ##

    nsteps=bf.count_timesteps()

    times=[]
    recs=[]

    for s in utils.progress(range(nsteps)):
        t,p=bf.read_timestep(s)
        sel=np.nonzero( p['id']==bad_part)[0]
        if len(sel):
            assert len(sel)==1,"what?"
            recs.append( p[sel[0]].copy() )
            times.append(utils.to_dt64(t))
    traj=np.array(recs)

    traj=add_cells(traj,grid)

    t=np.array(times)

    ##
    plt.figure(3).clf()
    ccoll=grid.plot_cells(values=grid.cells['z_bed'],cmap='jet',clim=[-20,0])

    plt.plot(traj['x'][:,0],
             traj['x'][:,1],
             'k-',zorder=1)
    scat=plt.scatter(traj['x'][:,0],
                     traj['x'][:,1],
                     20,
                     traj['x'][:,2],cmap='jet',zorder=2)
    scat.set_clim([-20,0])
    plt.axis('equal')
    plt.axis( (575943., 581457., 4211703., 4215158.) )

    # So this specific one starts at CCCSD, sloshes over to the ghost fleet,
    # and weasles its way onto the shoal, and sometimes out of the grid.

    ##

    plt.figure(4).clf()
    plt.plot(t,traj['x'][:,2],'g',label='Particle z')
    plt.plot(t,grid.cells['z_bed'][traj['cell']],'k',label='Cell z_bed')
    plt.legend()

    ##

    plt.figure(2).clf()

    ccoll=grid.plot_cells(values=grid.cells['z_bed'],cmap='jet',clim=[-20,0])
    scat=plt.scatter(part_obs['x'][good,0],
                     part_obs['x'][good,1],
                     20,
                     part_obs['x'][good,2],
                     cmap='jet')
    scat.set_clim([-20,0])

    plt.colorbar(ccoll)
    plt.axis('equal')


    ##
    # that returned 4.8M points.  with output for a single time step
    # via time_range, that becomes 41k.
    # 18 groups.
    # 30 days in, 5 particles/hour/group
    # 30 days*24 hr/day * 5 particles/grp/hr *18 grps
    # but time range includes 5 days,

    particles=part_obs


    conc0=particle_to_density(particles,grid,normalize='area')

    M=grid.smooth_matrix(f=0.5,dx='grid',A='grid',V='grid',K='scaled')

    conc=conc0
    for _ in range(20):
        conc=M.dot(conc)


    plt.figure(1).clf()
    fig,axs=plt.subplots(1,2,sharex=True,sharey=True,num=1)

    ax=axs[0]
    age_secs=(part_obs['obs_time']-part_obs['rel_time'])/np.timedelta64(1,'s')
    if 0: # scale by mass, color by age
        scat=ax.scatter(part_obs['x'][:,0],part_obs['x'][:,1],
                        part_obs['mass']/2e4,
                        age_secs/86400.,
                        cmap='jet')
        label='Age (days)'
    if 1: # color by log(mass)
        scat=ax.scatter(part_obs['x'][:,0],part_obs['x'][:,1],
                        10,
                        np.log(part_obs['mass'].clip(100,np.inf)),
                        cmap='jet')
        label='log(mass)'
    grid.plot_edges(lw=0.5,color='k',zorder=-1,ax=ax)

    plt.colorbar(scat,label=label,ax=ax)

    ax=axs[1]

    from matplotlib import colors

    ccoll=grid.plot_cells(values=conc.clip(1e-5,np.inf),ax=ax,lw=0.5,cmap='jet',
                          norm=colors.LogNorm(vmin=1e-5,vmax=1,clip=True))
    ccoll.set_edgecolor('face')

    plt.colorbar(ccoll,label='Conc',ax=ax)

    fig.savefig('sample-output.png',dpi=120)

    # plt.setp(axs,aspect='equal')


