import re
import time
import os
from collections import defaultdict
import logging as log

import xarray as xr
import numpy as np
from sqlite3 import dbapi2 as sql

from stompy import utils, memoize

import postprocess_v00 as post

def make_load_table(con,load_name,load_fn,msg=log.debug,clean=True,
                    stormwater_scale=1/0.33,
                    wastewater_scale=1/0.70):
    """
    bring concentrations into db as a temporary, in-memory
    table.
    This version assumes the old runs, where only a small
    number of stormwater sources are included.

    clean: if True, the memory database is started fresh.  otherwise, just
    add in the one table (replacing it if that specific table already exists)
    
    the new table is in the schema load, with table name = load_name
    """
    load_ds=xr.open_dataset(load_fn)
    
    curs=con.cursor()
    curs.execute("select name from ptm_group group by name")
    group_names=curs.fetchall()

    behavior_to_ws={'down50000':0.05,
                    'down5000':0.005,
                    'down500':0.0005,
                    'none':0.0,
                    'up500':-0.0005,
                    'up5000':-0.005,
                    'up50000':-0.05}

    # almost everything is stormwater, so just create a default
    # dict and explicitly name the non-stormwater.
    source_map=defaultdict(lambda:'stormwater')
    source_map['cccsd']='CCCSD'
    source_map['sunnyvale']='SUNN'
    source_map['fs']='FSSD'
    source_map['palo_alto']='PA'
    source_map['san_jose']='SJ'
    source_map['src000']='EBDA'
    source_map['src001']='EBMUD'
    source_map['src002']='SFPUC'

    # these we don't actually use
    source_map['petaluma']='SKIP'

    # These shouldn't be used, but including just to be sure
    # that if they somehow show up, they won't contaminate
    # stormwater.
    source_map['SacRiver']='DELTA'
    source_map['SJRiver']='DELTA'

    # tuples of group_name, particles/m3.
    # can omit rows with 0 concentration.
    conc_rows=[] 

    for group_name, in group_names:
        m=re.match(r'(.*)_(down\d+|up\d+|none)(_rel.*)?',group_name)
        source=m.group(1)
        behavior=m.group(2)
        rel_time=m.group(3) # may be missing
        w_s=behavior_to_ws[behavior]

        source_name=source_map[source]
        if source_name in ['DELTA','SKIP']:
            conc=0.0
        else:
            conc=load_ds.conc.sel(source=source_name,w_s=w_s).item()

        if source_name=='stormwater':
            conc*=stormwater_scale
        else:
            conc*=wastewater_scale

        # load netcdf is in particles/l, but we want to set the calculation
        # up as particle/m3. Updated 2019-11-17
        conc*=1000
        if conc>0.0:
            conc_rows.append( (group_name,conc) )
        msg(f"{source:15s}  {behavior:9s} => {source_name:15s} {conc:.4f}")

    # --- load that into a table
    if clean:
        try:
            curs.execute("DETACH load")
        except sql.OperationalError:
            print("memory table wasn't there.  no worries.")
        
    try:
        curs.execute("ATTACH ':memory:' as load")
    except sql.OperationalError:
        print("maybe the memory table was already attached?")

    curs.execute(f"""DROP TABLE if exists load.{load_name}""")

    curs.execute(f"""CREATE TABLE load.{load_name} (
       group_name text,
       part_per_m3 double);""")
    curs.executemany(f"""INSERT into load.{load_name} (group_name,part_per_m3) VALUES (?,?)""",
                    conc_rows)
    con.commit()
    load_ds.close()



class PtmSet(object):
    cache_dir="/opt2/sfb_ocean/ptm/all_source/queries"
    cache_ver="v002"
    
    # A sample run for this group of ptm runs, used to get
    # a proper grid
    base_ptm_run_dir="/opt2/sfb_ocean/ptm/all_source/20170715/w0.0"
    load_fns={'std':"../loads/plastic_loads-7classes-v03.nc",
              'nofiber':"../loads/plastic_loads-7classes-v03-nofiber.nc",
              'fiber':"../loads/plastic_loads-7classes-v03-fiber.nc",
              'fiber_bundle':"../loads/plastic_loads-7classes-v03-fiber_bundle.nc",
              'film':"../loads/plastic_loads-7classes-v03-film.nc",
              'foam':"../loads/plastic_loads-7classes-v03-foam.nc",
              'fragment':"../loads/plastic_loads-7classes-v03-fragment.nc",
              'sphere':"../loads/plastic_loads-7classes-v03-sphere.nc"
    }
    
    # new-ish runs
    # "/opt2/sfb_ocean/ptm/all_source_020/20170715/w0.0"
    # probable need to alter load_fns, or at least how they are setup above.
    
    # List of paths to sql databases
    databases=[]
    
    # for named z_filters, how thick the layer is
    # Old Run settings
    z_thickness=0.5
    stormwater_scale=1/0.33
    wastewater_scale=1/0.70

    def __init__(self,**kw):
        utils.set_keywords(self,kw)
        self.cons={} # map database path to connection

    @memoize.imemoize()
    def base_ptm_run(self):
        return post.PtmRun(run_dir=self.base_ptm_run_dir)
    @memoize.imemoize()
    def grid(self):
        return self.base_ptm_run().grid()
    @memoize.imemoize()
    def grid_template_ds(self):
        return self.grid().write_to_xarray()
    
    @memoize.imemoize()
    def poly(self):
        return self.grid().boundary_polygon()
    @memoize.imemoize()
    def Msmooth(self):
        return self.grid().smooth_matrix(f=0.5,dx='grid',A='grid',V='grid',K='scaled')
    def smooth(self,c): 
        # This may get refactored in some other way.
        M=self.Msmooth()
        for _ in range(20):
            c=M.dot(c)
        return c

    def db_to_con(self,db):
        if db not in self.cons:
            con=sql.connect(db)
            self.cons[db]=con
        return self.cons[db]

    def conc_query(self,**kw):
        # naive approach:
        ds=self.grid_template_ds().copy()
        conc2d=np.zeros(self.grid().Ncells(),np.float64)
        for db in self.databases:
            print("Trying: %s"%db)
            one_db=self.conc_query_one_db(db,**kw)
            conc2d[:] += one_db.conc2d.values
        ds['conc']=('face',),conc2d
        ds['conc'].attrs['units']='particles m-2'
        
        # specific to this call
        ds['conc'].attrs['z_filter']=kw.get('z_filter','none')
        ds['conc'].attrs['t_start']=kw.get('t_start',-1)
        ds['conc'].attrs['t_stop']= kw.get('t_stop',-1)
        ds['conc'].attrs['grp_filter']=kw.get('grp_filter','none')
        ds['conc'].attrs['max_age']=kw.get('max_age','none')
        ds['conc'].attrs['loads']=kw.get('loads','default')
                
        return ds
    
    def conc_query_one_db(self,db,t_start,t_stop,z_filter=None,
                          grp_filter="",max_age=None,loads='std'):
        """
        db: path to database file
        
        t_start: np.datetime64 for min value, inclusive of observed time
        t_stop:  np.datetime64 for max value, exclusive of observed time
        
        z_filter: None, 'all','surf','bed', or a sql clause.
        grp_filter: "" or a sql clause to match on grp.name
        """
        epoch_start=int( utils.to_unix( t_start ) )
        epoch_stop =int( utils.to_unix( t_stop ) )

        if z_filter is None or z_filter=="all":
            z_filter=""
        elif z_filter=='bed':
            z_filter=f"and loc.z_from_bed<{self.z_thickness}"
        elif z_filter=='surf':
            z_filter=f"and loc.z_from_surface>{-self.z_thickness}"
        else:
            print(f"Will use z_filter as provided: {z_filter}")
            
        if max_age is None:
            max_age_clause=""
        else:
            max_age_clause=f"and (loc.time-rel.time)<{int(max_age/np.timedelta64(1,'s'))}"

        query=f"""
           select min(loc.cell),sum(load.part_per_m3 * rel.volume / rel.count )
            from particle_loc as loc, particle as p, ptm_release as rel, 
                 ptm_group as grp, load.{loads} as load
            where loc.time>={epoch_start} and loc.time < {epoch_stop}
              and loc.particle_id=p.id
              and p.release_id=rel.id
              and rel.group_id=grp.id
              and grp.name=load.group_name 
              and loc.cell>=0
              {max_age_clause}
              {z_filter}
              {grp_filter}
            group by loc.cell"""

        ds=xr.Dataset()
            
        # enough metadata for caching
        ds['query']=(),query
        ds['db']=(),db
        ds['method']=(),"conc_query_one_db"
        
        # note that the choice of loads is in the query, but only
        # by name, not filename.
        cache_key=memoize.memoize_key([query,db,
                                       "conc_query_one_db"])
        cache_fn=os.path.join(self.cache_dir,f"query-results-{self.cache_ver}-{cache_key}.nc")
                                      
        if os.path.exists(cache_fn):
            # print("Reading from cache: %s"%cache_fn)
            ds=xr.open_dataset(cache_fn)
            ds.load()
            ds.close()
            return ds
        
        con=self.db_to_con(db)

        self.setup_loads(con,loads)

        curs=con.cursor()
        print("Query")
        print(query)
        print("DB: ",db)
        t=time.time()
        curs.execute(query)        
        data=np.array( curs.fetchall() )
        print("Query time to return data %.2fs"%(time.time()-t))
        
        cell_count=np.zeros(self.grid().Ncells(),np.float64)

        if len(data)>0:
            # to get a proper average, need to know how many time steps that covered
            curs.execute("""select count(1),
                                   DATETIME(min(time),'unixepoch'),
                                   DATETIME(max(time),'unixepoch')
                                   from 
                  (select distinct loc.time as time
                     from particle_loc as loc
                    where loc.time>=? and loc.time< ?)""",
              [epoch_start,epoch_stop])
            count,t_start_real,t_stop_real=curs.fetchone()
            cell_count[data[:,0].astype(np.int32)]=data[:,1]
            cell_2dconc=cell_count/self.grid().cells_area()/count
        else:
            cell_2dconc=0*cell_count # redundant, just to be clear.
            
        ds['conc2d']=('face',),cell_2dconc
        
        if cache_fn is not None:
            cache_folder=os.path.dirname(cache_fn)
            os.makedirs(cache_folder,exist_ok=True)
            ds.to_netcdf(cache_fn)
        return ds

    def setup_loads(self,con,load_name,table_name=None,**kw):
        table_name=table_name or load_name
        make_load_table(con,table_name,self.load_fns[load_name],
                        stormwater_scale=self.stormwater_scale,
                        wastewater_scale=self.wastewater_scale,
                        **kw)
    
class PtmSetNew(PtmSet):
    # for named z_filters, how thick the layer is
    z_thickness=0.25
    stormwater_scale=1.0 # or a bit bigger
    wastewater_scale=1/0.70
    cache_dir="/opt2/sfb_ocean/ptm/all_source_020b/queries"
    
