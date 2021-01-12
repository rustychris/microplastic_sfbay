"""
Classes for configurable plotting of MP concentration 
data in SF Bay.
"""
from stompy import utils
import numpy as np

from stompy.grid import unstructured_grid
import matplotlib.pyplot as plt
from matplotlib import colors
from stompy.plot import plot_wkb
import stompy.plot.cmap as scmap

cmap=scmap.load_gradient('turbo.cpt')

class BayConcFigure(object):
    figsize=(8.4,8)
    ax=None
    fig=None
    vmin=1e-4
    vmax=100.0
    zoom=(517521., 609000., 4139744., 4230000.)
    cmap=cmap
    cax_loc=[0.7,0.25,0.03,0.35]
    txt_loc=[0.65,0.7] # in ax coords
    cbar_label="Particles/m$^2$"
    cbar_args={} # don't modify - replace.
    draw_boundary=True
    fontsize=14
    extra_text=[]
    grid=None
    grid_poly=None
    num=None
    
    def __init__(self,ds,**kw):
        """
        ds: xr.Dataset with conc field over cells.
        should either have ugrid grid, or must specify grid via keyword.
        """
        utils.set_keywords(self,kw)
        self.ds=ds
        if self.grid is None:
            self.grid=unstructured_grid.UnstructuredGrid.from_ugrid(ds)
        if self.grid_poly is None:
            self.grid_poly=self.grid.boundary_polygon()
            
        conc=self.ds['conc'].values

        if self.fig is None:
            self.fig=plt.figure(figsize=self.figsize,num=self.num)
        if self.ax is None:
            self.ax=self.fig.add_subplot(1,1,1)
        
        self.ccoll=self.grid.plot_cells(values=conc.clip(self.vmin,self.vmax),
                                        cmap=self.cmap,norm=colors.LogNorm(vmin=self.vmin,vmax=self.vmax),
                                        edgecolor='face',lw=0.4,ax=self.ax)
        if self.draw_boundary:
            self.boundary=plot_wkb.plot_wkb(self.grid_poly,ax=self.ax,ec='k',lw=0.5,fc='none')
        if self.cax_loc is not None:
            self.cax=self.fig.add_axes(self.cax_loc) # Need refactor
            plt.colorbar(self.ccoll,cax=self.cax,label=self.cbar_label,extend='both',
                         **self.cbar_args)
        self.ax.axis('equal')
        self.ax.axis(self.zoom)
        self.ax.xaxis.set_visible(0)
        self.ax.yaxis.set_visible(0)
        self.fig.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.01)
        
        self.add_labels()
        
    def __del__(self):
        try:
            self.ds.close()
        except AttributeError:
            pass
    def add_labels(self):
        texts=self.behavior_label()
        texts+=self.average_label()
        texts+=self.date_label()
        texts+=self.age_label()
        texts+=self.extra_text
        self.ax.text(self.txt_loc[0],self.txt_loc[1],"\n".join(texts),
                     fontsize=self.fontsize,va='top',ha='left',transform=self.ax.transAxes)
    def behavior_label(self):
        # go from a list of groups to a label
        grp_filter=self.ds.conc.attrs['grp_filter']
        if grp_filter not in ["","none"]:
            return [self.ds.conc.attrs['grp_filter']]
        else:
            return []
        #if behavior=='none':
        #    label='Passive'
        #elif behavior.startswith('up'):
        #    w_mmps=float(behavior.replace('up',''))/1000.0
        #    label=f'Rise {w_mmps:.1f} mm/s'
        #elif behavior.startswith('down'):
        #    w_mmps=float(behavior.replace('down',''))/1000.0
        #    label=f'Settle {w_mmps:.1f} mm/s'
    def age_label(self):
        return [ 
            "Max age: %d days"%(self.ds.conc.attrs['max_age']/np.timedelta64(1,'D')) 
             ]    
    def date_label(self):
        return [self.date_label_ext(self.ds)]
    @classmethod
    def date_label_ext(cls,ds):
        t_start=ds.conc.attrs['t_start']
        t_stop =ds.conc.attrs['t_stop']
        
        def fmt_t(t): return utils.to_datetime(t).strftime('%Y-%m-%d')
        return f"{fmt_t(t_start)} – {fmt_t(t_stop)}"
    
    def average_label(self):
        return [self.average_label_ext(self.ds)]
    
    @classmethod
    def average_label_ext(cls,ds):
        if ds.conc.attrs['z_filter']=='bed':
            return "Near bed"
        elif ds.conc.attrs['z_filter']=='surf':
            return "Near surface"
        elif ds.conc.attrs['z_filter']=='all':
            return "Full depth"
        else:
            return "Vertical range: %s"%(ds.conc.attrs['z_filter'])

# Same idea but settings for coastal region.
# note that with the 15 day output, these get truncated
# and show little action in the ocean.
class CoastalConcFigure(BayConcFigure):
    figsize=(7,8.1)
    zoom=(345000., 613202., 4050000., 4230105.)
    cax_loc=[0.05,0.20,0.03,0.35]
    txt_loc=[0.05,0.13] # in ax coords
    
