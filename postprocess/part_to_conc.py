# Make a figure demonstrating the mapping of particles to concentrationo
from stompy.model.fish_ptm import ptm_tools
from stompy.grid import unstructured_grid
import matplotlib.pyplot as plt
import os
from stompy.plot import plot_wkb

import numpy as np

##

ptm_out="/home/rusty/src/sfb_ocean/ptm/all_sources/all_source_select_w_const/SacRiver_none_bin.out"
pbf=ptm_tools.PtmBin(ptm_out)

##

g=unstructured_grid.UnstructuredGrid.from_ugrid('/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/splice-merge-05-filled-edit70.nc')
M=g.smooth_matrix()

poly=g.boundary_polygon()

##

time,parts1=pbf.read_timestep(1000)
time,parts2=pbf.read_timestep(1100)

parts=np.concatenate( [parts1,parts2] )

cells=g.points_to_cells(parts['x'][:,:2],method='cells_nearest')
valid=cells>=0

cells=cells[valid]
parts=parts[valid]

conc_raw=np.bincount(cells,minlength=g.Ncells()) / g.cells_area()

##

#zoom=(577497.8704259364, 590985.558637775, 4214418.355955187, 4224534.122114066)
zoom=(574842.4818092306, 588330.1700210692, 4210498.496568622, 4220614.262727501)
## 
import contextlib

fig_dir='to_conc'
os.makedirs(fig_dir,exist_ok=True)

@contextlib.contextmanager
def F(n):
    plt.figure(n).clf()
    fig,ax=plt.subplots(1,1,num=n)
    g.plot_edges(color='k',lw=0.5,zorder=2,alpha=0.25,clip=zoom)
    plot_wkb.plot_wkb(poly,ax=ax,facecolor='none',edgecolor='k',lw=1)
    
    try:
        yield fig,ax
    finally:
        fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
        ax.axis('off')
        ax.axis('equal')
        ax.axis( zoom )
        fig.savefig(os.path.join(fig_dir,"frame%02d.png"%n))

with F(0) as (fig,ax):    
    ax.plot(parts['x'][:,0],
            parts['x'][:,1],
            'bo')

##     

clim=[0,3e-5]

with F(1) as (fig,ax):
    ccoll=g.plot_cells(values=conc_raw,cmap='CMRmap_r',clim=clim,clip=zoom)

##

for i,n in enumerate([2,10,40]):
    conc=conc_raw
    for _ in range(n):
        conc=M.dot(conc)

    with F(2+i) as (fig,ax):
        g.plot_cells(values=conc,cmap='CMRmap_r',clim=clim,clip=zoom)
