import os
import pandas as pd
import glob
import time
import six
from sqlite3 import dbapi2 as sql

import xarray as xr
import numpy as np
from stompy.grid import unstructured_grid
from stompy import utils
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from stompy.model.fish_ptm import ptm_tools
from stompy import memoize
##

grid=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/splice-merge-05-filled-edit70.nc")

run_dir="/home/rusty/src/sfb_ocean/ptm/all_sources/all_source_select_w_const"
groups=glob.glob(os.path.join(run_dir,"*_bin.out"))

pbf=ptm_tools.PtmBin(groups[0])

## 

def grid_to_cell_cache_fn(self,grid):
    grid_key=memoize.memoize_key( grid.nodes['x'][ grid.cells['nodes'].clip(-1) ] )
    return self.fn+".cells-"+grid_key+".nc"

def precompute_cells(self,grid,force=False):
    """
    The cached cell info is about 15% the size of the original
    bin.out file. ad-hoc binary and netcd yield about the same
    file size.

    grid: UnstructuredGrid. Will compute the cell index (or -1)
       for each particle at each time step.

    force: when False, use cached data when possible, otherwise
      recompute.

    returns the cached filename, a netcdf file.
    """
    cell_cache_fn=grid_to_cell_cache_fn(self,grid)

    n_steps=self.count_timesteps()
    dnums=[]
    xys=[]

    # loop once to gather all points
    for ts in utils.progress(range(n_steps)):
        dnum,parts=pbf.read_timestep(ts)
        if ts%100==0:
            print(f"{ts} {dnum} {len(parts)}")
        dnums.append(dnum)
        xys.append( parts['x'][:,:2].copy() )
    all_xy=np.concatenate(xys)

    # compute cells:
    t=time.time()
    # be sure to insert these as regular int
    all_cell=grid.points_to_cells(all_xy)
    elapsed=time.time() - t
    print("Python mapping time: %.3fs"%elapsed)

    ds=xr.Dataset()
    ds['cell']=('particle_loc',),all_cell.astype(np.int32)
    ds['dnum']=('time',),utils.to_dt64(np.array(dnums))
    counts=np.array( [len(xy) for xy in xys] )
    ds['count']=('time',),counts
    ds['offset']=('time',),np.cumsum(counts)-counts
    ds.to_netcdf(cell_cache_fn,mode='w')
        
    return cell_cache_fn

cell_cache_fn=precompute_cells(pbf,grid)

##

# what about a fairly intensive preprocessing step, that takes 
# advantage of the non-overlapping, convex polygon nature of
# the grid cells.

#grid=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/splice-merge-05-filled-edit70.nc")

# have grid and all_xy

# want to come up with a sequence of lines

centroids=grid.cells_centroid()

## 
js=np.random.choice(grid.Nedges(),size=50)
edge_norm=grid.edges_normals()

def partition_by_j(j,xy):
    p=grid.nodes['x'][ grid.edges['nodes'][j,0], :]
    # for now, not worrying about the edge cases
    return np.dot(xy-p,edge_norm[j])>0

def partition_score(j,xy):
    parts=partition_by_j(j,xy)
    return (parts.sum()/len(parts)-0.5)**2

scores= [partition_score(j,xy)
         for j in js]
j_sel=js[ np.argmin(scores) ]

##

# so we end up with a binary tree.
levels=20

tree_shape=2**(levels+1)

tree_node_dtype=[ ('level',np.int8),
                  ('j',np.int32),
                  ('vec',(np.float64,2)),
                  ('thresh',np.float64),
                  ('cell',np.int32)]
                  
tree=np.zeros( tree_shape, dtype=tree_node_dtype )

tree['cell']=-1
    
# the points we're trying to place
xy=centroids # test/build with centroids
# everyone starts at the root of the tree.
tree_locs=np.ones(len(xy),np.int32)

# This is a solid start, but the building process needs
# to take into account the full polygons.
# that of course will make the index building vastly slower.
# Come back to this later.  Not worth the investment right now.

for level in range(levels):
    counts=[]
    # iterate over the nodes at this level...
    for node in range(2**level, 2**(level+1)):
        tree['level'][node]=level
        parent=node//2
        
        if tree['cell'][parent]>=0:
            # already know that it's a leaf
            cells_for_node=[tree['cell'][parent]]
            # so just copy same test
            tree['cell'][node]=tree['cell'][parent]
            
            tree['j'][node]=tree['j'][parent]
            tree['vec'][node]=tree['vec'][parent]
            tree['thresh'][node]=tree['thresh'][parent]
            count=1
        else:
            # which points are on this node:
            cells_for_node=np.nonzero( tree_locs==node )[0]
            count=len(cells_for_node)
            
            if count==1:
                tree['j'][node]=tree['j'][parent]
                tree['vec'][node]=tree['vec'][parent]
                tree['thresh'][node]=tree['thresh'][parent]
                tree['cell'][node]=cells_for_node[0]
            elif count==0:
                assert False, "Shouldn't have a 0 child of a non-determined parent"
            else:
                # have to make a choice of a good subdividing line
                # get some candidate edges from a subset of cells
                sel_cells=np.random.choice(cells_for_node,size=min(50,len(cells_for_node)))
                js=[ j for c in sel_cells for j in grid.cell_to_edges(c)]
                scores=[partition_score(j,xy[cells_for_node])
                        for j in js]
                j=js[ np.argmin( scores ) ]
                p=grid.nodes['x'][grid.edges['nodes'][j,0]]

                tree['j'][node]=j
                tree['vec'][node]=edge_norm[j]
                tree['thresh'][node]=np.dot(edge_norm[j],p)
            
        counts.append(count)

    # here we'd do a bulk comparison on nodes, and update all tree locs.
    result= (xy*tree['vec'][tree_locs]).sum(axis=1) > tree['thresh'][tree_locs]
    tree_locs*=2 # everyone advances to the next level..
    tree_locs[result]+=1 # and some get the other leaf

    counts=np.array(counts)
    print(f"Level {level}: max count {counts.max()}  mean {counts.mean():.1f}")

## 

# just as a rough check on computational complexity
# there are 2**16 cells.
# so at a minimum, 16 partitioning steps.
# each time, all the points
t=time.time()
for j in js[:16]:
    partition_by_j(j,all_xy)
print("Lower bound on linear partitioning: %.3fs"%( time.time() - t))

# that's 700k points. when all of them go up against
# the same edge on each round, it's 0.146s

##

t=time.time()
jidxs=np.zeros(len(all_xy),np.int32)
for _ in range(25):
    ps=grid.nodes['x'][ grid.edges['nodes'][jidxs,0], :]
    # for now, not worrying about the edge cases
    delta=all_xy-ps
    norms=edge_norm[jidxs,:]
    dot=delta[:,0]*norms[:,0] + delta[:,1]*norms[:,1]
    result=dot>0
# 0.4s
print("Lower bound on linear partitioning: %.3fs"%( time.time() - t))

##
t=time.time()
tree_locs=np.ones(len(all_xy),np.int32)

for _ in range(levels):
    result= (all_xy*tree['vec'][tree_locs]).sum(axis=1) > tree['thresh'][tree_locs]
    tree_locs*=2 # everyone advances to the next level..
    tree_locs[result]+=1 # and some get the other leaf
lin_cells=tree['cell'][tree_locs]
print("Linear partitioning: %.3fs"%( time.time() - t))

##

t=time.time()
all_cell=grid.points_to_cells(all_xy)
print("Existing approach: %.3fs"%( time.time() - t))
