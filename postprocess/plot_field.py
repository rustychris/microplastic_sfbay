# Show the field data in relation to
#

import matplotlib.pyplot as plt 
import pandas as pd
from stompy import utils
utils.path("../field_data")
import plastic_data
from stompy.spatial import proj_utils
from stompy.grid import unstructured_grid
from stompy.plot import plot_wkb
##

g=unstructured_grid.UnstructuredGrid.from_ugrid("/home/rusty/src/sfb_ocean/suntans/grid-merge-suisun/spliced-bathy.nc")
poly=g.boundary_polygon()


## 
manta=plastic_data.manta_df

#   'CB4', 'CB5', 'CB6', 'CB7', 'CB9', 'CBNMS22', 'CBNMS23', 'CBNMS24',
#   'GFNMS25', 'GFNMS26', 'GFNMS27', 'GFNMS28', 'MBNMS29', 'LSB14',
#   'MBNMS30', 'MBNMS32', 'SB10', 'SB11', 'SB13', 'SPB2', 'SPB3',
#   'CB8', 'SB12', 'LABQA', 'LSB15', 'LSB16', 'MBNMS31', 'SFBay',
#   'SUB1'

manta_volumes=pd.read_excel("../field_data/Manta Volumes (1).xlsx",sheet_name=0,skiprows=1,
                            na_values=['None','Missing'])
manta_volumes['DATE']=pd.to_datetime(manta_volumes['DATE'])

grouped=manta_volumes.groupby('SAMPLE LOCATION').mean()
grouped['lat']=0.5*(grouped['LAT START']+grouped['LAT END'])
grouped['lon']=0.5*(grouped['LONG START']+grouped['LONG END'])
manta_ll=grouped.loc[:,['lat','lon']]

# 47 don't get a lat/lon
manta=pd.merge(manta,manta_ll,left_on='StationCode',right_index=True,how='left')

print("Stations missing lat/lon in manta: %s"%
      ", ".join(manta['StationCode'][ np.isnan(manta['lat']) ].unique() ))

##
grp=manta.groupby('SampleID')

manta_per_sample=pd.DataFrame()
manta_per_sample['volume_l']=grp['Volume'].first()

manta_per_sample['count_preblank'] = grp.size()

for fld in ['lat','lon','SampleDate','Season','SampleType']:
    manta_per_sample[fld] = grp[fld].first()

# 33.5
median_blank_count=manta_per_sample['count_preblank'][ manta_per_sample['SampleType']=='FieldBlank' ].median()

manta_per_sample['count']=np.maximum(0, (manta_per_sample['count_preblank'] - median_blank_count).values )
##

# Combine DUP
to_delete=[]
for v in manta_per_sample.index.values:
    if 'DUP' in v:
        print(v)
        non_dupe=v.replace('-DUP','')
        assert non_dupe in manta_per_sample.index.values
        manta_per_sample.loc[non_dupe,'count'] += manta_per_sample.loc[v,'count']
        manta_per_sample.loc[non_dupe,'volume_l'] += manta_per_sample.loc[v,'volume_l']
        to_delete.append(v)
manta_per_sample.drop(to_delete,inplace=True)

##
manta_per_sample['part_per_m3']= 1000*manta_per_sample['count'] / manta_per_sample['volume_l']

ll=manta_per_sample.loc[: , ['lon','lat']].values
xy=np.nan*np.ones( (len(manta_per_sample),2) )

valid=np.isfinite(ll[:,0])
xy[valid,:]=proj_utils.mapper('WGS84','EPSG:26910')(ll[valid])

manta_per_sample['x']=xy[:,0]
manta_per_sample['y']=xy[:,1]

manta_per_sample.to_csv('manta_summary.csv')

##


from matplotlib import colors

fig=plt.figure(1)
fig.set_size_inches([7,7],forward=True)
fig.clf()
fig,axs=plt.subplots(2,1,num=1)

spacing=0
for ax,offset,season,marker in zip(axs,[-spacing,spacing],
                                   ['Wet','Dry'],
                                   ['o','v']):
    sel=(manta_per_sample['Season']==season).values & (manta_per_sample['SampleType']=='Trawl')

    df=manta_per_sample[sel]
    
    scat=ax.scatter(df['x'][sel] + offset,
                    df['y'][sel],
                    70,
                    df['part_per_m3'],
                    norm=colors.LogNorm(vmin=1e-3,vmax=100),
                    marker=marker,
                    cmap='CMRmap_r',zorder=2)
    # no data / blank removal took it to zero.
    sel_nan=sel & (df['part_per_m3']<=0.0)
    scat_nan=ax.scatter(df['x'][sel_nan] + offset,
                        df['y'][sel_nan],
                        70,
                        color='none',edgecolors='k',
                        marker=marker,zorder=2)
    
    ax.plot([np.nan],[np.nan],marker=marker,ls='none',ms=7,color='k',label=season)
    
    ax.legend(loc='lower left')
    # ccoll=g.plot_cells(values=-g.cells['z_bed'],
    #                    norm=colors.LogNorm(vmin=2,vmax=4000),
    #                    cmap='gray_r',zorder=0,
    #                    lw=0.5,edgecolor='face',
    #                    ax=ax)

    # g.plot_boundary(lw=0.5,color='k',select_by='cells',ax=ax,zorder=1)
    plot_wkb.plot_polygon(poly,fc='0.85',ec='k',lw=0.5,ax=ax)

    ax.axis('equal')
    zoom=(441437, 613259, 4135548, 4233319)
    ax.axis(zoom)
    ax.xaxis.set_visible(0)
    ax.yaxis.set_visible(0)

cax=fig.add_axes([0.85,axs[1].get_position().ymin,
                  0.025,axs[0].get_position().ymax - axs[1].get_position().ymin])
plt.colorbar(scat,label='Particles / m$^3$',cax=cax)
fig.subplots_adjust(left=0.02,right=0.83)

fig.savefig('manta-blank_sub-particle_counts.png',dpi=150)

##

