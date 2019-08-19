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
manta_per_sample['area_m2']=1e6*grp['Area'].first() # comes in km^2

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
        manta_per_sample.loc[non_dupe,'count_preblank'] += manta_per_sample.loc[v,'count_preblank']
        manta_per_sample.loc[non_dupe,'volume_l'] += manta_per_sample.loc[v,'volume_l']
        manta_per_sample.loc[non_dupe,'area_m2'] += manta_per_sample.loc[v,'area_m2']
        to_delete.append(v)
manta_per_sample.drop(to_delete,inplace=True)

##
manta_per_sample['part_per_m3']= 1000*manta_per_sample['count'] / manta_per_sample['volume_l']
manta_per_sample['part_per_m2']= manta_per_sample['count'] / manta_per_sample['area_m2']

ll=manta_per_sample.loc[: , ['lon','lat']].values
xy=np.nan*np.ones( (len(manta_per_sample),2) )

valid=np.isfinite(ll[:,0])
xy[valid,:]=proj_utils.mapper('WGS84','EPSG:26910')(ll[valid])

manta_per_sample['x']=xy[:,0]
manta_per_sample['y']=xy[:,1]

manta_per_sample.to_csv('manta_summary.csv')

##

plt.figure(1).clf()
plt.plot( manta_per_sample.part_per_m3,
          manta_per_sample.part_per_m2,
          'g.')
plt.xlabel('per m$^3$')
plt.ylabel('per m$^2$')

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

# And the sediment data:
sed=plastic_data.sediment_df

# There is a Mass column, which appears to be the mass of sediment analyzed per SampleID.
#
# StationCode:
#array(['CB10', 'CB15', 'CB32', 'CB37', 'LSB02', 'LSB04', 'LSB06', 'SB051',
#       'SB056', 'SB074', 'SOSL16', 'SOSL40', 'SPB128', 'SPB15', 'SUB52',
#       'SUB53', 'TB101', 'TB102', 'LABQA', 'CB001S', 'SB002S'],

# presumably CB: Central Bay
#            LSB: Lower South Bay
#            SB: South Bay
#            SOSL: Southern Sloughs called out as Margins - Stormwater and Wastewater
#            SPB: San Pablo Bay
#            SUB: Suisun Bay (called out in Group1b as Margin - Wastewater)
#            TB: Tomales Bay
# Group1 includes Ambient, Margins, Reference.
# SOSL falls under Group1B "Margin - Stormwater and Wastewater"

# Overall, I think Group2 is the way to go.

# regardless, have to process them by SampleID first
grp=sed.groupby('SampleID')

sed_samples=pd.DataFrame()

sed_samples['count']=grp.size()
for fld in ['Mass','StationCode','Group1','Group1b','Group2']:
    sed_samples[fld]=grp[fld].first()

sed_samples['lat']=grp['ActualLatitude'].first()
sed_samples['lon']=grp['ActualLongitude'].first()
sed_samples['part_per_mass']=sed_samples['count'] / sed_samples['Mass']

##

grp=sed_samples.groupby('Group2')
sed_groups=pd.DataFrame()
sed_groups['part_per_mass']=grp['part_per_mass'].mean()
# sed_groups['group2']=grp['Group2'].first()
sed_groups['total_particles']=grp['count'].sum()
sed_groups['total_mass']=grp['Mass'].sum()
sed_groups['agg_part_per_mass']=sed_groups['total_particles']/sed_groups['total_mass']

sed_locs=pd.read_csv('sed_loc.csv').set_index('group2')

sed_groups=pd.merge(sed_groups,sed_locs,left_index=True,right_index=True)

## 

sed_groups.to_csv('sed_data_grouped.csv')

## 
# group by embayment?
# do ambient need to be split out from margin?
#  CB ambient 1.76 part/gram
#  SB ambient 2.70 part/gram
# CB margins: more like 3.8
# SB margins: 0.8 to 2.6

# overall, the story from the data looks like:
#   Southern Sloughs and LSB have the highest part/g.
# South Bay is not that high
# Central Bay, San Pablo Bay a bit higher than South Bay (?)
# Suisun more or less like South Bay.
# Tomales lower than most of those.

fig=plt.figure(1)
fig.clf()
ax=fig.add_subplot(1,1,1)
g.plot_cells(values=np.log10(-g.cells['z_bed'].clip(-np.inf,-0.1)),cmap='jet',ax=ax,
             alpha=0.2)
ax.axis('equal')

for idx,row in sed_locs.iterrows():
    xy=sed_plot_locs[k]
    ax.text(row['x'],row['y'],row['code'])
ax.plot(sed_locs['x'],sed_locs['y'],'ko')

ll=np.c_[sed_samples['lon'],
         sed_samples['lat']]
xy=proj_utils.mapper('WGS84','EPSG:26910')(ll)
sed_samples['x']=xy[:,0]
sed_samples['y']=xy[:,1]

for idx,s in sed_samples.iterrows():
    ll=[s['lon'],s['lat']]
    if not np.isfinite(ll[0]):
        print(f"Sample {idx} has no lat/lon")
        continue
    xy=ll2utm(ll)
    ax.text(xy[0],xy[1],f"{s['Group2']}: {idx}",fontsize=10)
    
ax.plot(sed_samples['x'],sed_samples['y'],'go',ms=4)
