# Show the field data in relation to
#
from matplotlib import colors
import numpy as np
import pdb

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

src_df=plastic_data.manta_df
manta_no_labqa=src_df[ src_df.StationCode!='LABQA' ]

##

# 2019-08-28: Fix SPB3 / Aug21 starting latitude per Carolynn's email.
manta_volumes=pd.read_excel("../field_data/Manta Volumes-fix20190828.xlsx",sheet_name=0,skiprows=1,
                            na_values=['None','Missing'])
manta_volumes['DATE']=pd.to_datetime(manta_volumes['DATE'])

##

# Rather than joining based on SAMPLE_LOCATION, need to be
# a bit smarter.
# first, there are some rows in manta_volumes that are not useful.
# we don't want blank data in there, and don't want CB7-1
# starts off with 74 rows.
manta_volumes1=manta_volumes[ manta_volumes['SAMPLE LOCATION']!='CB7-1'].copy()
# now 65 rows, matching the report of 65 field samples
# and standardize name of CB7
def vol_station_rename(src,dst):
    manta_volumes1.loc[ manta_volumes1['SAMPLE LOCATION']==src,'SAMPLE LOCATION'] = dst

vol_station_rename('CB7-2','CB7')
vol_station_rename('MBNMS-29','MBNMS29')
vol_station_rename('SFBay9-18','SFBay')

# standardize the blank naming:
blanks=manta_volumes1.TYPE.isin(['Manta - Blank','Manta - BLANK']) 
manta_volumes1.loc[blanks,'TYPE']='Manta - Blank'

print("This should match the report, with 58 regular, 8 blank, 7 dupe")
print(manta_volumes1.groupby('TYPE').size())

##

# rather than join by SAMPLE LOCATION, which is not unique,
# populated SampleIDs on manta_volumes, then join with that.

manta_volumes1['SampleID']=None
manta_volumes1['FibersYN']=None
manta_volumes1['Season']=None

def match_samples(recs):
    station_codes=recs['StationCode'].unique()
    rec_dates=recs['SampleDate'].unique()
    rec_seasons=recs['Season'].unique()
    fibers=recs['FibersYN'].unique()
    
    assert len(station_codes)==1
    assert len(rec_dates)==1
    assert len(fibers)==1
    assert len(rec_seasons)==1
    
    station_code=station_codes[0]
    rec_date=rec_dates[0]
    rec_season=rec_seasons[0]
    rec_fibers=fibers[0]
    
    sample_id=recs['SampleID'].values[0]

    match_loc=manta_volumes1['SAMPLE LOCATION']==station_code
    match_date=manta_volumes1['DATE']==rec_date

    # this could be a dupe.
    if 'DUP' in sample_id:
        match_type=manta_volumes1['TYPE']=='Manta - DUP'
    elif ('Blank' in sample_id) or ('BLANK' in sample_id):
        match_type=manta_volumes1['TYPE']=='Manta - Blank'
    else:
        match_type=manta_volumes1['TYPE']=='Manta'
        
    match = match_loc & match_date & match_type
    
    if match.sum()==0:
        print(f"No match from manta station {station_code} to volumes for date {rec_date}")
        # No match from manta station MBNMS29 to volumes for date 2017-11-17
        # 
        return False
    if match.sum()==1:
        if manta_volumes1.loc[match,'SampleID'].values[0] is not None:
            pdb.set_trace()
        manta_volumes1.loc[match,'SampleID'] = sample_id
        manta_volumes1.loc[match,'FibersYN'] = rec_fibers
        manta_volumes1['Season']=rec_season
        return True
    else:
        print(f"Station {station_code} and date {rec_date} and dup {sample_id} were not unique")
        pdb.set_trace()
        return False

# we're just trying to map sample_ids onto volumes, without regard to fibers/nofibers.
# so use the manta_no_labqa dataframe.
sample_ids=manta_no_labqa.groupby('SampleID').apply(match_samples)
assert np.all(sample_ids.values),"Some sample ids failed to match a volume"

##

# 2 manta volumes are missing a match to records.
# one of those is a blank from SPB3, 2017-11-17. 
# there are 212 particles from SPB3 on that date, but all from a trawl, no field blank.
#SPB3         2017-08-21     836
#             2017-11-17     212
# the report does say that a field blank was collected in wet weather at SPB3, and
# that fibers were not counted in it. 

no_match=manta_volumes1['SampleID'].isnull()
volumes_with_no_recs=manta_volumes1.loc[ no_match, ['SAMPLE LOCATION','DATE','TYPE']]
print(f"Locations in manta_volumes that got no records: ")
print(volumes_with_no_recs)
print("   these are field blanks with no fibers counted, and are assumed to have been only fibers.")

# So I'm going to assume that SPB3 is a field blank that was only fibers,
# and fibers weren't counted, and just manufacture a sampleID.
spb3=manta_volumes1['SAMPLE LOCATION']=='SPB3'
manta_volumes1.loc[ no_match & spb3, 'SampleID' ]='SPB3-Manta-Blank-17Nov17'
manta_volumes1.loc[ no_match & spb3, 'FibersYN' ]='N'

# The other one is CBNMS22 blank from wet weather
# it also is shown as not having fibers counted, and a blank is listed in the report.
cbnms22=manta_volumes1['SAMPLE LOCATION']=='CBNMS22'
manta_volumes1.loc[ no_match & cbnms22, 'SampleID' ]='CBNMS22-Manta-Blank-30Mar2018'
manta_volumes1.loc[ no_match & cbnms22, 'FibersYN' ]='N'

assert manta_volumes1['SampleID'].isnull().sum()==0,"Some manta volumes don't have a sample id"
assert manta_volumes1['FibersYN'].isnull().sum()==0,"Some manta volumes don't have a fibers flag"


assert manta_volumes1[ manta_volumes1.TYPE.isin(['Manta','Manta - DUP']) ]['LAT START'].isnull().sum()==0

# At this point manta_volumes is verified to have exactly one row for every sample, to
# have sample_ids which match the individual particles.
# each has at least a LAT START, though not necessarily a LAT END
# 

## 

df=manta_no_labqa
manta_std=df[ df['FibersYN']=='Y' ].copy()
print(f"Removing samples that didn't count fibers: {len(df)} => {len(manta_std)} particles")
df=manta_no_labqa
manta_nofiber=df[ df['Category_Final']!='Fiber' ].copy()
print(f"Removing fibers: {len(df)} => {len(manta_nofiber)} particles")

## 

# choose a representative point in the middle
manta_volumes1['lat']=0.5*( manta_volumes1['LAT START'] + manta_volumes1['LAT END'])
manta_volumes1['lon']=0.5*( manta_volumes1['LONG START'] + manta_volumes1['LONG END'])

missing=manta_volumes1['lat'].isnull() & manta_volumes1['TYPE'].isin(['Manta','Manta - DUP'])
manta_volumes1.loc[missing,'lat']=manta_volumes1.loc[missing,'LAT START']
manta_volumes1.loc[missing,'lon']=manta_volumes1.loc[missing,'LONG START']

missing2=missing & manta_volumes1['lat'].isnull()
assert missing2.sum()==0

##

grp_std=manta_std.groupby('SampleID')
grp_nofiber=manta_nofiber.groupby('SampleID')

##

manta_per_sample=manta_volumes1.set_index('SampleID').rename(columns={'VOLUME (m^3)':'volume_m3',
                                                                      'AREA (km^2)':'area_km2'})
# SI, please
manta_per_sample['volume_l']=manta_per_sample['volume_m3']*1000.0
manta_per_sample['area_m2'] =manta_per_sample['area_km2']*1e6

##


for fld,grp in [('count_preblank_std',grp_std),
                ('count_preblank_nofiber',grp_nofiber)]:
    # This honors the index on SampleID, leaving some rows NaN
    manta_per_sample[fld] = grp.size()
    # so set those to 0
    manta_per_sample.loc[manta_per_sample[fld].isnull(),fld]=0.0

## 

# counts in each sample per category, and unstack the categories to
# become columns

per_sample_per_cat=manta_no_labqa.groupby(['SampleID','Category_Final']).size().unstack(fill_value=0.0)
manta_per_sample2=pd.merge(manta_per_sample,per_sample_per_cat, left_index=True, right_index=True,
                           how='left')
# if there were no records, then manta_per_sample2 will get a NaN, but we'd like it to get
# 0.
for col in per_sample_per_cat.columns.values:
    print(col)
    missing=manta_per_sample2[col].isnull()
    manta_per_sample2.loc[missing,col]=0.0

##

# This is almost certainly to delete
# # wait -- I'm doing this processing at the category level, so why split
# # out std vs nofiber here??
# suffixes=[]
# for suffix,manta_set in [('_std',manta_std),
#                          ('_nofiber',manta_nofiber)]:
#     per_sample_per_cat=manta_set.groupby(['SampleID','Category_Final']).size().unstack(fill_value=0.0)
#     renames=dict([ (x,x+suffix) for x in per_sample_per_cat.columns.values])
#     per_sample_per_cat=per_sample_per_cat.rename(columns=renames)
#     manta_per_sample2=pd.merge(manta_per_sample2,per_sample_per_cat, left_index=True, right_index=True,
#                                how='left')
#     # if there were no records, then manta_per_sample2 will get a NaN, but we'd like it to get
#     # 0.
#     for col in renames.values():
#         print(col)
#         missing=manta_per_sample2[col].isnull()
#         manta_per_sample2.loc[missing,col]=0.0
#     suffixes.append(suffix)

##

# Becky suggests mean, and separate by category

cats=list(manta_no_labqa.Category_Final.unique())

adj_cats=[]
for cat in cats:
    if cat=='Fiber':
        # blank level for Fibers can only come from blanks where fibers were counted
        blank_p=(manta_per_sample2['TYPE']=='Manta - Blank')&(manta_per_sample2['FibersYN']=='Y')
    else:
        blank_p=(manta_per_sample2['TYPE']=='Manta - Blank')

    n_blanks=blank_p.sum()
    print(f"  {cat} {n_blanks} blank samples") # 8 total, but only 4 with fibers

    mean_blank_per_sample_cat = manta_per_sample2.loc[blank_p,cat].mean()
    print(f"    {cat:15}: average of {mean_blank_per_sample_cat:7.3f} in blanks")

    adj_cat=cat+'_adj' # name of column for adjust, post-blank removal count
    manta_per_sample2[adj_cat] = np.maximum(0, manta_per_sample2[cat]-mean_blank_per_sample_cat)
    adj_cats.append(adj_cat)

##


# This part does have to be done by suffix
fibers=manta_per_sample2.FibersYN=='Y'
manta_per_sample2['count_std'] = manta_per_sample2.loc[:, adj_cats ].values.sum(axis=1)
# but if fibers weren't counted it's not valid
manta_per_sample2.loc[~fibers,'count_std']= np.nan
# or count everything but fibers
non_fiber_cats=[c for c in adj_cats if c!='Fiber_adj']
manta_per_sample2['count_nofiber']= manta_per_sample2.loc[:, non_fiber_cats ].values.sum(axis=1)

##

# Could combine DUP, but maybe more interesting to keep them in there as a way
# to see some variability later.  If this code is re-enabled, it will have to
# be updated with the new structure where std and nofiber are in the same dataframe.
if 0:
    to_delete=[]
    for v in manta_per_sample2.index.values:
        if 'DUP' in v:
            print(v)
            non_dupe=v.replace('-DUP','')

            if non_dupe not in manta_per_sample2.index.values:
                # possible that the dupe counted fibers, but the non-dupe did not
                # and has been dropped.
                continue

            # need to be careful about combining a DUPE where one had fibers counted
            # and the other did not.  It does happen ('GFNMS26-Manta-12Sept2017')
            # the original sample does not have fibers, and the dupe does.
            if (manta_per_sample2.loc[v,'FibersYN']!=manta_per_sample2.loc[non_dupe,'FibersYN']):
                if 'nofiber' in version:
                    # No problem - we don't care about fibers anyway
                    pass
                else:
                    # Whichever sample
                    raise Exception("A sample and a dupe have differing FibersYN, and we *are* counting fibers")

            for to_add in ['count','count_preblank','volume_l','area_m2'] + adj_cats + cats:
                manta_per_sample2.loc[non_dupe,to_add] += manta_per_sample2.loc[v,to_add]
            to_delete.append(v)
    manta_per_sample2.drop(to_delete,inplace=True)

##

# 2019-11-08 Remove blanks. Not sure why they were left in before, but it makes processing more
# annoying.
blanks=manta_per_sample2.TYPE=='Manta - Blank'
manta_per_sample3=manta_per_sample2.loc[ ~blanks, :].copy()

##
for suffix in ['_std','_nofiber']:
    manta_per_sample3['part_per_m3'+suffix]= 1000*manta_per_sample3['count'+suffix] / manta_per_sample3['volume_l']
    manta_per_sample3['part_per_m2'+suffix]= manta_per_sample3['count'+suffix] / manta_per_sample3['area_m2']

    manta_per_sample3['part_per_m3_raw']= 1000*manta_per_sample3['count_preblank'+suffix] / manta_per_sample3['volume_l']
    manta_per_sample3['part_per_m2_raw']= manta_per_sample3['count_preblank'+suffix] / manta_per_sample3['area_m2']

##
ll=manta_per_sample3.loc[: , ['lon','lat']].values
xy=np.nan*np.ones( (len(manta_per_sample3),2) )

valid=np.isfinite(ll[:,0])
xy[valid,:]=proj_utils.mapper('WGS84','EPSG:26910')(ll[valid])

manta_per_sample3['x']=xy[:,0]
manta_per_sample3['y']=xy[:,1]

manta_per_sample3.to_csv(f'manta_summary-v02.csv')


## --------------------------------------------------------------------------------

if 1:
    plt.figure(1).clf()
    plt.plot( manta_per_sample3.part_per_m3,
              manta_per_sample3.part_per_m2,
              'g.')
    plt.xlabel('per m$^3$')
    plt.ylabel('per m$^2$')

    fig=plt.figure(1)
    fig.set_size_inches([7,7],forward=True)
    fig.clf()
    fig,axs=plt.subplots(2,1,num=1)

    spacing=0
    for ax,offset,season,marker in zip(axs,[-spacing,spacing],
                                       ['Wet','Dry'],
                                       ['o','v']):
        sel=(manta_per_sample3['Season']==season).values & (manta_per_sample3['SampleType']=='Trawl')

        df=manta_per_sample3[sel]

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

    fig.savefig(f'manta-blank_sub-particle_counts-{version}.png',dpi=150)

##

# And the sediment data:
sed=maybe_remove_fiber(plastic_data.sediment_df)

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
assert all( grp['field_sample_p'].first() == grp['field_sample_p'].last() )
sed_samples['field_sample_p']=grp['field_sample_p'].first()

# Blank derating:
#  which samples are blanks?
#  Group1==LABQA, or StationCode=LABQA
#   there are two LABQA samples.
#   but I think one of those is bad?
# sediment chapter says 3 lab blanks, one field blank.
# the second Lab Blank had no particles, so I don't have record of it here.
# there is a field blank called 17MMP-S-LSB04-MP-FB

# Diana's stormwater chapter lists the blank counts/morphology: 
# blank_per_cat_per_sample={'Fiber':40.3,
#                           'Fragment':2.0,
#                           'Film':0.3}

# This gets the same values.
blank_samples=sed[ (sed['StationCode']=='LABQA') | (sed['SampleID']=='17MMP-S-LSB04-MP-FB') ]
# LabBlank-1: 18 particles, size fraction > 500um.
# LabBlank-2: 0 particles.  NOTE: manually included in count below
# LabBlank-3: 1 particle, size fraction 125 to 355.
n_blanks=4 # 3 lab blanks, 1 field.
blank_per_cat_per_sample=blank_samples.groupby('Category_Final').size() / n_blanks
blank_per_cat_per_sample.rename('blank_rate',inplace=True)

##

sed['derate']=1.0

#for cat,blank_rate in blank_per_cat_per_sample.iteritems():
sample_cat_counts=sed.groupby(['SampleID','Category_Final']).size().rename('count')
merge_blank=pd.merge(sample_cat_counts,blank_per_cat_per_sample,left_on='Category_Final',right_index=True)
merge_blank['count_adj']=np.maximum(0,merge_blank['count']-merge_blank['blank_rate'])
adj_counts=merge_blank.groupby('SampleID')['count_adj'].sum()
sed_samples['count_adj']=adj_counts

##
sed_samples['part_per_mass_raw']=sed_samples['count'] / sed_samples['Mass']
sed_samples['part_per_mass']=sed_samples['count_adj'] / sed_samples['Mass']

##

# Write out the per-sample data, but limit to field samples
sed_samples[ sed_samples.field_sample_p ].to_csv(f'sed_samples-{version}.csv')


##
grp=sed_samples[ sed_samples['field_sample_p'] ].groupby('Group2')
sed_groups=pd.DataFrame()
sed_groups['part_per_mass']=grp['part_per_mass'].mean()
sed_groups['part_per_mass_raw']=grp['part_per_mass_raw'].mean()
# sed_groups['group2']=grp['Group2'].first()
sed_groups['total_particles']=grp['count'].sum()
sed_groups['total_mass']=grp['Mass'].sum()
sed_groups['agg_part_per_mass']=sed_groups['total_particles']/sed_groups['total_mass']

sed_locs=pd.read_csv('sed_loc.csv').set_index('group2')

sed_groups=pd.merge(sed_groups,sed_locs,left_index=True,right_index=True)

## 

sed_groups.to_csv(f'sed_data_grouped-{version}.csv')

## 
# group by embayment?
# do ambient need to be split out from margin?
# [these numbers are from before the blank adjustment]
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
    ax.text(row['x'],row['y'],row['code'])
ax.plot(sed_locs['x'],sed_locs['y'],'ko')

ll2utm=proj_utils.mapper('WGS84','EPSG:26910')
ll=np.c_[sed_samples['lon'],
         sed_samples['lat']]
xy=ll2utm(ll)
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
