"""
v05: try to refactor, process stormwater and effluent at the same time.
process each category independently, sum at the end.
"""
import os
import re

import six
from stompy import utils
import xarray as xr
import numpy as np
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO)
logging.root.setLevel(logging.INFO)

try:
    cwd=os.path.dirname(__file__)
except NameError:
    cwd="." # assume interactive use from same directory
    
utils.path(os.path.join(cwd, "../field_data"))

import plastic_data

##

version='v05'

# Assemble two primary dataframes
# part_df: all particles
# sample_df: all samples, with a many:1 mapping to stations

# trim down the source dataframes before merging
slim_fields=['SampleID','Category_Final','StationCode','pathway','w_s','field_sample_p']
storm_slim=plastic_data.storm_df.loc[:,slim_fields].copy()
waste_slim=plastic_data.effluent_df.loc[:,slim_fields].copy()

# For stormwater, omit the "Bottle Blank" per discussion with Diana
storm_slim=storm_slim[storm_slim['SampleID']!='Bottle Blank']

# All stormwater field samples get combined for a single stormwater distribution
# Can I do this later?  If I do it here, I'll have to manually fix the number
# of samples later on
# storm_slim.loc[ storm_slim.field_sample_p,'StationCode']='storm'

# All of the per-particle data
particles=pd.concat( [storm_slim,waste_slim] )

# There are some 'nan' category particles -- just a handful, but don't count them.
categories=np.array(['Fiber', 'Fiber Bundle', 'Film', 'Foam', 'Fragment', 'Sphere'] )

## 

# Per-sample data:
# For each SampleID in particles, and potentially additional SampleIDs that had no
# particles, I want to know the sample volume.

# These numbers from table 1 in Stormwater Chapter 05-10-2019.docx
# StationCode seems to be just a more descriptive name for what is in
# Watershed
# The keys here match up with StationCode
storm_sample_volumes_l={
    'Line 12F below PG&E station':25,
    'Line 12J at mouth to 12K':67,
    'Refugio Ck at Tsushima St':54,
    'Colma Ck at S. Linden Blvd':197,
    'Rodeo Creek at Seacliff Ct. Pedestrian Br.':63,
    'Line 12K at Coliseum Entrance':295,
    'Guadalupe River at Highway 101':138,
    'Line12MColWay':2*68, # 2 because there is a field dupe with the same volume
    'Line12AShell':115,
    'Meeker Slough':67,
    'FIELDQA':0,
    'MMP-Storm-SB-SM':114,
    'Coyote Creek':114,
    'LABQA':0
}
storm_blank_types=['FIELDQA','LABQA']

# The sample volumes, from Wastewater Chapter 05-09-2019.docx
# each entry maps a regular expression to the sample volume.
# all sample_ids matching the regex should be combined.
waste_sample_volume_patts={
    '20170907.*-CCC?SD-.*':3244,
    '20171206.*-CCCSD-.*':916,
    
    '20170831.*-EBDA-.*':3914,
    '20170926.*-EBDA-.*':7885,
    
    '20170822.*-EBMUD-.*':3483, # looks reasonable, SampleTime=A
    '2017(0926|1020).*-EBMUD-.*':3405, # looks reasonable, SampleTime==B -- combines the split 21h/3h samples

    '20170823.*-FSSD-.*':5954, 
    '20170907.*-FSSD-.*':8150,
    
    'LabBlank.*':0,

    # collected on 2 days:
    # 7/20: included duplicate
    #       A,B: 355A,125A ==> 'both real'
    #            355B,125B ==> DUPLICATE
    # Attached a y-splitter, two flow meters, so 7/20 A vs B should have different
    # volumes, and DL verified it's just luck that 9887 came up twice.
    # 8/01: regular sample
    '20170720.*PA-.*B':9887, # Good - SampleTime=='A'
    '20170801.*-PA-.*': 9887, # good- SampleTime=='B' verified this is the same #, by luck

    # some confusion as to whether this is a field dupe, or if 7/20 125B/355B is officially
    # the field dupe.  But there's no difference in the collection, so just use both.
    '20170720.*-PA-.*A':12313, # reported 12313, field DUPE - SampleTime=='FDUP'

    '20171106.*-SFPUC-\d+$':2859, # ignore -Blank
    '20171107.*-SFPUC-\d+$':2263, # ignore -Blank

    '20171107.*-SFPUC-.*-Blank':0,
    
    '20170810.*-SJ-.*':7586,
    '20170919.*-SJ-.*':4868,
    
    '20170919.*-SUNN-.*':1440,
    '20171017.*-SUNN-.*':2880,
}

grp=particles.groupby(['StationCode','SampleID','pathway','field_sample_p'])
recs=[]
for (station_code,sample_id,pathway,field_sample_p) in grp.groups:
    recs.append(dict(station_code=station_code,
                     sub_sample_id=sample_id,
                     pathway=pathway,
                     field_sample_p=field_sample_p) )
sub_samples=pd.DataFrame(recs)

full_sample_ids=[] # when multiple sampleIDs correspond to different fractions of the same source water.

for idx,row in sub_samples.iterrows():
    if row['pathway']=='stormwater':
        if row['field_sample_p']:
            # had been this, but it discards the number of samples
            # full_sample_ids.append('stormwater_field')
            full_sample_ids.append(row['station_code'])
        else:
            full_sample_ids.append(row['sub_sample_id'])
    elif row['pathway']=='effluent':
        matches= [ patt
                   for patt in waste_sample_volume_patts.keys()
                   if re.match(patt,row['sub_sample_id']) is not None ]
        assert len(matches)==1
        # volumes_l.append( waste_sample_volume_patts[matches[0]] )
        full_sample_ids.append(matches[0])
    else:
        assert False
sub_samples['full_sample']=full_sample_ids

# len(sub_samples['full_sample'].unique()) => 25
# that's 1 stormwater combined sample
# 17 wastewater samples
# 7 blanks

# 
full_samples=sub_samples.groupby(['pathway','station_code','full_sample']).size().reset_index()
full_samples=full_samples.set_index('full_sample')
full_samples.rename(columns={0:'n_sub_samples'},inplace=True)

# sample volume is a property of a full_sample, rather than a sub_sample
full_samples['volume_l']=np.nan
# This is where the stormwater samples are all directed to a single source
full_samples.loc[full_samples.pathway=='stormwater', 'station_code']='stormwater'

for k in storm_sample_volumes_l:
    full_samples.loc[k,'volume_l']=storm_sample_volumes_l[k]
        
for patt in waste_sample_volume_patts:
    vol=waste_sample_volume_patts[patt]
    if vol!=0.0:
        full_samples.loc[patt,'volume_l']=vol
    # blank samples are left with NaN volume.

# And 'station' denotes how full_samples map to outputs, 

# TODO: make sure the handling of multiple sub-samples for wastewater
#   is correct.  There are currently 43 "samples" listed for effluent,
#   but only 19 sample volumes: 2 blanks, 7 sites with 2 each, 1 site with 3 

# Calculate blank rates per category and pathway
blank_recs=[]

for pathway in ['effluent','stormwater']:
    # How many blank were collected?
    sub_sample_sel=(  (sub_samples.pathway==pathway)
                      & (~sub_samples['field_sample_p']) )
    n_blank_samples=sub_sample_sel.sum()
        
    for category in categories:
        part_sel=(particles.pathway==pathway)&(~particles['field_sample_p'])&(particles['Category_Final']==category)
        n_blank_particles=part_sel.sum()
        rec=dict(pathway=pathway,category=category,n_particle=n_blank_particles,n_sample=n_blank_samples)
        blank_recs.append(rec)
blank_rates=pd.DataFrame(blank_recs).set_index(['pathway','category'])
blank_rates['rate']=blank_rates['n_particle']/blank_rates['n_sample']

# With the stormwater numbers omitting Bottle Blank, this now matches the
# previous numbers

# Specify what w_s values will be used for the groups in the PTM
w_s_centers=np.array([0.05,0.005,0.0005,
                      0,
                      -0.05,-0.005,-0.0005])
w_s_centers.sort()

ds=xr.Dataset.from_dataframe(blank_rates).rename({'rate':'blank_rate',
                                                  'n_sample':'n_blank_samples',
                                                  'n_particle':'n_blank_particles'})

ds['w_s']=('w_s',),w_s_centers
ds.w_s.attrs['name']="Center of bin, settling velocity, positive down, m/s"

# source names from samples, omitting the blank samples, 
ds['source']=('source',),full_samples[full_samples['volume_l'].notnull()]['station_code'].unique()

# For each station this volume gets us from particle count to concentration
ds['total_volume']=('source',), np.nan*np.ones(ds.dims['source'],np.float64)
# this number of samples tells how to derate those concentrations
ds['n_samples']=('source',), np.nan*np.ones(ds.dims['source'],np.float64)

# Iterate over samples to get total_volume and n_samples, omitting
# blanks
grp=full_samples[full_samples.volume_l.notnull()].groupby('station_code')
df_grp=pd.DataFrame()
df_grp['volume_l']=grp['volume_l'].sum()
df_grp['count_full']=grp.size()
df_grp['count_sub']=grp['n_sub_samples'].sum()

# Do these numbers make sense? Seem to.

#                 volume_l  count_full  count_sub
#   station_code                                 
#   CCCSD           4160.0           2        4.0 YES - that's wet/dry x {125,355}
#   EBDA           11799.0           2        4.0 YES - same
#   EBMUD           6888.0           2        6.0 YES - EBMUD had one sample that spanned 2 days. 
#   FSSD           14104.0           2        4.0 YES
#   PA             32087.0           3        6.0 YES - A,B, field dupe.
#   SFPUC           5122.0           2        4.0 YES (properly drops the blanks)
#   SJ             12454.0           2        4.0 YES
#   SUNN            4320.0           2        4.0 YES
#   stormwater      1317.0          12       15.0 YES (see below)

# Stormwater: there are 15 sub-samples.  The 3 extras are from
# 'Line 12J at mouth to 12K' has the two size classes split out
# 'Line12MColWay' has a field dupe, and also has a weird, very small count, add'l sample.

derate_per_full_sample=False

for idx,rec in df_grp.iterrows():
    ds.total_volume.loc[idx]=rec['volume_l']
    if derate_per_full_sample:
        ds.n_samples.loc[idx]=rec['count_full']
    else:
        ds.n_samples.loc[idx]=rec['count_sub']

# n_samples is tricky here.
# It is only used to figure out how much to derate the observations.
# The basic question is whether to treat a single sample that was processed in
# two sieves as being susceptible to 1x or 2x the blank rates.
        
## 
# counts where w_s is known
ds['count_w_s']=('source','category','w_s'),np.zeros( (ds.dims['source'],ds.dims['category'],ds.dims['w_s']), np.int32)
# and where it isn't known
ds['count_no_w_s']=('source','category'),np.zeros( (ds.dims['source'],ds.dims['category']), np.int32)


particles['source']=particles['StationCode']
particles.loc[particles.pathway=='stormwater','source']='stormwater'

particles_with_w_s=particles[ particles['Category_Final'].notnull()
                              & particles['field_sample_p']
                              & particles['w_s'].notnull() ].copy()
particles_no_w_s=particles[ particles['Category_Final'].notnull()
                            & particles['field_sample_p']
                            & particles['w_s'].isnull() ]

particles_with_w_s['w_s_binned']=utils.nearest_val(ds.w_s.values,
                                                   particles_with_w_s['w_s'])

count_w_s=particles_with_w_s.groupby(['source','Category_Final','w_s_binned']).size()

for (source,cat,w_s_bin),row in count_w_s.iteritems():
    ds['count_w_s'].loc[source,cat,w_s_bin] = row

df_count_no_w_s=particles_no_w_s.groupby(['source','Category_Final']).size()

for (source,cat),row in df_count_no_w_s.iteritems():
    ds['count_no_w_s'].loc[source,cat] = row

## 

# Scale up counts and get raw concentration
count_with_w_s=ds.count_w_s.sum(dim='w_s') # (source,category)
count_no_w_s=ds['count_no_w_s']
count_total=count_no_w_s+count_with_w_s
# w_s_scale_upA=count_total/count_with_w_s # (source,category)
w_s_scale_up=xr.DataArray(np.where( count_total>0, count_total/count_with_w_s, 0.0),
                           dims=['source','category'])

# there are many cases where a source,category combination has no samples
# with w_s, but luckily there are not cases where a pathway,category combination
# has only samples without w_s.  So the scaling is always valid.

# How does this compare to previous conc_raw?
ds['conc_raw'] = ds.count_w_s * w_s_scale_up / ds.total_volume

if 0: # print values for comparison to previous
    print("------------- NEW CODE v05 --------------")
    print(ds.conc_raw.sum(dim='category').to_dataframe().unstack('w_s'))

    loads_v04=xr.open_dataset('plastic_loads-7classes-v04None.nc')
    print("------------- CODE v04 --------------")
    print(loads_v04.conc_raw.to_dataframe().unstack('w_s'))

    print("------------- NEW CODE v05 --------------")
    print(ds.conc_raw.sum(dim='category').sum(dim='w_s').to_dataframe())

    print("------------- CODE v04 --------------")
    print(loads_v04.conc_raw.sum(dim='w_s').to_dataframe())

    # Overall-- total counts in stormwater goes down, matching the small
    # increase in sample volume.
    # WWTP counts match exactly.
    # distribution over w_s changes a good bit, overall shifting towards more
    # sinking, due to the change in how w_s samples are scaled up.  Quite
    # different, but no indication of a bug.

sps=[]
for s in ds.source.values:
    if s=='stormwater':
        sps.append('stormwater')
    else:
        sps.append('effluent')
        
ds['source_pathway']=('source',), sps

# broadcast to source,category
blank_rates=ds['blank_rate'].loc[ ds.source_pathway ]

# Maybe?
# for example conc_raw says we sampled 100l from EBDA, found 10 passive fibers =>
# conc_raw[ebda,fiber,0.0] = 0.1 particle / l
# across all w_s, ebda had 80 fibers...
# the average wastewater blank found 20 fibers => blank_rates[ebda,fiber]=20
# the ebda data came from 2 samples: n_samples[ebda]=2
# so across all w_s, the blank concentration is 20 * 2 / 100
ds['conc_blank']= blank_rates*ds.n_samples/ds.total_volume
# Verified in dev that these are mostly positive.

ds['blank_derate']=(ds.conc_raw.sum(dim='w_s')-ds.conc_blank)/ds.conc_raw.sum('w_s')
# when a source,category combination has no field particles and no blank particles,
# blank_derate is nan. give it the benefit of the doubt and call it 1.0.
ds['blank_derate'].values=np.nan_to_num(ds['blank_derate'].values.clip(0),nan=1.0)
ds['conc']=ds.conc_raw * ds.blank_derate

## The rest is cleanup, tidying for output
# Drop the nan source
source_sel=np.array( [isinstance(s,str) for s in ds.source.values] )
ds=ds.isel(source=source_sel)
# cast the objects to unicode string
for obj_fld in ['pathway','category','source']:
    ds[obj_fld].values = ds[obj_fld].values.astype('U')

# add some metadata to help with poor memory...
ds.conc.attrs['units']='particle l-1'
ds.conc_raw.attrs['units']='particle l-1'
ds.conc.attrs['description']="Particle concentration adjusted for blank contamination"
ds.blank_rate.attrs['units']='particles per blank sample'
ds.total_volume.attrs['units']='l'

# Print a nice table:
df_out=ds.to_dataframe().unstack()
print(str(df_out))

out_fn=f'plastic_loads-7classes-{version}.nc'
os.path.exists(out_fn) and os.unlink(out_fn)
ds.to_netcdf(out_fn)

# df_out.to_excel(out_fn.replace('.nc','.xlsx'))
