import os
import six
from stompy import utils
import xarray as xr
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)
logging.root.setLevel(logging.INFO)

try:
    cwd=os.path.dirname(__file__)
except NameError:
    cwd="." # assume interactive use from same directory
    
utils.path(os.path.join(cwd, "../field_data"))

import plastic_data

six.moves.reload_module(plastic_data) # DEV
##

#version='v03-nofiber'
version='v03'

## 

def maybe_remove_fiber(df):
    if 'nofiber' in version:
        slim=df[ df['Category_Final']!='Fiber' ]
        print(f"Removing fibers: {len(df)} => {len(slim)} particles")
        return slim
    else:
        return df

# Specify what w_s values will be used for the groups in the PTM
w_s_centers=np.array([0.05,0.005,0.0005,
                      0,
                      -0.05,-0.005,-0.0005])

w_s_centers.sort()

ds=xr.Dataset()

ds['w_s']=('w_s',),w_s_centers
ds.w_s.attrs['name']="Center of bin, settling velocity, positive down, m/s"

##

ds['source']=('source',), ['stormwater',
                           'CCCSD',
                           'EBDA',
                           'EBMUD',
                           'FSSD',
                           'PA',
                           'SFPUC',
                           'SJ',
                           'SUNN' ]

ds['conc']=('w_s','source'),np.nan*np.ones( (ds.dims['w_s'],ds.dims['source']), np.float64)
ds['conc_raw']=('w_s','source'),np.nan*np.ones( (ds.dims['w_s'],ds.dims['source']), np.float64)


## 
df=maybe_remove_fiber(plastic_data.storm_df)


# Goal is to come up with a table of
# particles/L binned by settling velocity.

# Limit to field samples

# StationCode seems to be just a more descriptive name for what is in
# Watershed

# These numbers from tabel 1 in Stormwater Chapter 05-10-2019.docx
sample_volume_l={
    'Line 12F below PG&E station':25,
    'Line 12J at mouth to 12K':67,
    'Refugio Ck at Tsushima St':54,
    'Colma Ck at S. Linden Blvd':197,
    'Rodeo Creek at Seacliff Ct. Pedestrian Br.':63,
    'Line 12K at Coliseum Entrance':295,
    'Guadalupe River at Highway 101':138,
    'Line12MColWay':68,
    'Line12AShell':115,
    'Meeker Slough':67,
    'FIELDQA':0,
    'MMP-Storm-SB-SM':114,
    'Coyote Creek':114,
    'LABQA':0
}

# 1317 l
total_sample_volume_l=np.sum( list(sample_volume_l.values()))

# sel=(df['field_sample_p'] & np.isfinite(df['w_s'])).values
# had been limiting this to particles with w_s, but that was too early
field=df[df['field_sample_p']].copy()

# 1072 left, out of 12,525 original.
##

# Adjust particle weights by blanks.
storm_blank_types=['FIELDQA','LABQA']

blanks=df[ (df['SampleID']!='Bottle Blank') & df['SampleType_AW'].isin(storm_blank_types) ].copy() 

n_blank_samples=len(blanks['SampleID'].unique())
n_field_samples=len(field['SampleID'].unique())

print("--- Stormwater Blanks ---")
if 'LABQA' in storm_blank_types:
    print(f"  Considering Field & Lab Blanks, {n_blank_samples:3} distinct blank samples")
else:
    # 2
    print(f"  Considering Field Blanks only, {n_blank_samples:3} distinct blank samples")
    
# 15
print(f"                                 {n_field_samples:3} distinct field samples")

# Assume that each sample accumulates the same counts of erroneous particles.
blank_category_per_sample=blanks.groupby('Category_Final').size() / n_blank_samples
print("Per-blank sample, per-category counts")
print(blank_category_per_sample)
print()

# All particles start with weight of 1.0
field['derated']=1.0

for cat,blank_freq in blank_category_per_sample.iteritems():
    print(f"Stormwater, category={cat}.  Per-blank count {blank_freq}")
    cat_sel=field['Category_Final']==cat    
    n_field_per_cat=cat_sel.sum()
    print(f"    field count over {n_field_samples:3} samples: {n_field_per_cat}, {n_field_per_cat/n_field_samples:.2f} per sample")
    real_derate=(n_field_per_cat - n_field_samples*blank_freq) / n_field_per_cat
    derate=max(0.0,real_derate)
    print(f"    de-rating: {derate:.3f} ({real_derate:.3f})")
    print()
    
    field.loc[ cat_sel, 'derated']=derate
##

# simplest approach:
#   Assume that all stations are sampling the same water (i.e. a single distibution
#   over all sites).  Then it doesn't matter at which station a particle was counted.

bin_conc=np.zeros(len(w_s_centers),np.float64)
bin_conc_raw=np.zeros(len(w_s_centers),np.float64)
# field includes w/ and w/o w_s here --
field_valid=np.isfinite(field['w_s'].values)
# but bin_choice now reflects just the field samples with a w_s
bin_choice=utils.nearest(w_s_centers,field['w_s'].values[field_valid])

# how many non-QA samples?  Both StationCode='LABQA' and 'FIELDQA' have
# SampleGross_AW='QA'.
# => 12362
fraction_with_w_s = field_valid.sum() / len(field)

for bin_idx,particles in utils.enumerate_groups(bin_choice):
    # Here the particles are de-rated based on the per-category blank counts
    # Lots going on here:
    #                   the derated particle 'weight' for *all* particles
    #                                           just the ones with w_s
    #                                                       of those, just this bin of w_s
    #                                                                           to a concentration
    #                                                                                                   scale up to account for
    #                                                                                                   how many particles don't have w_s
    bin_conc[bin_idx] = field['derated'].values[field_valid][particles].sum() / total_sample_volume_l / fraction_with_w_s
    # But do the non-de-rated calc as well:
    bin_conc_raw[bin_idx] = len(particles) / total_sample_volume_l / fraction_with_w_s

# So adjust the bin concentrations to reflect the total number of particles

for w_s, conc, conc_raw in zip(w_s_centers,bin_conc,bin_conc_raw):
    print(f" w_s ~ {w_s: .5f}m/s    conc ~ {conc:.3f} particles/l, raw ~ {conc_raw:.3f}")
print(f"  Total                      {bin_conc.sum():.3f} particles/l, raw ~ {bin_conc_raw.sum():.3f}")

##

# fill in the dataset
ds['conc'].sel(source='stormwater').values[:]= bin_conc
ds['conc_raw'].sel(source='stormwater').values[:]= bin_conc_raw

##
# Output:
#
#  w_s ~ -0.05000m/s    conc ~ 0.245 particles/l, raw ~ 0.245
#  w_s ~ -0.00500m/s    conc ~ 0.514 particles/l, raw ~ 0.517
#  w_s ~ -0.00050m/s    conc ~ 0.543 particles/l, raw ~ 0.569
#  w_s ~  0.00000m/s    conc ~ 0.016 particles/l, raw ~ 0.018
#  w_s ~  0.00050m/s    conc ~ 1.097 particles/l, raw ~ 1.217
#  w_s ~  0.00500m/s    conc ~ 6.195 particles/l, raw ~ 6.252
#  w_s ~  0.05000m/s    conc ~ 0.568 particles/l, raw ~ 0.569
#   Total                      9.179 particles/l, raw ~ 9.386

##

# And wastewater:
df=maybe_remove_fiber(plastic_data.effluent_df)

# StationCode: 'CCCSD', 'EBDA', 'EBMUD', 'FSSD', 'LABQA', 'PA', 'SFPUC', 'SJ', 'SUNN'
# SampleID: includes date, station and sieve size.
# SampleTypeCode: 'Field', 'LabBlank', 'FieldBlank'
#    - select on this

# The sample volumes, from Wastewater Chapter 05-09-2019.docx
# each entry maps a regular expression to the sample volume.
# all sample_ids matching the regex should be combined.
import re
sample_volume_patts={
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

# Make sure that all of the sample_IDs line up with exactly one pattern
for sample_id in df.SampleID.unique():
    matches= [ patt for patt in sample_volume_patts.keys() if  re.match(patt,sample_id) is not None ]
    if len(matches)!=1:
        print(f"Got {len(matches)} matches for sample id {sample_id}")
        print(matches)
        raise Exception("Matches are off")
##

# Summarize the blank counts by category
include_lab_blanks=True

if include_lab_blanks:
    blanks=df[ ~df['field_sample_p'] ].copy() # use field and lab blanks
else:
    blanks=df[ (df['StationCode']!='LABQA') & (~df['field_sample_p']) ].copy() # only use field blanks

field=df[df['field_sample_p']].copy()
    
n_blank_samples=len(blanks['SampleID'].unique())
n_field_samples=len(field['SampleID'].unique())

print("--- Effluent Blanks ---")
if include_lab_blanks:
    print(f"  Considering Field & Lab Blanks, {n_blank_samples:3} distinct blank samples")
else:
    print(f"  Considering Field Blanks only, {n_blank_samples:3} distinct blank samples")
    
# 15
print(f"                                 {n_field_samples:3} distinct field samples")

# Assume that each sample accumulates the same counts of erroneous particles.
blank_category_per_sample=blanks.groupby('Category_Final').size() / n_blank_samples
print("Per-blank sample, per-category counts")
print(blank_category_per_sample)
print()

##

# Derate
# All particles start with weight of 1.0
# The derating has to be done per-station.
field['derated']=1.0

stationcodes=field.StationCode.unique()

for station in stationcodes:
    field_station=field[(field.StationCode==station)].copy()
    total_count=len(field_station)
    field_valid=np.isfinite(field_station.w_s)
    bins=utils.nearest(w_s_centers,field_station[field_valid].w_s.values)
    w_s_count=field_valid.sum()

    # what was the total volume sampled for this station?
    sample_volume_l=0
    sample_ids=field_station.SampleID.unique()
    for patt in sample_volume_patts.keys(): # does this volume go with this station?
        for samp_id in sample_ids: # yes, iff 1 or more sample_ids match
             if re.match(patt,samp_id) is not None:
                 sample_volume_l+=sample_volume_patts[patt]
                 break
    print(f"{station:20} total sample volume {sample_volume_l:7} l")

    # Derate
    n_field_samples=len(sample_ids)
    
    for cat,blank_freq in blank_category_per_sample.iteritems():
        print(f"Effluent, category={cat}.  Per-blank count {blank_freq:.3f}")
        cat_sel=field_station['Category_Final']==cat    
        n_field_per_cat=cat_sel.sum() # total number of field samples at this station of this category
        print(f"    At {station}, field count over {n_field_samples:3} samples: {n_field_per_cat}, {n_field_per_cat/n_field_samples:.2f} per sample")
        if n_field_per_cat>0:
            real_derate=(n_field_per_cat - n_field_samples*blank_freq) / n_field_per_cat
        else:
            real_derate=0.0
        derate=max(0.0,real_derate)
        print(f"    de-rating: {derate:.3f} ({real_derate:.3f})")
        print()
        field_station.loc[ cat_sel, 'derated']=derate

    bin_conc_raw=np.zeros(len(w_s_centers))
    bin_conc=np.zeros(len(w_s_centers))

    # Adjustment factor for how many samples have w_s vs total counted
    adjust=total_count / w_s_count
    assert sample_volume_l>0.0
    for bin_idx,particles in utils.enumerate_groups(bins):
        bin_conc_raw[bin_idx] = len(particles) * adjust/sample_volume_l
        bin_conc[bin_idx]=field_station['derated'].values[field_valid][particles].sum() * adjust/sample_volume_l 

    ds.conc_raw.sel(source=station).values[:] = bin_conc_raw
    ds.conc.sel(source=station).values[:] = bin_conc
    assert np.all(np.isfinite(bin_conc))

# StationCode: 'CCCSD', 'EBDA', 'EBMUD', 'FSSD', 'LABQA', 'PA', 'SFPUC', 'SJ', 'SUNN'

ds.conc.attrs['units']='particle l-1'

# Print a nice table:
df_out=ds.to_dataframe().unstack()
print(str(df_out))

##


##

out_fn=f'plastic_loads-7classes-{version}.nc'
os.path.exists(out_fn) and os.unlink(out_fn)
ds.to_netcdf(out_fn)

##

##
df_out.to_excel(out_fn.replace('.nc','.xlsx'))

##


#                  conc                                                    \
#  w_s          -0.0500   -0.0050   -0.0005    0.0000    0.0005    0.0050   
#  source                                                                   
#  CCCSD       0.017457  0.034913  0.011870  0.000000  0.013965  0.006983   
#  EBDA        0.006716  0.007233  0.007233  0.000000  0.010333  0.004650   
#  EBMUD       0.015244  0.052359  0.011267  0.000000  0.014581  0.011930   
#  FSSD        0.000000  0.001407  0.000469  0.000000  0.003752  0.003049   
#  PA          0.000000  0.000698  0.000349  0.000000  0.003840  0.004538   
#  SFPUC       0.027176  0.090588  0.026170  0.000000  0.026170  0.023150   
#  SJ          0.004212  0.011433  0.002407  0.000903  0.003911  0.003009   
#  SUNN        0.000992  0.000992  0.001984  0.000000  0.011905  0.011905   
#  stormwater  0.245147  0.513660  0.543446  0.016483  1.096585  6.195458   
#  
#                        conc_raw                                          \
#  w_s           0.0500   -0.0500   -0.0050   -0.0005    0.0000    0.0005   
#  source                                                                   
#  CCCSD       0.002793  0.017457  0.034913  0.011870  0.000000  0.013965   
#  EBDA        0.001550  0.006716  0.007233  0.007233  0.000000  0.010333   
#  EBMUD       0.001326  0.015244  0.052359  0.011267  0.000000  0.014581   
#  FSSD        0.000469  0.000000  0.001407  0.000469  0.000000  0.003752   
#  PA          0.000175  0.000000  0.000698  0.000349  0.000000  0.003840   
#  SFPUC       0.001007  0.028026  0.093421  0.026988  0.000000  0.022836   
#  SJ          0.000301  0.004212  0.011433  0.002407  0.000903  0.003911   
#  SUNN        0.000000  0.000992  0.000992  0.001984  0.000000  0.011905   
#  stormwater  0.567930  0.245169  0.516607  0.569143  0.017512  1.217091   
#  
#                                  
#  w_s           0.0050    0.0500  
#  source                          
#  CCCSD       0.006983  0.002793  
#  EBDA        0.004650  0.001550  
#  EBMUD       0.011930  0.001326  
#  FSSD        0.003049  0.000469  
#  PA          0.004538  0.000175  
#  SFPUC       0.014532  0.001038  
#  SJ          0.003009  0.000301  
#  SUNN        0.011905  0.000000  
#  stormwater  6.251819  0.569143  



# Note that these are concentrations, so not loads...
# Interesting conclusion that on the buoyant end, the gap is not nearly so large,
# esp. comparing SFPUC, CCCSD, and EBMUD to stormwater.  Factor of 10, but
# much closer than on the settling end where stormwater is 6 particles/l.
# and sfpuc at -0.005 is within factor of 5 of stormwater
