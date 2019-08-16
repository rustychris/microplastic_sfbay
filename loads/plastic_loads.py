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

# Specify what w_s values will be used for the groups in the PTM

# w_s_centers=[0.06,0.02,0.006,0.002,0.0006,
#              0,
#              -0.06,-0.02,-0.006,-0.002,-0.0006]

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


## 
df=plastic_data.storm_df


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

sel=((df['SampleGross_AW']=='Field') & (np.isfinite(df['w_s']))).values

# 1088 left, out of 12,525 original.

# simplest approach:
#   Assume that all stations are sampling the same water (i.e. a single distibution
#   over all sites).  Then it doesn't matter at which station a particle was counted.

bin_conc=np.zeros(len(w_s_centers),np.float64)
bin_choice=utils.nearest(w_s_centers,df.iloc[sel,:]['w_s'].values)

for bin_idx,particles in utils.enumerate_groups(bin_choice):
    bin_conc[bin_idx] = len(particles) / total_sample_volume_l

##

# That only accounts for the particles that had enough information to get a w_s.
# From the chapter, what do we know about particles getting counted vs. particles
# having enough info to get w_s?

# how many non-QA samples?  Both StationCode='LABQA' and 'FIELDQA' have
# SampleGross_AW='QA'.
# => 12362
total_non_qa_particles= (df['SampleGross_AW']=='Field').sum()

fraction_with_w_s=sel.sum() / total_non_qa_particles

# So adjust the bin concentrations to reflect the total number of particles
adj_bin_conc = bin_conc/fraction_with_w_s

for w_s, conc in zip(w_s_centers,adj_bin_conc):
    print(f" w_s ~ {w_s: .5f}m/s    conc ~ {conc:.3f} particles/l")
print(f"  Total                      {adj_bin_conc.sum():.3f} particles/l")

##

# fill in the dataset
ds['conc'].sel(source='stormwater').values[:]= adj_bin_conc


##
# Output:
# 
#  w_s ~ -0.05000m/s    conc ~ 0.250 particles/l
#  w_s ~ -0.00500m/s    conc ~ 0.552 particles/l
#  w_s ~ -0.00050m/s    conc ~ 0.707 particles/l
#  w_s ~  0.00000m/s    conc ~ 0.017 particles/l
#  w_s ~  0.00050m/s    conc ~ 1.147 particles/l
#  w_s ~  0.00500m/s    conc ~ 6.160 particles/l
#  w_s ~  0.05000m/s    conc ~ 0.552 particles/l
#   Total                      9.386 particles/l

##

# And wastewater:
# 
df=plastic_data.effluent_df

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
stationcodes=df.StationCode.unique()

for station in stationcodes:
    if station in ['LABQA']: # station ids we don't care about
        continue
    df1=df[ df.StationCode==station ]
    total_count=len(df1)
    df2=df1[ np.isfinite(df1.w_s) ]
    bins=utils.nearest(w_s_centers,df2.w_s.values)
    w_s_count=len(df2)
    
    # what was the total volume sampled for this station?
    sample_volume_l=0
    sample_ids=df1.SampleID.unique()
    for patt in sample_volume_patts.keys(): # does this volume go with this station?
        for samp_id in sample_ids: # yes, iff 1 or more sample_ids match
             if re.match(patt,samp_id) is not None:
                 sample_volume_l+=sample_volume_patts[patt]
                 break
    print(f"{station:20} total sample volume {sample_volume_l:7} l")

    bin_conc=np.zeros(len(w_s_centers))
    for bin_idx,particles in utils.enumerate_groups(bins):
        bin_conc[bin_idx] = len(particles) / sample_volume_l

    # adjust for the fraction of particles with w_s, and that these 8 POTWs
    # are 70% of the overall discharge
    # 2019-08-16: no longer adjust for the 70% figure here.
    #adj_bin_conc = (1.00/0.70) * bin_conc * total_count / w_s_count
    adj_bin_conc = bin_conc * total_count / w_s_count
    ds.conc.sel(source=station).values[:] = adj_bin_conc

# StationCode: 'CCCSD', 'EBDA', 'EBMUD', 'FSSD', 'LABQA', 'PA', 'SFPUC', 'SJ', 'SUNN'

# ds depends on w_s_centers, so don't write it out

ds.conc.attrs['units']='particle l-1'

##

out_fn='plastic_loads-7classes-v02.nc'
os.path.exists(out_fn) and os.unlink(out_fn)
ds.to_netcdf(out_fn)

##

# Print a nice table:
df_out=ds.to_dataframe().unstack()
print(str(df_out))

df_out.to_excel(out_fn.replace('.nc','.xlsx'))

# This is slightly out-dated -- see the .xlsx for a better
# view.

#OLD ------------------- End point ---------------
#OLD                 conc                                                            
#OLD w_s          -0.0500   -0.0050   -0.0005    0.0000    0.0005    0.0050    0.0500
#OLD source                                                                          
#OLD CCCSD       0.025709  0.049513  0.019043  0.000000  0.018091  0.009522  0.003809
#OLD EBDA        0.009796  0.010496  0.011895  0.000000  0.013994  0.005598  0.002099
#OLD EBMUD       0.021121  0.077138  0.016530  0.000000  0.020203  0.015611  0.001837
#OLD FSSD        0.000000  0.002287  0.001307  0.000000  0.005226  0.003593  0.000653
#OLD PA          0.000000  0.001600  0.001143  0.000000  0.005028  0.005714  0.000229
#OLD SFPUC       0.038278  0.133972  0.038278  0.000000  0.035544  0.030075  0.001367
#OLD SJ          0.005753  0.017670  0.003287  0.001233  0.005342  0.003698  0.000411
#OLD SUNN        0.001167  0.005836  0.004669  0.000000  0.014006  0.014006  0.000000
#OLD stormwater  0.243128  0.599137  0.642553  0.017366  1.172225  6.156353  0.555722

# Note that these are concentrations, so not loads...
# Interesting conclusion that on the buoyant end, the gap is not nearly so large,
# esp. comparing SFPUC, CCCSD, and EBMUD to stormwater.  Factor of 10, but
# much closer than on the settling end where stormwater is 6 particles/l.
# and sfpuc at -0.005 is within factor of 5 of stormwater
