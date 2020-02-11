"""
v02: try to refactor, process stormwater and effluent at the same time.
process each category independently, sum at the end.
"""
import os
import re

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

##

#version='v03-nofiber'
version='v04'

def trim_category(df,cat_version):
    if cat_version is None:
        return df
    if 'nofiber' in cat_version:
        slim=df[ df['Category_Final']!='Fiber' ]
        print(f"Removing fibers: {len(df):5} => {len(slim):5} particles")
        return slim.copy()
    for cat,official in [('fiber','Fiber'),
                         ('fiber_bundle','Fiber Bundle'),
                         ('film','Film'),
                         ('foam','Foam'),
                         ('fragment','Fragment'),
                         ('sphere','Sphere') ]:
        if cat == cat_version:
            slim=df[ df['Category_Final']==official ]
            print(f"Choosing only {cat}: {len(df):5} => {len(slim):5} particles")
            return slim.copy()
    assert False

# 5 particles with nan category??
categories=plastic_data.storm_df['Category_Final'].unique()
cat_df=pd.DataFrame({'category':categories})

def calc_blank_count_per_category(blanks,n_samples=1,column='Category_Final'):
    df=pd.DataFrame({'category':categories}).set_index('category')
    df['count']=0

    for key,grp in blanks.groupby(column).groups.items():
        df.loc[key,'count']=len(grp)
    df['rate']=df['count']/n_samples
    return df

# Specify what w_s values will be used for the groups in the PTM
w_s_centers=np.array([0.05,0.005,0.0005,
                      0,
                      -0.05,-0.005,-0.0005])
w_s_centers.sort()

# These numbers from tabel 1 in Stormwater Chapter 05-10-2019.docx
# StationCode seems to be just a more descriptive name for what is in
# Watershed
sample_volumes_l={
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

# The sample volumes, from Wastewater Chapter 05-09-2019.docx
# each entry maps a regular expression to the sample volume.
# all sample_ids matching the regex should be combined.
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

    
for cat_version in [
        None, # standard
        # "nofiber",
        "fiber",
        "fiber_bundle",
        "film",
        "foam",
        "fragment",
        "sphere"
]:
    print(f"----------------------{cat_version}----------------")
    ds=xr.Dataset()
    ds['w_s']=('w_s',),w_s_centers
    ds.w_s.attrs['name']="Center of bin, settling velocity, positive down, m/s"

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

    df_all=plastic_data.storm_df
    storm_blank_types=['FIELDQA','LABQA']

    blanks_all=df_all[ (df_all['SampleID']!='Bottle Blank') & df_all['SampleType_AW'].isin(storm_blank_types) ].copy() 
    field_all=df_all[df_all['field_sample_p']].copy()
    
    n_blank_samples = len(blanks_all['SampleID'].unique())
    n_field_samples = len(field_all['SampleID'].unique())

    # trim the blank/field dataframes instead 
    # df=trim_category(plastic_data.storm_df,cat_version)
    
    # Goal is to come up with a table of
    # particles/L binned by settling velocity.
    # Limit to field samples

    # 1317 l
    total_sample_volume_l=np.sum( list(sample_volumes_l.values()))
    print(f"Total sample volume, stormwater: {total_sample_volume_l} l")

    # sel=(df['field_sample_p'] & np.isfinite(df['w_s'])).values
    # had been limiting this to particles with w_s, but that was too early
    field=trim_category(field_all,cat_version) 
    blanks=trim_category(blanks_all,cat_version)

    # Had been counting n_blank, n_field here, but that may lose some sample ids
    # esp. when it's a rare category
    #n_blank_samples=len(blanks['SampleID'].unique())
    #n_field_samples=len(field['SampleID'].unique())

    print("--- Stormwater Blanks ---")
    if 'LABQA' in storm_blank_types:
        print(f"  Considering Field & Lab Blanks, {n_blank_samples:3} distinct blank samples")
    else:
        # 2
        print(f"  Considering Field Blanks only, {n_blank_samples:3} distinct blank samples")

    # 15
    print(f"                                  {n_field_samples:3} distinct field samples")

    # Assume that each sample accumulates the same counts of erroneous particles.
    # In the current setup, this could be done outside the loop.  But in case the different
    # versions start to include different blank handling, this stays inside the loop.
    blank_rates=calc_blank_count_per_category(blanks, n_blank_samples)

    # All particles start with weight of 1.0
    field['derated']=1.0

    for cat,blank_rate in blank_rates['rate'].iteritems():
        print(f"Stormwater, category={cat}.  Per-blank count {blank_rate}")
        cat_sel=field['Category_Final']==cat    
        n_field_per_cat=cat_sel.sum()
        print(f"    field count over {n_field_samples:3} samples: {n_field_per_cat}, {n_field_per_cat/n_field_samples:.2f} per sample")
        if n_field_per_cat>0:
            real_derate=(n_field_per_cat - n_field_samples*blank_rate) / n_field_per_cat
            derate=max(0.0,real_derate)
            print(f"    de-rating: {derate:.3f} ({real_derate:.3f})")
        else:
            print(f"    no field occurrences, de-rating irrelevant")
            derate=1.0
        print()

        field.loc[ cat_sel, 'derated']=derate

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

    for w_s, conc, conc_raw in zip(w_s_centers,bin_conc,bin_conc_raw):
        print(f" w_s ~ {w_s: .5f}m/s    conc ~ {conc:.3f} particles/l, raw ~ {conc_raw:.3f}")
    print(f"  Total                      {bin_conc.sum():.3f} particles/l, raw ~ {bin_conc_raw.sum():.3f}")


    # fill in the dataset
    ds['conc'].sel(source='stormwater').values[:]= bin_conc
    ds['conc_raw'].sel(source='stormwater').values[:]= bin_conc_raw

    # Output: (I think this is for std)
    #
    #  w_s ~ -0.05000m/s    conc ~ 0.245 particles/l, raw ~ 0.245
    #  w_s ~ -0.00500m/s    conc ~ 0.514 particles/l, raw ~ 0.517
    #  w_s ~ -0.00050m/s    conc ~ 0.543 particles/l, raw ~ 0.569
    #  w_s ~  0.00000m/s    conc ~ 0.016 particles/l, raw ~ 0.018
    #  w_s ~  0.00050m/s    conc ~ 1.097 particles/l, raw ~ 1.217
    #  w_s ~  0.00500m/s    conc ~ 6.195 particles/l, raw ~ 6.252
    #  w_s ~  0.05000m/s    conc ~ 0.568 particles/l, raw ~ 0.569
    #   Total                      9.179 particles/l, raw ~ 9.386

    # And wastewater:

    # StationCode: 'CCCSD', 'EBDA', 'EBMUD', 'FSSD', 'LABQA', 'PA', 'SFPUC', 'SJ', 'SUNN'
    # SampleID: includes date, station and sieve size.
    # SampleTypeCode: 'Field', 'LabBlank', 'FieldBlank'
    #    - select on this

    # Make sure that all of the sample_IDs line up with exactly one pattern, but before
    # trimming to subset of particles
    for sample_id in plastic_data.effluent_df.SampleID.unique():
        matches= [ patt
                   for patt in sample_volume_patts.keys()
                   if re.match(patt,sample_id) is not None ]
        if len(matches)!=1:
            print(f"Got {len(matches)} matches for sample id {sample_id}")
            print(matches)
            raise Exception("Matches are off")
    
    df_all=plastic_data.effluent_df
   
    # Summarize the blank counts by category
    include_lab_blanks=True

    if include_lab_blanks:
        blanks_all=df_all[ ~df_all['field_sample_p'] ].copy() # use field and lab blanks
    else:
        blanks_all=df_all[ (df_all['StationCode']!='LABQA') & (~df_all['field_sample_p']) ].copy() # only use field blanks

    field_all=df_all[df_all['field_sample_p']].copy()

    n_blank_samples=len(blanks_all['SampleID'].unique())
    n_field_samples=len(field_all['SampleID'].unique())

    print("--- Effluent Blanks ---")
    if include_lab_blanks:
        print(f"  Considering Field & Lab Blanks, {n_blank_samples:3} distinct blank samples")
    else:
        print(f"  Considering Field Blanks only, {n_blank_samples:3} distinct blank samples")

    # 15
    print(f"                                 {n_field_samples:3} distinct field samples")

    print(" Field, effluent: ",end="")
    field=trim_category(field_all,cat_version)
    print("Blanks, effluent: ",end="")
    blanks=trim_category(blanks_all,cat_version)
    
    # Assume that each sample accumulates the same counts of erroneous particles.
    blank_rates=calc_blank_count_per_category(blanks, n_blank_samples)

    print("Per-blank sample, per-category counts")
    print(blank_rates)
    print()

    # Derate
    # All particles start with weight of 1.0
    # The derating has to be done per-station.
    field['derated']=1.0

    stationcodes=field_all.StationCode.unique()

    for station in stationcodes:
        # this might be empty
        field_station=field[field.StationCode==station].copy()
        total_count=len(field_station)
        field_valid=np.isfinite(field_station.w_s)
        bins=utils.nearest(w_s_centers,field_station[field_valid].w_s.values)
        w_s_count=field_valid.sum()

        # what was the total volume sampled for this station?
        # here we have to go back to field_all, since some sample IDs
        # for this station may have zero counts
        sample_volume_l=0
        field_all_station=field_all[ field_all.StationCode==station ]
        sample_ids=field_all_station.SampleID.unique()
        for patt in sample_volume_patts.keys(): # does this volume go with this station?
            for samp_id in sample_ids: # yes, iff 1 or more sample_ids match
                 if re.match(patt,samp_id) is not None:
                     sample_volume_l+=sample_volume_patts[patt]
                     break
        print(f"{station:20} total sample volume {sample_volume_l:7} l")

        # Derate
        n_field_samples=len(sample_ids)

        for cat,blank_rate in blank_rates['rate'].iteritems():
            print(f"Effluent, category={cat}.  Per-blank count {blank_rate:.3f}")
            cat_sel=field_station['Category_Final']==cat    
            n_field_per_cat=cat_sel.sum() # total number of field particles at this station of this category
            print(f"    At {station}, field count over {n_field_samples:3} samples: {n_field_per_cat}, {n_field_per_cat/n_field_samples:.2f} per sample")
            if n_field_per_cat>0:
                # i.e. 8 particles over 2 samples, with a blank rate of 0.5/sample
                #            (8               - 2*0.5) / 8
                # so I want my 8 particles to get counted as 8-2*0.5 = 7 particles because
                # I'm estimating 1 particle is contamination. And thus I weigh each of those
                # original 8 particles by 7/8 so they sum to 7 instead of 8.
                real_derate=(n_field_per_cat - n_field_samples*blank_rate) / n_field_per_cat
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

    out_fn=f'plastic_loads-7classes-{version}{cat_version}.nc'
    os.path.exists(out_fn) and os.unlink(out_fn)
    ds.to_netcdf(out_fn)

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

##

# Sanity check of new per-category adjusted loads vs. previous
# std and no fiber, and current std, and no fiber.
ds_v04_nof=xr.open_dataset('plastic_loads-7classes-v04nofiber.nc')
ds_v04_std=xr.open_dataset('plastic_loads-7classes-v04None.nc')

ds_v03_nof=xr.open_dataset('plastic_loads-7classes-v03-nofiber.nc')
ds_v03_std=xr.open_dataset('plastic_loads-7classes-v03.nc')

## 
ds_sum=None
for fn in [
        'plastic_loads-7classes-v04fiber_bundle.nc',
        'plastic_loads-7classes-v04fiber.nc',
        'plastic_loads-7classes-v04film.nc',
        'plastic_loads-7classes-v04foam.nc',
        'plastic_loads-7classes-v04fragment.nc',
        'plastic_loads-7classes-v04sphere.nc']:
    ds_v04_single=xr.open_dataset(fn)
    if ds_sum is None:
        ds_sum=ds_v04_single.copy()
    else:
        assert np.all( ds_sum.w_s.values==ds_v04_single.w_s.values)
        assert np.all( ds_sum.source.values==ds_v04_single.source.values)
        ds_sum['conc'].values += ds_v04_single.conc.values
        ds_sum['conc_raw'].values += ds_v04_single.conc_raw.values
    ds_v04_single.close()
        
##
fld='conc_raw'
plt.figure(1).clf()
fig,axs=plt.subplots(2,3,num=1)
for ax,conc in [ (axs[0,0], ds_v04_std[fld].values ),
                 (axs[0,1], ds_sum[fld].values ),
                 (axs[0,2], ds_v03_std[fld].values ) ]:
    ax.imshow(np.log10(conc.clip(1e-4)),clim=[-4,0])

for ax,conc in [ (axs[1,0], ds_v04_std[fld].values - ds_sum[fld].values ),
                 (axs[1,1], ds_sum[fld].values - ds_v03_std[fld].values ),
                 (axs[1,2], ds_v03_std[fld].values - ds_v04_std[fld].values ) ]:
    ax.imshow( conc,cmap='seismic',clim=[-1e-4,1e-4])

# okay - so v03 and v04 match up exact.    
# but the sum of the individuals do not match up.
# and the error falls to both sides of the line.
# error exists in conc_raw, too.  so not a derating issue.
# the issue crops up for stormwater and wastewater.
# so debug just for stormwater

# There is a difference here, but I think it is in some sense real.
# The w_s distributions are scaled up based on the fraction of particles
# with an estimated w_s.
# this is after multiple categories are aggregated.
# The question then is whether it is more appropriate to scale up the
# bin distributions per category or in sum.
# it's entirely possible that particles of one category are more likely to
# get a valid w_s than particles in another category.
# if that's the case, then it's better to take the per-category numbers
# and scale up the w_s distributions per-category.

plt.figure(2).clf()
plt.loglog( ds_v04_std[fld].values.ravel(),
           ds_sum[fld].values.ravel(), 'g.')
plt.loglog( [1e-5,1],[1e-5,1],'k-',lw=0.5)


##

#   ----------------------None----------------
#   Total sample volume, stormwater: 1317 l GOOD
#   --- Stormwater Blanks ---
#     Considering Field & Lab Blanks,   4 distinct blank samples GOOD 
#                                      15 distinct field samples GOOD
#   Stormwater, category=Fiber.  Per-blank count 36.75 GOOD
#       field count over  15 samples: 4794, 319.60 per sample GOOD
#       de-rating: 0.885 (0.885) GOOD
#   
#   Stormwater, category=Fragment.  Per-blank count 1.25 GOOD
#       field count over  15 samples: 7309, 487.27 per sample GOOD
#       de-rating: 0.997 (0.997) GOOD
#   
#   Stormwater, category=Foam.  Per-blank count 0.0 GOOD
#       field count over  15 samples: 46, 3.07 per sample GOOD
#       de-rating: 1.000 (1.000) GOOD
#   
#   Stormwater, category=Film.  Per-blank count 0.5 GOOD
#       field count over  15 samples: 113, 7.53 per sample GOOD
#       de-rating: 0.934 (0.934) GOOD
#   
#   Stormwater, category=Sphere.  Per-blank count 0.0 GOOD
#       field count over  15 samples: 25, 1.67 per sample GOOD
#       de-rating: 1.000 (1.000) GOOD
#   
#   Stormwater, category=nan.  Per-blank count 0.0 GOOD
#       field count over  15 samples: 0, 0.00 per sample GOOD
#       no field occurrences, de-rating irrelevant GOOD
#   
#   Stormwater, category=Fiber Bundle.  Per-blank count 0.25 GOOD
#       field count over  15 samples: 71, 4.73 per sample GOOD
#       de-rating: 0.947 (0.947) GOOD
#   
#    w_s ~ -0.05000m/s    conc ~ 0.245 particles/l, raw ~ 0.245
#    w_s ~ -0.00500m/s    conc ~ 0.514 particles/l, raw ~ 0.517
#    w_s ~ -0.00050m/s    conc ~ 0.543 particles/l, raw ~ 0.569
#    w_s ~  0.00000m/s    conc ~ 0.016 particles/l, raw ~ 0.018
#    w_s ~  0.00050m/s    conc ~ 1.097 particles/l, raw ~ 1.217
#    w_s ~  0.00500m/s    conc ~ 6.195 particles/l, raw ~ 6.252
#    w_s ~  0.05000m/s    conc ~ 0.568 particles/l, raw ~ 0.569
#     Total                      9.179 particles/l, raw ~ 9.386

# per category:
fiber:          3.222   3.640
fiber bundle:   0.051   0.054
film:           0.080   0.086
foam:           0.035   0.035
fragment:       5.535   5.550
sphere:         0.019   0.019
# sums to 8.942 adusted
# 9.384 raw.  that's actually almost exact. scales out to error of 2.64
# particles across the full sample volume.
# there are 5 nan category particles in total.

#   ----------------------fiber----------------
#   Total sample volume, stormwater: 1317 l GOOD
#   Choosing only fiber: 12362 =>  4794 particles 
#   Choosing only fiber:   156 =>   147 particles 
#   --- Stormwater Blanks ---
#     Considering Field & Lab Blanks,   4 distinct blank samples GOOD
#                                      15 distinct field samples GOOD
#   Stormwater, category=Fiber.  Per-blank count 36.75
#       field count over  15 samples: 4794, 319.60 per sample GOOD
#       de-rating: 0.885 (0.885) GOOD
#   
#   Stormwater, category=Fragment.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#   Stormwater, category=Foam.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#   Stormwater, category=Film.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#   Stormwater, category=Sphere.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#   Stormwater, category=nan.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#   Stormwater, category=Fiber Bundle.  Per-blank count 0.0
#       field count over  15 samples: 0, 0.00 per sample
#       no field occurrences, de-rating irrelevant
#   
#    w_s ~ -0.05000m/s    conc ~ 0.000 particles/l, raw ~ 0.000
#    w_s ~ -0.00500m/s    conc ~ 0.000 particles/l, raw ~ 0.000
#    w_s ~ -0.00050m/s    conc ~ 0.227 particles/l, raw ~ 0.257
#    w_s ~  0.00000m/s    conc ~ 0.021 particles/l, raw ~ 0.023
#    w_s ~  0.00050m/s    conc ~ 2.416 particles/l, raw ~ 2.730
#    w_s ~  0.00500m/s    conc ~ 0.558 particles/l, raw ~ 0.630
#    w_s ~  0.05000m/s    conc ~ 0.000 particles/l, raw ~ 0.000
#     Total                      3.222 particles/l, raw ~ 3.640
#   
