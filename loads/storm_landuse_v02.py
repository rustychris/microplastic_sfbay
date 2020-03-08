"""
v02: use hydrology_matchup to come with a better landuse distribution
"""

import pandas as pd
from stompy.spatial import wkb2shp


# Looking into land-use and how to better apportion stormwater loads.
# Source data: 12 watersheds
# Flows are already set, so just need to know
# (a) how do land-use types map to particle load
# (b) what is the land-use composition of each model inflow

# One issue will be how to deal with particle categories.
# -- does the existing analysis do this? Or is the land-use
#    factor just for total particle load?
# qualitatively, look at fig 2.3:
#   Line 12M had the greatest abundance, was more dominated by
#  fragment than any other site, and has a similar amount of
#  industrial to other sites.

# transcription of table 2.3 
landuse_coeffs=[('industrial',62),
                ('transportation',10),
                ('commercial',5),
                ('residential',1),
                ('ag_open',0.1), # includes open-space
]
# easier to index 
landuse_coeffs=pd.Series(dict(landuse_coeffs))

                
# The RWSM Pollutant Spreadsheet model includes
# a tab with a wide range of land-use categories,
# calculated to % WS area (plus other metrics)
# The closest matches to the above categories are:
# industrial: New industrial + Old industrial
# transportation: New transportation + Old transportation
# commercial: {New,Old} commercial {high,low}
# residential: {new,old} residential {high,med,low,rural}
# agriculture: {ag, open compacted, open uncompacted}

# poke around to understand any overlaps or missing
# Should use the % WS area, to match the fit from the
# report.
df=pd.read_excel('rwsm/pollutant/Pollutant Model/Pollutant Spreadsheet Model Calculations - Region.xlsx')

## 
LU_cols=[c for c in df.columns
         if '% WS Area' in c]

##

sum1=(df['LU Agriculture % WS Area']
      +df['LU Null % WS Area']
      +df['LU OpenCompacted % WS Area']
      +df['LU OpenUncompacted % WS Area']
      +df['LU SourceArea % WS Area'] # what is this?
      +df['LU Water % WS Area']
      +df['LU WaterRunoff % WS Area'] # and this??
      +df['LU newCommercialHigh % WS Area']
      +df['LU newCommercialLow % WS Area']
      +df['LU newIndustrial % WS Area']
      +df['LU newResidentialHigh % WS Area']
      +df['LU newResidentialLow % WS Area']
      +df['LU newResidentialMed % WS Area']
      +df['LU newResidentialRural % WS Area']
      +df['LU newTransportation % WS Area']
      )
sum2=( df['LU oldCommercialHigh % WS Area']
       + df['LU oldCommercialLow % WS Area']
       + df['LU oldIndustrial % WS Area']
       + df['LU oldResidentialHigh % WS Area']
       + df['LU oldResidentialLow % WS Area']
       + df['LU oldResidentialMed % WS Area']
       + df['LU oldResidentialRural % WS Area']
       + df['LU oldTransportation % WS Area']
)
sum_100=sum1+sum2

# This is ==100.0 for all watersheds.
df_lump=pd.DataFrame()
df_lump['industrial']=( df['LU oldIndustrial % WS Area']
                        + df['LU newIndustrial % WS Area'] )
df_lump['transportation']=( df['LU oldTransportation % WS Area']
                            + df['LU newTransportation % WS Area']
                            + df['LU SourceArea % WS Area'] # assuming this is airports
)


df_lump['commercial']=( df['LU oldCommercialHigh % WS Area']
                   + df['LU oldCommercialLow % WS Area']
                   + df['LU newCommercialHigh % WS Area']
                   + df['LU newCommercialLow % WS Area'] )
df_lump['residential']=(df['LU oldResidentialHigh % WS Area']
                   + df['LU oldResidentialLow % WS Area']
                   + df['LU oldResidentialMed % WS Area']
                   + df['LU oldResidentialRural % WS Area']
                   + df['LU newResidentialHigh % WS Area']
                   + df['LU newResidentialLow % WS Area']
                   + df['LU newResidentialMed % WS Area']
                   + df['LU newResidentialRural % WS Area'])
df_lump['ag_open']=( df['LU Agriculture % WS Area']
                +df['LU OpenCompacted % WS Area']
                +df['LU OpenUncompacted % WS Area'] )

total=df_lump.industrial+df_lump.transportation+df_lump.commercial+df_lump.residential+df_lump.ag_open

# These are relatively small, but unclear how to deal with them.
# omit, and scale up the remainder to sum to 1.0
pct_missing=( df['LU Null % WS Area']
              +df['LU Water % WS Area']
              +df['LU WaterRunoff % WS Area'] # and this??
              )

scale=1./total
for col in ['residential','ag_open','industrial','transportation','commercial']:
    df_lump[col] *= scale

# sanity check    
total_adj=df_lump.industrial+df_lump.transportation+df_lump.commercial+df_lump.residential+df_lump.ag_open

##

# So how hard is it to map the sources I have in the model to landuse from this?
# df_ord=df.sort_values('Tot. Area (km2)')
# print(df_ord.loc[:,['Watershed','Tot. Area (km2)']].values)

# That maps the watershed names.
# So why am I bothering to truncate these names, since I'm working
# with SUNTANS now?  Not sure, but I am.
# an additional translation step to go from the inflow names in hydro_match
# to the names that are actually specified to suntans via BC.nc

# These are just the PTM sources
sun_to_watershed=[
    ("Alameda_Creek","Alameda Creek"),
    ("Guadalupe_Slo","Guadalupe Slough"),
    ("San_Pablo_Cre",'San Pablo Creek'),
    ("Stevens_Creek",'Stevens Creek'),
    ("Matadero_and_",'Matadero and Adobe Creek'), 
    ("Petaluma_Rive",'Petaluma River'),
    ("Suisun_Slough",'Suisun Slough'),
    ("Coyote_Creek_",'Coyote Creek, Santa Clara'),
    # One of these is the Grizzly side, the other is labeled Sac side,
    # but is more like Denverton.
    # potential RWSM matches: GrizzlyIsland,
    ("Montezuma_Slo",'Montezuma Slough Grizzly'),
    ("Montezuma_Slo_1ser",'Montezuam Slough Sac'), 
    ("Sonoma_Creek",'Sonoma Creek'),
    ("Sulphur_Sprin", 'Sulphur Springs Creek'),
    ("San_Francisqu",'San Francisquito'),
    ("Napa_River",'Napa River'),
    # Near Mallard Slough
    # from https://www.cccleanwater.org/watersheds/watersheds-in-contra-costa-county,
    # maybe Willow or Mt Diablo? not a choice
    # http://cocowaterweb.org/wp-content/uploads/Watershed-Atlas.pdf
    # no additional names.
    ("unnamed07",'unnamed07'),
    ("Glen_Echo_Cre",'Glen Echo Creek'),
    ("Old_Alameda_C",'Old Alameda Creek'),
    ("San_Leandro_C",'San Leandro Creek'),
    ("unnamed08",'unnamed08'), 
    ("Guadalupe_Riv",'Guadalupe River'),
    ("Pacheco_Creek",'Pacheco Creek'), 
    ("San_Lorenzo_C",'San Lorenzo Creek'),
    ("Steinberger_S", 'Steinberger Slough')
]

ptm_to_watershed={a:b for a,b in sun_to_watershed}
watershed_to_ptm={b:a for a,b in sun_to_watershed}

##
# These are the sources that have PTM data, less the POTWs and Delta:

hydro_match=wkb2shp.shp2geom('../../sfb_ocean/suntans/grid-merge-suisun/hydrology_matchup.shp')
## 
storm_ptm_to_rwsm={}

for i in range(len(hydro_match)):
    inflow=hydro_match['inflow'][i]
    if inflow not in watershed_to_ptm:
        # Not a PTM source
        continue
    short=watershed_to_ptm[inflow]
    if short not in storm_ptm_to_rwsm:
        storm_ptm_to_rwsm[short]=[]
    name=hydro_match['newName2'][i]
    assert name.strip()
    storm_ptm_to_rwsm[short].append( name )

for k in storm_ptm_to_rwsm:
    assert k
    print(f"{k}")
    for r in storm_to_rwsm[k]:
        print(f"  {r}")

## 

# I have 72-ish inflows in the SF Bay model.
# 21+ of those have PTM sources.
# The RWSM has 219 watersheds.
# I now have a many:1 mapping of RWSM watershed to SUNTANS inflow
# Also watershed_inflow_locations has the watershed area for each SUNTANS inflow.

# Ultimately I need particle concentration, hourly, per PTM source, per settling class.
# 
# Option I:
#   For each SUNTANS inflow, determine amp_rwsm_to_sun factor based on the
#   ratio of the SUNTANS watershed area and the assigned RWSM watershed area(*).
#   Also determine the landuse composition.
#   To account for the mismatch between PTM source and SUNTANS sources, determine
#   a single amp_sun_to_ptm factor based on the total watershed area in SUNTANS
#   sources vs subset of those with PTM sources (*).
#   This approach is simplest, doesn't require scanning the BC files.

# Option II:
#   The RWSM data isn't for the same period, so net flows cannot be compared
#   between RWSM and suntans.  Calculate amp_rwsm_to_sun as above, using watershed
#   areas.
#   Sum the SUNTANS flows, per inflow, for the full simulation period.
#   Calculate amp_sun_to_ptm based on these flows, rather than watershed areas.
#   This seems like a minimal improvement with a significant additional annoyance.

# Option III:
#   Like I, but rather than a single amp_sun_to_ptm based on suntans areas, use
#   RWSM areas.


##

# For each of the PTM sources, aggregate the land use fractions
# from the RWSM regions, scale by landuse_coeffs, and write out to
# csv

df_lump['area']=df['Tot. Area (km2)']
df_lump_by_name=df_lump.set_index( df['Watershed'])

## 
recs=[]
for ptm_src,rwsm_watersheds in storm_ptm_to_rwsm.items():
    break
    rec=dict(source=ptm_src)
    
    rec['rwsm_regions']="|".join(rwsm_watersheds)

    rwsm_lumped = df_lump_by_name.loc[ rwsm_watersheds, :]

    areas=rwsm_lumped['area']
    lu_types=['industrial','transportation','commercial','residential','ag_open']

    fracs=rwsm_lumped.loc[:,lu_types].values * rwsm_lumped['area'].values[:,None]
    fracs=fracs.sum(axis=0)
    fracs=fracs/fracs.sum()
    net_coeff=(fracs*landuse_coeffs[lu_types].values).sum()
    rec['net_coeff']=net_coeff
    recs.append(rec)

## 
stormwater_concs=pd.DataFrame(recs)

stormwater_concs.to_csv('stormwater_concs-v01.csv',index=False)
