import pandas as pd
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

# These are the sources that have PTM data, less the POTWs and Delta:
# 
storm_ptm_to_rwsm=[
    ("Alameda_Creek",['AlamedaCreek']),
    ("Guadalupe_Slo",['CalabazasCreek','SanTomas']), 
    ("San_Pablo_Cre",['SanPabloCreek']),
    ("Stevens_Creek",['StevensCreek']),
    ("Matadero_and_",['MataderoCreek','AdobeCreek'] ), 
    ("Petaluma_Rive",['PetalumaRiver']),
    ("Suisun_Slough",['SuisunSlough']),
    ("Coyote_Creek_",['LowerCoyoteCreekbelowAndersonDam']),
    # One of these is the Grizzly side, the other is labeled Sac side,
    # but is more like Denverton.
    # potential RWSM matches: GrizzlyIsland,
    #      # 
    ("Montezuma_Slo",['GrizzlyIsland']), # just a guess
    ("Montezuma_Slo_1ser",['SuisunIslands']), # just a guess
    ("Sonoma_Creek",['SonomaCreek']),
    ("Sulphur_Sprin", ['GoodyearSlough']), # ish
    ("San_Francisqu",['SanFrancisquitoCreek']),
    ("Napa_River",['NapaRiver']),
    # Near Mallard Slough
    # from https://www.cccleanwater.org/watersheds/watersheds-in-contra-costa-county,
    # maybe Willow or Mt Diablo? not a choice
    # http://cocowaterweb.org/wp-content/uploads/Watershed-Atlas.pdf
    # no additional names.
    ("unnamed07",['KirkerCreek']), # More like Willow Ck, but close enough?
    ("Glen_Echo_Cre",['GlenEchoCreek']),
    ("Old_Alameda_C",['WardandZeileCreeks']), # Ward is one of the main inputs to Old Alameda
    ("San_Leandro_C",['SanLeandroCreekBelowLakeChabot']),
    ("unnamed08",['GoodyearSlough']), # ish 
    ("Guadalupe_Riv",['GuadalupeRiver']),
    ("Pacheco_Creek",['WalnutCreek']), # ish...
    ("San_Lorenzo_C",['SanLorenzoCreek']),
    ("Steinberger_S", ['RedwoodCkandArroyoOjodeAguaCk','CordillerasCreek','PulgasCreek'])
]

# 21 stormwater sources

##

# For each of the PTM sources, aggregate the land use fractions
# from the RWSM regions, scale by landuse_coeffs, and write out to
# tsv

df_lump['area']=df['Tot. Area (km2)']
df_lump_by_name=df_lump.set_index( df['Watershed'])

## 
recs=[]
for ptm_src,rwsm_watersheds in storm_ptm_to_rwsm:
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
    
stormwater_concs=pd.DataFrame(recs)

stormwater_concs.to_csv('stormwater_concs-v00.csv',index=False)
