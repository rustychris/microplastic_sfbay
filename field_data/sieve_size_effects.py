"""
Look at particle size and sieve # across effluent, stormwater, and 
manta data in order to estimate size-selection bias.
"""
import pandas as pd
import plastic_data
import seaborn as sns
import matplotlib.pyplot as plt

## 
from plastic_data import effluent_df,storm_df,combined,manta_df, sediment_df


##

#  effluent_df.columns
#  Index(['StationCode', 'SampleTime', 'SampleID', 'Category_Final',
#         'PlasticType_Final', 'Color_Final', 'Treatment', 'SampleDate',
#         'SampleTypeCode', 'MatrixName', 'MethodName', 'AnalyteName',
#         'LabSampleID', 'Length.mm', 'Width.mm', 'Raman.ID', 'Remove_Flag',
#         'ScreenSize', 'Volume.L'],
#        dtype='object')

#  storm_df.columns
#  Index(['StationCode', 'Watershed', 'SampleID', 'SizeFraction_AW',
#         'SampleGross_AW', 'SampleType_AW', 'Number', 'LabSampleID', 'Length.mm',
#         'Width.mm', 'PlasticType_Final', 'Category_Final', 'Color_Final',
#         'Raman.ID', 'Remove_Flag'],
#        dtype='object')
#  

# common columns:
#  StationCode
#  Category_Final
#  PlasticType_Final
#  Color_Final
#  Length.mm
#  Width.mm
#  Raman.ID
#

##

# combined['pathway'] ~ ['effluent', 'stormwater', 'sediment', 'manta', 'fish']
# careful about non-field samples.  check field_sample_p
# Seems that only effluent has ScreenSize.
# And 1700 particles have ScreenSize=125, and 1724 particles have ScreenSize=355
# 86 particles with ScreenSize=0.
# On the surface that suggests 50.4% efficiency for manta.

# Effluent, Stormwater, and Sediment all have SizeFraction column
# For Effluent, same as ScreenSize.
# For stormwater, sediment, it's all over the place.

# stormwater: 
#   SizeFraction
#   125                582
#   125-355um          862
#   125um-355um       6678

#   355                188
#   355-500um          258
#   355um-500um       1285
#   500-1mm            276
#   500um-1mm         1229
#   >1                  22
#   >1m                 22
#   >1mm               761

#   125um-->1mm          7
#   combined (N/A)      59

# Combine these to estimate total count < 355, vs. > 355.
#  < 355um: 582+862+6678 = 8122
#  > 355um: 188+258+1285+276+1229+22+22+761 = 4041
# So the stormwater data suggests manta would be 4041/(4041+8122) = 0.332 efficient.

# For sediment data:
# SizeFraction
# 125 µm - 355 µm      68
# 125 μm-355 μm       320
# 125-355            3440
# ==> 68 + 320 + 3440 = 3828 below manta size.

# 1.00mm               11
# 1.0mm                 1
# 355 µm - 500 µm      31
# 355 μm-500 μm       115
# 355-500             873
# 355um               108
# 356 μm-500 μm         1
# 500 µm - 1 mm        19
# 500 μm-1 mm          77
# 500-1               167
# 500um               123
# 500um                 6
# >1 mm                44
# >1mm                 96
# >500                469
# >500                  1

# ==> sums to 2142
# 2142 / (2141 + 3828) ==> 0.359

# 125-500              52 # n/a

##

size_frac_to_manta_p={125:False,355:True,
                      0:None,'125um-355um':False,
                      '355um-500um':True,'500um-1mm':True,'>1mm':True,
                      '125':False,np.nan:None,'>1m':True,
                      'combined (N/A)':None,'500-1mm':True,'125-355um':False,
                      '355-500um':True, '>1':True, '355':True,
                      '125um-->1mm':None, '>500':True, '355-500':True,
                      '125-355':False, '500-1':True, '500 μm-1 mm':True,
                      '355 μm-500 μm':True,
                      '125 μm-355 μm':False, '>500 ':True, '1.00mm':True,
                      '500um':True, '355um':True, '>1 mm':True,
                      '125-500':None,'500 µm - 1 mm':True, '355 µm - 500 µm':True,
                      '125 µm - 355 µm':False,'356 μm-500 μm':True, '1.0mm':True, '500um ':True,
                      '1mm':True, '125um':False, '1 mm':True, '500':True,
                      '212':False, '126':False}

combined['manta_p']=combined['SizeFraction'].map(size_frac_to_manta_p)

##

# Oddly, this shows that the manta samples, even though collected in 355um
# mesh, most samples from manta were sieve at smaller size class.


print(combined.groupby(['pathway','manta_p']).size())

nonfiber=combined[ ~combined['Category_Final'].isin(['Fiber','Fiber Bundle']) ]
print(nonfiber.groupby(['pathway','manta_p']).size())


floaters=combined[ combined['w_s']<0 ]
floaters.groupby(['pathway','manta_p']).size()

# 

