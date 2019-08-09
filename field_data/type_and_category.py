import six
import pandas as pd
import plastic_data

##
from plastic_data import effluent_df,storm_df,combined


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

types=combined.groupby(['Category_Final','PlasticType_Final']).size()

# types=types.unstack(0,fill_value=0)

types=types.sort_values('combined',ascending=False)

table=types[:50].to_string()

