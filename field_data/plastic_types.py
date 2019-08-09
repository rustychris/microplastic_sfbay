import six
import pandas as pd
import plastic_data

##
six.moves.reload_module(plastic_data)
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

types=combined.groupby(['pathway','PlasticType_Final'])['pathway'].count()

types=types.unstack(0,fill_value=0)

types['combined']=types['stormwater'] + types['effluent']
types=types.sort_values('combined',ascending=False)

table=types[:50].to_string()

fig=plt.figure(1)
fig.clf()
fig.text( 0.01, 0.98, table, va='top', family='monospace' )

fig.savefig("plastic-types.pdf")
fig.savefig("plastic-types.png")

##

with open('types-table.tex', 'wt') as fp:
    fp.write(types[:50].to_latex())

# plus manual edits in types-table-doc.tex
