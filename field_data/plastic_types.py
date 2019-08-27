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

# Select real and dupe samples, no blanks
valid_storm=(combined['pathway']=='stormwater')&(
    (combined['SampleType_AW']=='Field')
    | (combined['SampleType_AW']=='FieldDupe') )

valid_eff=(combined['pathway']=='effluent')&(
    combined['SampleTypeCode']=='Field')

valid=combined[ valid_storm|valid_eff ].copy()
##

types=valid.groupby(['pathway','PlasticType_Final'])['pathway'].count()

types=types.unstack(0,fill_value=0)

types['combined']=types['stormwater'] + types['effluent']
types=types.sort_values('combined',ascending=False)

table=types[:50].to_string()

## 
fig=plt.figure(1)
fig.clf()
fig.text( 0.01, 0.98, table, va='top', family='monospace' )

fig.savefig("plastic-types.pdf")
fig.savefig("plastic-types.png")

##

with open('types-table.tex', 'wt') as fp:
    fp.write(types[:50].to_latex())

# plus manual edits in types-table-doc.tex
##


# 2019-08-25 -- with the new data, stormwater and effluent combined.
# also removed blanks from these counts.


#                                        combined



# Not Identified                             7975

# No density for these 6:

# Anthropogenic (unknown base)                511

# Unknown                                     191

# Stearates, Lubricants, Waxes                143
# Anthropogenic (synthetic)                    91
# Inorganic natural material                   17
# Styrene copolymer                             4


# All of these have a density:

# Rubber                                     5863
# Polyethylene                                265
# Anthropogenic (cellulosic)                  102
# Cellulosic                                   42
# Organic natural material                     16
# Polyester                                    99
# Cotton                                       67
# Cellulose acetate                            62
# Polypropylene                                62
# Acrylic                                      49
# Nylon                                        22
# Ethylene/vinyl acetate copolymer             22
# Polyvinyl chloride                           22
# Glass                                        17
# Polyethylene terephthalate                   15
# Polyurethane                                 13
# Wool                                         11
# Polyvinyl butyral                             8
# Polystyrene                                   8
# Paint                                         8
# Methyl vinyl ether copolymers                 6 # 
# Acrylonitrile butadiene styrene               5 # (ABS)
# Polyethylene co-acrylic acid                  5
# Polyethylene/polypropylene copolymer          5
# Polytetrafluoroethylene                       3

# 25 so far in the table.

# more that have densities, but aren't called out in the
# table in the report.

# Polyvinyl acetate                             3
# Phenolic resin                                2

# Need to figure out what these are below:


# wide range is possible.  leave undefined.
# Silicone                                      2


# this is water soluble, so probably not complete information.
# Polyvinyl alcohol                             2


# Have not gotten around to these
# Polystyrene/acrylic copolymer                 2
# Poly(Aryletherketone)                         1
# Polyvinyl ether                               1
# Polyacrolein                                  1
# Polyethylenimine                              1
# Polycarbonate                                 1
# Fluoroelastomer                               1
# Unknown Potentially Rubber                    1
# Polyether block amide                         1
