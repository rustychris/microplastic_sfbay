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

plt.figure(1).clf()
fig,axs=plt.subplots(2,1,num=1)

axl=axs[0]
axw=axs[1]

bins=np.logspace(-2.5,1.7,50)
base=np.zeros(len(bins)-1)

all_l=[]
all_w=[]

for src_df in [storm_df,effluent_df]:
    lengths=src_df['Length.mm'].values
    lvalid=np.isfinite(lengths)
    all_l.append(lengths[lvalid])
    
    widths=src_df['Width.mm'].values
    wvalid=np.isfinite(widths)
    all_w.append(widths[wvalid])

lcounts,lbins,lcolls=axl.hist( all_l, bins=bins,stacked=True)
axl.set_xscale('log')
axl.set_xlabel('Length.mm')
axl.set_ylabel('Count')

wcounts,wbins,wcolls=axw.hist( all_w, bins=bins, stacked=True)
axw.set_xscale('log')
axw.set_xlabel('Width.mm')
axw.set_ylabel('Count')

axl.legend(lcolls,["Storm","Effluent"])

fig.savefig("storm_effluent_dimensions.png",dpi=150)

