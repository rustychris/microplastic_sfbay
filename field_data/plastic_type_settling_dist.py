"""
Plot where different plastics tend to fall on the settling velocity spectrum
"""
import pandas as pd
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np
import common

##
from matplotlib import gridspec
gs=gridspec.GridSpec(5,1)

# 4 panel version of that
#plt.close(2)
fig=plt.figure(2)
fig.set_size_inches([8.4,6.0],forward=True)
fig.clf()

small=5e-4

bins=10**np.linspace(np.log10(small),-0.5,50)
# 2019-09-06: Don pointed out a discrepancy in the plots, and it traces back to here.
# used to use all samples.
#  w_s=combined['w_s']
#  w_s_val=w_s[np.isfinite(w_s)]
# It had already been using the right subset for the bars, but the histogram used
# the wrong set.
sel=combined['pathway'].isin(['stormwater','effluent']) & combined['field_sample_p']
valid=combined[sel].copy()
w_s=valid['w_s']
w_s_val=w_s[np.isfinite(w_s)]

# HERE -
# how hard is it to weights these according to loads?

w_s_mapped=common.map_bilog(w_s_val)
all_bins=np.r_[ -bins[::-1], bins ]
all_bins_mapped=common.map_bilog(all_bins)

ax=fig.add_subplot(gs[-1,:])
ax.cla()

ax.hist(w_s_mapped,bins=all_bins_mapped,color='0.5')

ax.axis(xmin=all_bins_mapped[0],xmax=all_bins_mapped[-1])

xticks=[common.map_bilog(x)
        for x in [-0.05,-0.005,0.0,0.005,0.05]]
ax.set_xticks(xticks)

labels=[common.mapped_label(y) for y in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_xlabel('Rising/Sinking Velocity (m/s)')
plt.setp(ax.get_yticklabels(),visible=0)
ax.axvline(0,ls='--',color='k',lw=0.5,zorder=-2)

common.set_bold_labels(ax,y=-0.11)

ax.set_ylabel('Abundance')

ax_typ=fig.add_subplot(gs[:-1,:],sharex=ax)
plt.setp(ax_typ.get_xticklabels(),visible=0)

# choose the 20 most common types
# split foams out separately

my_types=valid['PlasticType_Final'].copy()
foams=valid['Category_Final']=='Foam'

# condense several of the "unknown" to one
all_unknown=my_types.isin(['Not Identified','Unknown','Anthropogenic (unknown base)','Anthropogenic (synthetic)'])
my_types[all_unknown] ='Unknown'
my_types[foams]=my_types[foams] +" (foam)"

valid['my_type']=my_types
# For the purposes of showing w_s distribution, group the types slightly differently.
types=valid.groupby('my_type').size()
types=types.sort_values(ascending=False)

count=0
for idx,typ in enumerate(types.index):
    if count>=26: break
    print(typ)
    subgroup=valid.iloc[ (valid.my_type==typ).values, :]
    sub_w_s=subgroup['w_s'].values
    sub_w_s=sub_w_s[np.isfinite(sub_w_s)]
    if not len(sub_w_s):
        print(f"Type {typ} has no settling velocities.  Skip")
        continue
    sub_w_s=common.map_bilog(sub_w_s)
    sub_25,sub_75=np.percentile(sub_w_s,[25,75])

    if sub_25==sub_75:
        sub_25-=0.01
        sub_75+=0.01
    # plot some dots or something to show distribution of w_s.
    y=count
    ax_typ.plot([sub_25,sub_75],[y,y],lw=3,
                solid_capstyle='round')
    ax_typ.text(sub_75,y," "+typ,ha='left',va='center',fontsize=9)
    plt.setp(ax_typ.spines.values(),visible=0)
    ax_typ.xaxis.set_visible(0)
    ax_typ.yaxis.set_visible(0)
    count+=1

fig.subplots_adjust(top=0.99,left=0.1,hspace=0.01)
## 
fig.savefig('plastic_type_settling_dist.png',dpi=150)
