"""
Plot where different categories fall on the settling velocity spectrum
"""
import pandas as pd
import six
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np


##
from matplotlib import gridspec
gs=gridspec.GridSpec(3,1)

plt.close(1)
fig=plt.figure(1)
fig.set_size_inches([8.4,4.8],forward=True)
fig.clf()

import common
six.moves.reload_module(common)
from common import map_bilog, unmap_bilog, mapped_label, default_small
small=default_small

bins=10**np.linspace(np.log10(small),-0.5,50)
w_s=combined['w_s']
w_s_val=w_s[np.isfinite(w_s)]

    
w_s_mapped=map_bilog(w_s_val)
all_bins=np.r_[ -bins[::-1], bins ]
all_bins_mapped=map_bilog(all_bins)

ax=fig.add_subplot(gs[-1,:])
ax.cla()

ax.hist(w_s_mapped,bins=all_bins_mapped,color='0.5')

ax.axis(xmin=all_bins_mapped[0],xmax=all_bins_mapped[-1])

xticks=[map_bilog(x)
        for x in [-0.05,-0.005,0.0,0.005,0.05]]
ax.set_xticks(xticks)

labels=[mapped_label(y) for y in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_xlabel('Velocity (m/s)')
plt.setp(ax.get_yticklabels(),visible=0)
ax.axvline(0,ls='--',color='k',lw=0.5,zorder=-2)

ax.text(0.0,-0.1,"$Floats$",transform=ax.transAxes,va='top',ha='left')
ax.text(1,-0.1,"$Sinks$",transform=ax.transAxes,va='top',ha='right')

ax_typ=fig.add_subplot(gs[:-1,:],sharex=ax)
plt.setp(ax_typ.get_xticklabels(),visible=0)

# not actually types, but categories
types=combined.groupby('Category_Final')['pathway'].count()
types=types.sort_values(ascending=False)

labels=[]
w_s=[]
for idx,typ in enumerate(types.index):
    print(typ)
    subgroup=combined.iloc[ (combined.Category_Final==typ).values, :]
    sub_w_s=subgroup['w_s'].values
    sub_w_s=sub_w_s[np.isfinite(sub_w_s)]
    sub_w_s=map_bilog(sub_w_s)
    w_s.append(sub_w_s)
    labels.append(typ)
boxes=ax_typ.boxplot(w_s,vert=False,labels=labels,
                     patch_artist=True,
                     medianprops=dict(visible=0),
                     whis=[5,95],
                     boxprops=dict(lw=0),
                     flierprops=dict(marker='.',markersize=3))

# I want nicer colors
cyc=plt.rcParams['axes.prop_cycle']
for box,sty in zip(boxes['boxes'],cyc):
    box.set_facecolor(sty['color'])

    
plt.setp(ax_typ.spines.values(),visible=0)
ax_typ.xaxis.set_visible(0)

fig.savefig('category_settling_dist.png',dpi=150)

##

# spheres are dominated by polyethylene, but significant number of glass
