"""
Plot where different plastics tend to fall on the settling velocity spectrum
"""
import pandas as pd
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np


##
from matplotlib import gridspec
gs=gridspec.GridSpec(3,1)

# 4 panel version of that
plt.close(1)
fig=plt.figure(1)
fig.set_size_inches([8.4,4.8],forward=True)
fig.clf()

small=5e-4

bins=10**np.linspace(np.log10(small),-0.5,50)
w_s=combined['w_s']
w_s_val=w_s[np.isfinite(w_s)]


def map_bilog(x):
    return np.sign(x)*np.log10(np.abs(x).clip(small,np.inf)/small)
def unmap_bilog(y):
    return np.sign(y)*small*10**np.abs(y)

from matplotlib.ticker import ScalarFormatter
sf=ScalarFormatter()
sf._usetex=True
def mapped_label(y):
    if y==0.0:
        s="\pm " + sf.format_data(small)
    else:
        s=sf.format_data(unmap_bilog(y))
    return "$"+s+"$"

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

ax.text(0.02,0.95,"Rises",transform=ax.transAxes,va='top')
ax.text(0.92,0.95,"Sinks",transform=ax.transAxes,va='top',ha='right')


ax_typ=fig.add_subplot(gs[:-1,:],sharex=ax)
plt.setp(ax_typ.get_xticklabels(),visible=0)

# choose the 20 most common types
types=combined.groupby('PlasticType_Final')['pathway'].count()
types=types.sort_values(ascending=False)

for idx,typ in enumerate(types.index[:20]):
    print(typ)
    subgroup=combined.iloc[ (combined.PlasticType_Final==typ).values, :]
    sub_w_s=subgroup['w_s'].values
    sub_w_s=sub_w_s[np.isfinite(sub_w_s)]
    sub_w_s=map_bilog(sub_w_s)
    sub_25,sub_75=np.percentile(sub_w_s,[25,75])

    if typ=='Unknown':
        # w_s is just from the Foam unknowns.
        typ="Unknown foam"

    if sub_25==sub_75:
        sub_25-=0.01
        sub_75+=0.01
    # plot some dots or something to show distribution of w_s.
    ax_typ.plot([sub_25,sub_75],[idx,idx],lw=3,
                solid_capstyle='round')
    ax_typ.text(sub_75,idx," "+typ,ha='left',va='center')
    plt.setp(ax_typ.spines.values(),visible=0)
    ax_typ.xaxis.set_visible(0)
    ax_typ.yaxis.set_visible(0)

fig.subplots_adjust(top=0.95,left=0.1)    
fig.savefig('plastic_type_settling_dist.png',dpi=150)
