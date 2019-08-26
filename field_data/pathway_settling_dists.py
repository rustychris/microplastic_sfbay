"""
Plot potential difference in distribution of w_s for
wastewater vs. stormwater.
"""
import pandas as pd
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np
import common

##
small=common.default_small

# 4 panel version of that
plt.figure(1).clf()
fig,axs=plt.subplots(2,1,num=1)
fig.set_size_inches([8.4,4.8],forward=True)

bins=10**np.linspace(np.log10(small),-0.5,50)
all_bins=np.r_[ -bins[::-1], bins ]
all_bins_mapped=common.map_bilog(all_bins)

for cat_i,(df,cat) in enumerate([ (storm_df,'Stormwater') ,
                                  (effluent_df,'Effluent')] ):
    w_s=df['w_s']
    w_s_val=w_s[np.isfinite(w_s)]
    N=len(w_s_val)

    w_s_mapped=common.map_bilog(w_s_val)
    
    # percentages
    color="0.5"
    
    ax=axs[cat_i]
    ax.hist(w_s_mapped,bins=all_bins_mapped,color=color)
    ax.text(0.03,0.95,cat,fontsize=12,
            va='top',transform=ax.transAxes)

    ax.axis(xmin=all_bins_mapped[0],xmax=all_bins_mapped[-1])
    
    plt.setp(ax.yaxis.get_ticklabels(),visible=0)
    ax.set_ylabel('Distribution')

    xticks=[common.map_bilog(x)
            for x in [-0.05,-0.005,
                      # 0.0,
                      0.005,0.05]]
    ax.set_xticks(xticks)
    labels=[common.mapped_label(y) for y in ax.get_xticks()]
    ax.set_xticklabels(labels)
    ax.set_yticks([])
    ax.set_xlabel('Rising/Sinking Velocity (m/s)')
    ax.axvline(0,color='k',lw=0.8,dashes=[8,16])

txtargs=dict(style='italic',fontweight='bold',fontsize=12)
axs[1].text(0.0,-0.06,"Floats",transform=ax.transAxes,va='top',ha='left',**txtargs)
axs[1].text(1,-0.06,"Sinks",transform=ax.transAxes,va='top',ha='right',**txtargs)
axs[1].text(0.5,-0.06,"Passive",transform=ax.transAxes,va='top',ha='center',**txtargs)
nticks=len(axs[1].get_xticklabels())

axs[0].xaxis.set_visible(0)

for ax in axs:
    ax.spines['top'].set_visible(0)
    ax.spines['right'].set_visible(0)

fig.subplots_adjust(left=0.06,right=0.94,top=0.96,wspace=0.10,hspace=0.08)

## 
fig.savefig("pathway_settling_dist2.png",dpi=150)



