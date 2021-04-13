from matplotlib.ticker import FuncFormatter
import plastic_data
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
import common
##

manta_df=plastic_data.manta_df
sed_df=plastic_data.sediment_df

##

# First, how different are the w_s distributions?
plt.figure(1).clf()
fig,axs=plt.subplots(2,1,sharex=True,num=1)
fig.set_size_inches((7,6),forward=True)

bins=np.linspace(-3,3,30)

ret=axs[0].hist( common.map_bilog(manta_df.w_s[manta_df.w_s.notnull()]),bins=bins)
ret=axs[1].hist( common.map_bilog(sed_df.w_s[sed_df.w_s.notnull()]),bins=bins)

axs[1].xaxis.set_major_formatter(FuncFormatter(lambda s,i: common.mapped_label(s)))

axs[0].text(0.02,0.94,'Surface' ,transform=axs[0].transAxes,va='top',fontweight='bold',fontsize=14)
axs[1].text(0.02,0.94,'Sediment',transform=axs[1].transAxes,va='top',fontweight='bold',fontsize=14)

fig.subplots_adjust(left=0.12,right=0.98,top=0.98,bottom=0.15)
fig.set_size_inches([6,4],forward=True)
axs[1].set_xlabel('Settling / rise velocity')
axs[0].set_ylabel('Count')
axs[1].set_ylabel('Count')
fig.savefig('settling-manta_vs_sed.png',dpi=200)

##

# surface has some sinking particles, and sediment has some buoyant particles.
# I guess the immediate question is whether the buoyant particles in the
# sediment are significantly smaller?
# if larger surface-bound particles are fragmenting before sedimentation,
# expect to see only smaller particles in the sediment.
# But that is not distinguishable from the idea that the smaller, closer to
# neutral particles in the surface are the ones that can become sinkers
# and deposit.

# The modeling data are suggesting that particles entering the system include
# some particles near the surface but only for a couple days.


# What about size and category distributions of Bay vs. Coast?
# assuming coastal sources are negligible...
# within the Bay->Coast transit time, if biofilms remove particles from the surface
# the size distribution in the coast should be larger than in the Bay.
# if fragmentation is significant, we should see smaller particles in the ocean.
# Likewise, may seem some shift in category.

def valid(x):
    return x[x.notnull()]

def plot_length_comparison(manta_bay,manta_coast,ax,legend=True):
    bins=np.log10(np.logspace(-2,2,40))

    bay=np.log10(valid(manta_bay['Length.mm']))
    sns.distplot(bay,bins=bins,ax=ax,label='Bay (n=%d)'%len(bay))
    coast=np.log10(valid(manta_coast['Length.mm'])) 
    sns.distplot(coast,bins=bins,ax=ax,label='Coast (n=%d)'%len(coast))
    if legend:
        ax.legend(loc='upper left')

    sf=ScalarFormatter()
    sf._usetex=True
    ax.xaxis.set_ticks(np.arange(-2,3))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda s,i:'$'+sf.format_data(10.**s)+'$'))
    ax.set_xlabel('Length (mm)')

def fig_length_comparison(manta_bay,manta_coast,num):
    # First, overall dimension comparison.
    plt.figure(num).clf()
    fig,ax=plt.subplots(1,1,num=num)
    fig.set_size_inches([5,2.5],forward=True)
    plot_length_comparison(manta_bay,manta_coast,ax=ax)
    fig.subplots_adjust(bottom=0.23,top=0.90)
    return fig,ax

if 0:
    fig,ax=fig_length_comparison(
        manta_bay=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Bay') ],
        manta_coast=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Sanctuary') ],
        num=2)
    ax.set_title('All particles')

    for num,category in enumerate(manta_df.Category_Final.unique()):
        fig,ax=fig_length_comparison(
            manta_bay=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Bay') & (manta_df.Category_Final==category) ],
            manta_coast=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Sanctuary') & (manta_df.Category_Final==category)],
            num=3+num)
        ax.set_title('Category: %s'%category)

# This shows that
#  across all particles, the coast has larger particles.
#  fragments, fibers have the same length distribution (fibers are slightly smaller in the coast)
#  films and foams are larger in the coast.
#  spheres are slightly larger in the coast, but it's probably noise.
#  the fiber results may be misleading as a lot of that is probably contamination.
#  the film and fragment results are encouraging. films have a lot of area, so likely
#  to be fouled and sink, compared to fragments which have some bulk.
#  hard to say about foams.

# This suggests that loss in the Bay prefers small films and foams, and
# there is no loss of fragments.

# The size distributions are closer to log-normal than normal.

# Summary plot:
## 
cats=['Fiber',
      # 'Fiber Bundle' -- too few samples
      'Fragment', 'Foam', 'Film', 'Sphere']

plt.figure(20).clf()
fig,ax2d=plt.subplots(3,2,sharex=True,num=20)
axs=ax2d.ravel()

for num,category in enumerate(cats):
    plot_length_comparison(
        manta_bay=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Bay') & (manta_df.Category_Final==category) ],
        manta_coast=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Sanctuary') & (manta_df.Category_Final==category)],
        legend=False,
        ax=axs[num])
    axs[num].text(0.02,0.96,category,ha='left',va='top',transform=axs[num].transAxes)
    
plot_length_comparison(
    manta_bay=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Bay') ],
    manta_coast=manta_df[ manta_df.field_sample_p & (manta_df.Group1=='Sanctuary')],
    legend=False,
    ax=axs[-1])
axs[-1].text(0.02,0.96,'Total',ha='left',va='top',transform=axs[-1].transAxes)

for ax in ax2d[:-1,:].ravel():
    ax.set_xlabel('')
fig.subplots_adjust(top=0.95,left=0.07,right=0.97)

fig.savefig('manta-bay_vs_coast-length_distribution.png',dpi=200)
