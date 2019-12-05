"""
Estimate undercounts in manta data due to difference in size.
"""

import plastic_data
from plastic_data import manta_df, storm_df, effluent_df, combined
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib import ticker

##

def cap(s):
    return s[0].upper() + s[1:]

# For starters, just look at size distributions
plt.figure(1).clf()
fig,axs=plt.subplots(4,1,num=1,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    values=combined.loc[sel,'Length.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-1.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Length (mm)')
axs[0].set_title('All Particles')

# That shows effluent and stormwater having similar distributions,
# and manta is shifted right somewhat.


##

plt.figure(2).clf()
fig,axs=plt.subplots(4,1,num=2,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    values=combined.loc[sel,'Width.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-2.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Width (mm)')
axs[0].set_title('All Particles')

# That shows effluent and stormwater having similar distributions,
# and manta is shifted right somewhat.

##

plt.figure(3).clf()
fig,axs=plt.subplots(4,1,num=3,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    sel=sel & np.isfinite(combined.w_s) & (combined.w_s<0)
    values=combined.loc[sel,'Length.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-2.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Length (mm)')
axs[0].set_title('Buoyant Particles')

##

# That shows effluent and stormwater having similar distributions,
# and manta is shifted right somewhat.
plt.figure(4).clf()
fig,axs=plt.subplots(4,1,num=4,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    sel=sel & np.isfinite(combined.w_s) & (combined.w_s<0)
    values=combined.loc[sel,'Width.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-2.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Width (mm)')

axs[0].set_title('Buoyant Particles')

## 
plt.figure(5).clf()
fig,axs=plt.subplots(4,1,num=5,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    # sel=sel & np.isfinite(combined.w_s) & (combined.w_s<0)
    sel=sel & (combined.Category_Final!='Fiber')
    values=combined.loc[sel,'Length.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-2.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Length (mm)')
axs[0].set_title('Non-fiber Particles')

## 
plt.figure(6).clf()
fig,axs=plt.subplots(4,1,num=6,sharex=True)

for ax,pathway in zip(axs,
                      ['effluent', 'stormwater', 'manta']):
    sel=(combined.pathway==pathway)
    # sel=sel & np.isfinite(combined.w_s) & (combined.w_s<0)
    sel=sel & (combined.Category_Final!='Fiber')
    values=combined.loc[sel,'Width.mm']
    values=values[np.isfinite(values)]
    values=values[(values<10) & (values>0)]
    values=np.log10(values)
    
    sns.distplot(values,ax=ax,bins=np.linspace(-2.5,1,50))
    
    ax.text(0.01,0.9,cap(pathway),transform=ax.transAxes,va='top')
    if ax!=axs[-1]:
        ax.set_xlabel('')

    ax.axvline(np.log10(0.335),color='k',lw=0.5)

    sns.kdeplot(values,ax=axs[-1],label=cap(pathway))
    
axs[-1].xaxis.set_major_formatter( ticker.FuncFormatter(lambda x,i: f"{10**x:.3f}"))
axs[-1].set_xlabel('Width (mm)')
axs[0].set_title('Non-fiber Particles')

##

# Assume
#  - majority of particles are from stormwater
#  - the actual abundance is "independent" of size.  this is a terrible
#    assumption, but otherwise not sure how to back out the differences.
#  - specifically, that the distributions above 0.4mm should be the same

for dim in ['Length.mm','Width.mm']:
    for pathway in ['effluent', 'stormwater', 'manta']:
        sel=(combined.pathway==pathway)
        sel=sel & (combined.Category_Final!='Fiber')
        values=combined.loc[sel,dim]
        values=values[np.isfinite(values)]

        frac_400 = (values>0.400).sum() / len(values)
        print(f"{pathway:12}: {dim:9}  above 0.4mm  {frac_400:.3f}   below 0.4mm {1-frac_400:.3f}")

# So say the length is the effective parameter
# stormwater has 0.564 particles above 400
# manta has 0.900.
# then manta is missing x particles below 400..

## 
#mlarge/(mlarge+msmall) = 0.900
#mlarge/(mlarge+msmall+mlost) = 0.564
mlost_l=0.900/0.564 - 1 # 0.596

# or by width: 2.3
mlost_w=0.665/0.201 - 1

# By length 


