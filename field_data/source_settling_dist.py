import pandas as pd
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np

w_s=combined['w_s']

w_s_val=w_s[np.isfinite(w_s)]

plt.figure(1).clf()
fig,axs=plt.subplots(2,1,sharex=True,num=1)
fig.set_size_inches([6.4,4.8],forward=True)

bins=10**np.linspace(-3,-0.5,50)
axs[0].hist(w_s_val[w_s_val>0.0],bins=bins,color='#1f77b4')
axs[1].hist(-w_s_val[w_s_val<0.0],bins=bins,color='#d62728')

axs[0].axis(xmin=bins[0],xmax=bins[-1])
axs[0].set_xscale('log')
axs[1].set_xlabel('Velocity (m/s)')

axs[0].text(0.02,0.95,"Sinking",va='top',transform=axs[0].transAxes)
axs[1].text(0.02,0.95,"Rising", va='top',transform=axs[1].transAxes)

fig.savefig("source_settling_dist.png",dpi=150)

##

