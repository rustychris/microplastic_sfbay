# See Tokai et al 2021

# logistic parameters for 1mm net, giving catch efficiency (retention
# probability) as a function of longest dimension, excluding fibers.
a=-7.27
b=3.67

def rl_1mm(l,m=1.0):
    return np.exp(a+b*l/m)/(1+np.exp(a+b*l/m))


l=np.linspace(0,6,200)

import matplotlib.pyplot as plt
plt.figure(1).clf()

plt.plot(l,rl_1mm(l,1.0),label="1mm")
plt.plot(l,rl_1mm(l,0.335),label="335 $\mu$m")
plt.plot(l,rl_1mm(l,0.125),label="125 $\mu$m")

plt.legend(loc='lower right') 
plt.xlabel('Longest axis length (mm)')
plt.ylabel('Retention probability')
plt.savefig('tokai-retention-probability.png')

##

# What does this look like when summed over the stormwater
# size dist?

# Stormwater and wastewater both used 355um and 125um sieves.
from stompy import utils
utils.path("../field_data")
import plastic_data

##

# Stormwater samples
#  - distribution over the two sieve sizes
#  - Overall distribution of size, and implication when compared to catch efficiency
#    curves

# SizeFraction, SizeFraction_Final,
# Length.mm

# Ignore blanks
storm_field=plastic_data.storm_df[ plastic_data.storm_df.field_sample_p ]

# storm_field.groupby('SizeFraction_Final').size()
#  SizeFraction_Final
#  125um-355um    8326
#  355um-500um    1724
#  500um-1mm      1505
#  >1mm            801

# Catch efficiency of 355, relative to 125:
# (1724 + 1505 + 801) / (1724 + 1505 + 801 + 8326)
# = 0.326
#  
effluent_field=plastic_data.effluent_df[ plastic_data.effluent_df.field_sample_p ]
# effluent_field.groupby('SizeFraction_Final').size()
# SizeFraction_Final
# 125    1684
# 355    1702
# Relative catch efficiency 0.503

##
storm=storm_field.copy()
# Stormwater, compare with expected curves
storm['size_class']=np.where( storm['SizeFraction_Final']=='125um-355um',
                              125,355 )
# 216 have zero length, and many have no length at all
# I'm losing 7681 particles here due to nan length
#              isnan
#  size_class  Length.mm
#  125         False        2278
#              True         6048
#  355         False        2403
#              True         1633
# so almost 75% of the 125um size class is missing length,
# while only 33% of the 355um size class is missing length.


sel=storm['Length.mm']>0.0
storm_with_L=storm[ storm['Length.mm']>0.0 ]
# 4465 remain

plt.figure(2).clf()

fig,axs=plt.subplots(2,1,num=2,sharex=True)
# Reasonably log-normal
#axs[0].hist(np.log10(storm_lengths),bins=40,density=True)
from scipy.stats import gaussian_kde
log_lengths=np.linspace(-1.5,2.0)

# So fit gaussian kde to log of length
kernel_all=kde.gaussian_kde(np.log10(storm['Length.mm']))
kernel_355=kde.gaussian_kde(np.log10(storm['Length.mm'][storm.size_class==355]))
kernel_125=kde.gaussian_kde(np.log10(storm['Length.mm'][storm.size_class==125]))
Nall=len(storm['Length.mm'])
N355=(storm['size_class']==355).sum()
N125=(storm['size_class']==125).sum()

axs[0].plot(log_lengths, kernel_all(log_lengths),label='KDE all storm')
axs[0].plot(log_lengths, kernel_355(log_lengths),label='KDE 355um storm')
axs[0].plot(log_lengths, kernel_125(log_lengths),label='KDE 125um storm')

axs[0].legend()

# Compare the predicted ratio and the observed ratio
# Not quite right. these have all been normalized to be densities
# kernel
axs[1].plot(log_lengths, kernel_355(log_lengths)/kernel_all(log_lengths), label='Observed rel. $\eta$')
