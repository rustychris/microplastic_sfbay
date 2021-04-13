"""
Using updated data with weathered particles from Waldschlager, Born, Cowger, Gray, Schuttrump

Estimate some adjustment factors to apply to the previous estimates.

The paper suggests that films have the greatest deviation
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
## 

df=pd.read_csv("ws_data_2021.csv")

# calc_settling gets a record, and wants it to have:
#   'Category_Final': Foam,
cat_map={'Foam':'Foam',
         'Film':'Film',
         'Fragment':'Fragment',
         'Pellet':'Sphere'}

df['Category_Final']=df['Shape'].map(cat_map) # 'Foam', 'Film', 'Fragment', 'Pellet'

# data already has density, so pass that in and
df=df.rename({'Density  [g/cm³]':'density'},axis=1)
df['Length.mm']=df['a  [mm]']
df['Width.mm'] =df['b  [mm]']

##
import six
import calc_settling
six.moves.reload_module(calc_settling)

my_ws=[ calc_settling.record_to_ws(rec,adjust_for_weathered=False)
        for idx,rec in df.iterrows() ]

# convert to same sign convention, and m/s
real_ws= 0.01*df['Mean  Velocity '].values * np.where( df['Direc.'].values=='rising', -1, 1) 

df['my_ws']=my_ws
df['real_ws']=real_ws

##

import statsmodels.formula.api as smf

# Break that out by shape, fit linear regression

plt.figure(10).clf()
fig,axs=plt.subplots(2,2,num=10)

for i,shape in enumerate(df['Category_Final'].unique()):
    ax=axs.ravel()[i]
    
    df_cat=df[ df['Category_Final']==shape ]
    df_cat=df_cat.sort_values('my_ws')

    if shape=='Film':
        # with intercept, BIC=-180.2.
        # without intercept, BIC=-176.1
        # Use intercept.
        formula='real_ws ~ my_ws' 
    else:
        formula='real_ws ~ my_ws - 1'
    res=smf.ols(formula=formula,data=df_cat).fit()
    print(shape)
    print(res.summary())
    print()

    ax.plot( df_cat['my_ws'], df_cat['real_ws'],'g.')

    ax.plot( df_cat['my_ws'],res.predict(),'k-')

    ax.axis('equal')
    ax.axhline(0,color='k',lw=0.3)
    ax.axvline(0,color='k',lw=0.3)
    # drop dtype
    txts=[shape] + str(res.params).split("\n")[:-1]
    ax.text(0.05,0.95,"\n".join(txts),transform=ax.transAxes,
            va='top')
    ax.set_xlabel('my_ws (W&S 2019)')
    ax.set_ylabel('real ws')

fig.tight_layout()
fig.savefig('adjust_settling.png',dpi=200)

# Considering sinking and rising together
# Foam: factor of 0.41
# Film: factor of 0.64, but spread from 0.358 to 0.925 (95%)
# Fragment: factor of 0.318
# Sphere: factor of 0.4525

# Only film has good coverage on both rising and settling,
# so allow it to have an intercept.
# everybody else must go through 0.

# All set.
# Make adjustments, create new datasets, plow on...

