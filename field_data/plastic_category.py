import six
import pandas as pd
import plastic_data

##
six.moves.reload_module(plastic_data)
from plastic_data import effluent_df,storm_df,combined

types=combined.groupby(['pathway','Category_Final'])['pathway'].count()

types=types.unstack(0,fill_value=0)

types['combined']=types['stormwater'] + types['effluent']
types=types.sort_values('combined',ascending=False)

table=types.to_string()

fig=plt.figure(1)
fig.clf()
fig.text( 0.01, 0.98, table, va='top', family='monospace' )

fig.savefig("plastic-category.pdf")
fig.savefig("plastic-category.png")

##

