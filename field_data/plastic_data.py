import os
import pandas as pd
import calc_settling
import numpy as np

##
try:
    cwd=os.path.dirname(__file__)
except NameError:
    cwd="." # assume interactive run is from this directory.

effluent_df=pd.read_csv(os.path.join(cwd,"2019-04-08_Datafile_Effluent.csv"))
effluent_df['pathway']='effluent'
storm_df=pd.read_csv(os.path.join(cwd,"2019-04-08_Data_Stormwater.csv"))
storm_df['pathway']='stormwater'

for df in [effluent_df,storm_df]:
    w_s=[calc_settling.record_to_ws(row)
         for i,row in df.iterrows()]
    w_s=np.array(w_s)
    df['w_s']=w_s

combined=pd.concat([storm_df,effluent_df],sort=False)

