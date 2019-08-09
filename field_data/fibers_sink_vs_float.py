import pandas as pd
from plastic_data import effluent_df,storm_df,combined
import matplotlib.pyplot as plt
import numpy as np

##
def disp_sel(sel,label):
    ws=combined['w_s'][sel]
    ws_median=np.nanmedian(ws)
    
    print(f"{label}: {sel.sum()} samples {np.isfinite(ws).sum()} valid ws, median {ws_median}")

disp_sel( (combined['PlasticType_Final']=='Cellulose acetate') & (combined['Category_Final']=='Fiber'),
          "Cellulose acetate fibers")
          
disp_sel( (combined['PlasticType_Final']=='Polyester') & (combined['Category_Final']=='Fiber'),
          "Polyester fibers")


