import six
import numpy as np
import pandas as pd
import plastic_data
from plastic_data import effluent_df,storm_df,combined

## get 1900 w_s, from 16000 samples.
w_s=combined['w_s'].values

w_s_val=w_s[np.isfinite(w_s)]

w_s_rise=w_s_val[w_s_val<0]
w_s_sink=w_s_val[w_s_val>0]

pcts=[5,10,25,50,75,90,95]
print( np.c_[pcts,
             np.percentile(w_s_rise,pcts)] )

# [[ 5.00e+00 -5.84772906e-02]
#  [ 1.00e+01 -3.94178692e-02]
#  [ 2.50e+01 -2.08655867e-02]
#  [ 5.00e+01 -6.00456526e-03]
#  [ 7.50e+01 -2.09302333e-03]
#  [ 9.00e+01 -8.55967114e-04]
#  [ 9.50e+01 -5.70287643e-04]]

# 0.06

print( np.c_[pcts,
             np.percentile(w_s_sink,pcts)] )

# [[5.00000000e+00 1.36054781e-03]
#  [1.00000000e+01 1.73768666e-03]
#  [2.50000000e+01 2.79297779e-03]
#  [5.00000000e+01 8.03437305e-03]
#  [7.50000000e+01 1.38249748e-02]
#  [9.00000000e+01 2.23237786e-02]
#  [9.50000000e+01 3.13591741e-02]]

# assuming it's "nicer" to have symmetric
# velocities.
# 0.06 m/s
# 0.02 m/s
# 0.006 m/s
# 0.002 m/s
# 0.0006 m/s
# plus passive
# that's a total of 11 classes
