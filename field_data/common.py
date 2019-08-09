import numpy as np

default_small=5e-4

def map_bilog(x,small=default_small):
    return np.sign(x)*np.log10(np.abs(x).clip(small,np.inf)/small)

def unmap_bilog(y,small=default_small):
    return np.sign(y)*small*10**np.abs(y)

from matplotlib.ticker import ScalarFormatter
sf=ScalarFormatter()
sf._usetex=True
def mapped_label(y):
    if y==0.0:
        # s="\pm " + sf.format_data(small)
        s="Passive"
    else:
        s=sf.format_data(unmap_bilog(y))
    return "$"+s+"$"
