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


def set_bold_labels(ax,y=-0.06):
    """
    Common styling for x axis in bi-log plots.
    Puts a bold 'Sinks', 'Passive', 'Floats'
    on the xaxis, and hides the middle tick label
    """
    txtargs=dict(style='italic',fontweight='bold',fontsize=12)
    ax.text(0.0,y,"Floats",transform=ax.transAxes,va='top',ha='left',**txtargs)
    ax.text(1,y,"Sinks",transform=ax.transAxes,va='top',ha='right',**txtargs)
    ax.text(0.5,y,"Passive",transform=ax.transAxes,va='top',ha='center',**txtargs)
    nticks=len(ax.get_xticklabels())
    ax.get_xticklabels()[nticks//2].set_visible(0)
    #ax.xaxis.set_visible(0)
