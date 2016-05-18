import numpy as np

def lextrem(a,comp='max'):
    """
    Output local extrema. Regions of equal values with different values
    surrounding are not considered extrema! Values on the edges of the array are
    not considered extrema either.
    """
    cmpfunc={'max':np.greater,'min':np.less}[comp]
    mab=np.r_[False,cmpfunc(a[1:],a[:-1])]&np.r_[cmpfunc(a[:-1],a[1:]),False]
    mai=np.transpose(np.nonzero(mab))
    ma=a[mai]
    return (ma,mai)
