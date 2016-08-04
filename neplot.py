# Some tools for plotting

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

class PowerNormalize(colors.Normalize):
    def __init__(self, clr_gamma, vmin=None, vmax=None, clip=False):
        self.clr_gamma=clr_gamma
        colors.Normalize.__init__(self,vmin,vmax,clip)
    def __call__(self,value,clip=None):
        return np.ma.masked_array(((value-self.vmin)
                /(self.vmax-self.vmin))**self.clr_gamma)
