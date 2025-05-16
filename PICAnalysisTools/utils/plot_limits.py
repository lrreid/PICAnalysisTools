"""
Define functions for calculating plot limits
"""


import numpy as np
from PICAnalysisTools.utils.rounding import roundup, rounddn


def plt_limits(array, rounding):
    
    plt_min = rounddn(np.min(array), rounding)
    plt_max = roundup(np.max(array), rounding)
    
    return plt_min, plt_max

def plt_limits_absolute(array, rounding):
    
    plt_max = roundup( np.max(abs(array)) , rounding)

    return plt_max


def plt_limits_log(array, min_offset = 0, max_offset = 0):
    
    if np.min(array) == 0:
        plt_min = 1                                 # order of magnitude of min
    else:
        order_min = np.floor(np.log10(np.min(array))) # order of magnitude of min
        plt_min   = rounddn(np.min(array),  10**(order_min+min_offset))
        
    order_max = np.floor(np.log10(np.max(array))) # order of magnitude of max
    plt_max   = roundup(np.max(array), 10**(order_max+max_offset))
    
    return plt_min, plt_max