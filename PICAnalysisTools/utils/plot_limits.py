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
    # The name here is confusing. I made this for plotting on a log scale but I don't exclusivley use it for this.
    
    if np.min(array) == 0:
        plt_min = 1                                 # order of magnitude of min
    elif np.min(array) < 0:
        order_min = np.floor(np.log10(np.max(-1*array))) # order of magnitude of max
        plt_min   = -1*roundup(np.max(-1*array), 10**(order_min+min_offset))
    else:
        order_min = np.floor(np.log10(np.min(array))) # order of magnitude of min
        plt_min   = rounddn(np.min(array),  10**(order_min+min_offset))
        
    order_max = np.floor(np.log10(np.max(array))) # order of magnitude of max
    plt_max   = roundup(np.max(array), 10**(order_max+max_offset))
    
    return plt_min, plt_max


def plt_limits_log_absolute(array, offset = 0):

    order = np.floor(np.log10(np.max(abs(array)))) # order of magnitude of max
    plt_max   = roundup(np.max(abs(array)), 10**(order+offset))

    return plt_max