"""
Define functions for calculating plot limits
"""


import numpy as np
from rounding import roundup, rounddn


def plt_limits(array, rounding):
    
    plt_min = rounddn(np.min(array), rounding)
    plt_max = roundup(np.max(array), rounding)
    
    return plt_min, plt_max

def plt_limits_absolute(array, rounding):
    
    plt_max = roundup( np.max(abs(array)) , rounding)

    return plt_max


def plt_limits_log(array):
    
    if np.min(array) == 0:
        plt_min = 0                                 # order of magnitude of min
    else:
        order_min = np.floor(np.log10(np.min(array))) # order of magnitude of min
        # plt_min   = 10**(order_min)
        plt_min = rounddn(np.min(array),  10**(order_min+0))
        
    order_max = np.floor(np.log10(np.max(array))) # order of magnitude of max

    #plt_max = 10**(order_max+1)
    plt_max = roundup(np.max(array), 10**(order_max-1))
    
    return plt_min, plt_max