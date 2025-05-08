
"""

Date created: 01/05/2025
Authors: Lewis R Reid
"""

import numpy as np
from rounding import roundup, rounddn

def get_bins(array, r_res, r_round):
    # array can be single arrary of numbers array = x or list of arrays array = [x, y, z]

    r_max  = roundup(np.max(array), r_round)
    r_min  = rounddn(np.min(array), r_round)
    R_Bins = np.arange(r_min, r_max+r_res, r_res)

    return r_min, r_max, R_Bins

def get_bins_absolute(array, r_res, r_round):
    # array can be single arrary of numbers array = x or list of arrays array = [x, y, z]

    r_max  = roundup(np.max(abs(array)), r_round)
    R_Bins = np.arange(-r_max, r_max+r_res, r_res)

    return -r_max, r_max, R_Bins