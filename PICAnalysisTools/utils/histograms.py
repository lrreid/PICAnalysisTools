"""
Functions for creating histograms.

These functions have not been tested and are not implemented into other functions.

TO DO:
    - Add option for no weights. 
    - Write doc strings

"""

# import numpy as np
from numpy import histogram, histogram2d
from PICAnalysisTools.utils import binning


def get_histogram_1d(data, bin_res, data_round, weights, equal_bins = False):

    if equal_bins == False:
        Bin_Min, Bin_Max, Bins = binning.get_bins(data, bin_res, data_round)
    elif equal_bins == True: # makes Bin_Min = -Bin_Max
         Bin_Min, Bin_Max, Bins = binning.get_bins_absolute(data, bin_res, data_round)

    hist_1d = histogram(data, bins=Bins, weights=weights)

    return hist_1d, Bin_Min, Bin_Max, Bins


def get_histogram_2d(data_x, bin_res_x, data_round_x, data_y, bin_res_y, data_round_y, weights, equal_bins_x = False, equal_bins_y = False):

    if equal_bins_x == False:
        Bin_Min_x, Bin_Max_x, Bins_x = binning.get_bins(data_x, bin_res_x, data_round_x)
    elif equal_bins_x == True: # makes Bin_Min = -Bin_Max
         Bin_Min_x, Bin_Max_x, Bins_x = binning.get_bins_absolute(data_x, bin_res_x, data_round_x)

    if equal_bins_y == False:
        Bin_Min_y, Bin_Max_y, Bins_y = binning.get_bins(data_y, bin_res_y, data_round_y)
    elif equal_bins_y == True: # makes Bin_Min = -Bin_Max
         Bin_Min_y, Bin_Max_y, Bins_y = binning.get_bins_absolute(data_y, bin_res_y, data_round_y)

    hist_2d, _, _ = histogram2d(data_x, data_y, bins=(Bins_x, Bins_y), weights=weights )

    return hist_2d, Bin_Min_x, Bin_Max_x, Bins_x, Bin_Min_y, Bin_Max_y, Bins_y



