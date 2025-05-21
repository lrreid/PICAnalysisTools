'''
Statistical calculations

TO DO:
    - Create 1D version of D4S centroid calculation.

'''


import numpy as np

def w_std( a, weights ):                                # Define weighted standard deviation calculation
    average = np.average(a, weights=weights)
    variance = np.average((a - average) ** 2, weights=weights)
    return( average, np.sqrt(variance) )


def D4S_centroid(image): # Find centroid of 2D array using D4-sigma method
    v, h = image.shape

    p = np.sum(image)                       # float avoids integer overflow
    # find the centroid
    hh = np.arange(h, dtype=float)          # float avoids integer overflow
    vv = np.arange(v, dtype=float)          # ditto
    xc = sum(np.dot(image, hh))/p
    yc = sum(np.dot(image.T, vv))/p
    
    return(xc, yc)


def find_nearest(array, target):
    # What does this function do if it finds two equally good candidates?
    array = np.asarray(array)
    idx = (np.abs(array - target)).argmin()
    return idx, array[idx]