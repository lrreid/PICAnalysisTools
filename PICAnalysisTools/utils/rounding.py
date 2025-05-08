# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:42:33 2022

@author: lrreid

Rounding functions to round up and down to the nearest arbitrary number
"""

import numpy as np

def roundup(x,n):                           # Round up to nearest n
    return int(np.ceil(x / n)) * n

def rounddn(x,n):                           # Round down to nearest n
    return int(np.floor(x / n)) * n


def normalise(array):                       # Normalise array
    peak = np.max(abs(array))
    norm = array/peak
    return(norm)


def find_nearest(array, target):
    array = np.asarray(array)
    idx = (np.abs(array - target)).argmin()
    return idx, array[idx]