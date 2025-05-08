# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 09:51:39 2021

@author: lrreid

define colormap with white background
"""

from matplotlib.colors import ListedColormap
from numpy import asarray, linspace
from matplotlib.pyplot import get_cmap

def wbc(cmap , thresh , thresh_position = 'lower'):

    # get colormap and make array
    color_map = get_cmap(cmap , 100)
    cmap_array = asarray(color_map(linspace(0, 1, 100)))

    if thresh_position == 'upper':
        for i in range(thresh):
            cmap_array[99 - i , 3] = float(i) / thresh
    else:
        for i in range(thresh):
            cmap_array[i , 3] = float(i) / thresh

    return ListedColormap(cmap_array)


def cmap_white(cmap):
    # get colormap and make array
    color_map  = get_cmap(cmap , 255)
    cmap_array = asarray(color_map(linspace(0, 1, 255)))

    cmap_array[0 , :] = 0

    return ListedColormap(cmap_array)
