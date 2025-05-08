'''
particle selection

To Do
    - Add option to choose units?
'''


import numpy as np


def radial_selection(R_Select_min, R_Select_max, x, y, z, ux, uy, uz, w):
    R = np.sqrt(x**2 + y**2)                                            # Find distance for origin of each macroparticle (m)

    arg_max_radial = np.argwhere(R > R_Select_max)                      # Get indices of all macroparticles > r_max
    arg_min_radial = np.argwhere(R < R_Select_min)                      # Get indices of all macroparticles > r_max

    indices_raidal = np.concatenate((arg_max_radial, arg_min_radial))   # Combine lists of indices to remove from data set

    # delete all macroparticles outwith selection
    x_radial_select = np.delete(x, indices_raidal)
    y_radial_select = np.delete(y, indices_raidal)
    z_radial_select = np.delete(z, indices_raidal)

    ux_radial_select = np.delete(ux, indices_raidal)
    uy_radial_select = np.delete(uy, indices_raidal)
    uz_radial_select = np.delete(uz, indices_raidal)
    w_radial_select  = np.delete(w,  indices_raidal)

    return x_radial_select, y_radial_select, z_radial_select, ux_radial_select, uy_radial_select, uz_radial_select, w_radial_select
