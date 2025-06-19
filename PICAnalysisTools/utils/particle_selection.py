'''
particle selection

To Do
    - Add option to choose units?
'''


import numpy as np


def radial_selection(R_Select_min, R_Select_max, x, y, z, ux, uy, uz, w):
    """
    Choose radial selection of particle species macroparticles according to distance from the origin in a PIC simulation with cylindrical geometery

    Parameters
    ----------
    R_Select_min : float
        Minimum radius of particles within selection (m)
    R_Select_max : float
        Maximum radius of particles within selection (m)
    x : array_like
        x-coordinates of particle positions (m)
    y : array_like
        y-coordinates of particle positions (m)
    z : array_like
        z-coordinates of particle positions (m)
    ux : array_like
        x-coordinates of particle momenta (px/mc)
    uy : array_like
        y-coordinates of particle momenta (py/mc)
    uz : array_like
        z-coordinates of particle momenta (pz/mc)
    w : array_like
        Macroparticle weights

    Returns
    -------
    x_radial_select: array_like
        x-coordinates of particle positions within radial selection limits (m)
    y_radial_select: array_like
        y-coordinates of particle positions within radial selection limits (m)
    z_radial_select: array_like
        z-coordinates of particle positions within radial selection limits (m)
    ux_radial_select: array_like
        x-coordinates of particle momenta within radial selection limits (px/mc)
    uy_radial_select: array_like
        y-coordinates of particle momenta within radial selection limits (py/mc)
    uz_radial_select: array_like
        z-coordinates of particle momenta within radial selection limits (pz/mc)
    w_radial_select: array_like
        Macroparticle weights within radial selection limits
    """

    R              = np.sqrt(x**2 + y**2)                               # Find distance for origin of each macroparticle (m)
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
