"""
Define standard distributions useful for laser and plasma simulations and calculations.
Super-Gauss derivations for conversions between FWHM, 1/e^2 etc. in bound notebook, 12/04/2023.
"""

import numpy as np


def Gauss_1D_sigma(Amp, centroid, sigma, x_data):
    Gauss  = Amp * np.exp( -1 * ( ((x_data-centroid)**2)/(2*(sigma**2))  )    )
    FWHM   = (2*np.sqrt(2*np.log(2))*sigma)
    return Gauss, FWHM

def Gauss_1D_FWHM(Amp, centroid, FWHM, x_data):
    sigma  = FWHM/(2*np.sqrt(2*np.log(2)))
    Gauss  = Amp * np.exp( (-4*np.log(2)) * ( ((x_data-centroid)**2)/(FWHM**2) )    )
    return Gauss, sigma

def superGauss_1D(Amp, x0, sigma, order, x_data):
    Super_Gauss_1D = Amp * np.exp(-1*( ((x_data-x0)**2) / (2*sigma**2) )**order )

    rms_radius = np.sqrt(2)*sigma*(np.log(np.sqrt(np.exp(1))))**(1/(2*order))
    FWHM       = 2*np.sqrt(2)*sigma*(np.log(1/0.5))**(1/(2*order))
    e2_radius  = np.sqrt(2)*sigma*(np.log(np.exp(1)**2))**(1/(2*order))
    
    return Super_Gauss_1D, rms_radius, FWHM, e2_radius

def superGauss_1D_FWHM(Amp, x0, FWHM, order, x_data):
    Super_Gauss_1D = Amp * np.exp(-np.log(2)*( ((4*((x_data-x0)**2)) / (FWHM**2))  )**order )
    e2_radius      = (FWHM/2)*( (np.log( np.exp(1)**2)/np.log(2))**(1/(2*order)) ) 
    rms_radius     = (FWHM/2)*( (np.log( np.sqrt(np.exp(1)))/np.log(2))**(1/(2*order)) ) 

    return Super_Gauss_1D, rms_radius, e2_radius

def superGauss_width_conversion_sigma(order, sigma, norm_int):

    radius = np.sqrt(2)*sigma*(np.log(1/norm_int))**(1/(2*order))

    return radius

def superGauss_width_conversion_FWHM(order, FWHM, norm_int):

    radius = (FWHM/2)*( (np.log( 1/norm_int)/np.log(2))**(1/(2*order)) ) 

    return radius

