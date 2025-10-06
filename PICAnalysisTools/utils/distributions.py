"""
Define standard distributions useful for laser and plasma simulations and calculations.
Super-Gauss derivations for conversions between FWHM, 1/e^2 etc. in bound notebook, 12/04/2023.

To Do:
    - add docstrings
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

def superGauss_beam_conversions_FWHM(order, FWHM, show_prints:bool = False):

    rms_radius     = (FWHM/2)*( (np.log( np.sqrt(np.exp(1)))/np.log(2))**(1/(2*order)) ) 
    e2_radius      = (FWHM/2)*( (np.log( np.exp(1)**2)/np.log(2))**(1/(2*order)) ) 
    one_pc_contour = (FWHM/2)*( (np.log( 1/0.01  )/np.log(2))**(1/(2*order)) ) 
    clear_aperture = (FWHM/2)*( (np.log( 1/0.001 )/np.log(2))**(1/(2*order)) ) 

    if show_prints is True:
        print("For order %d and FWHM: %d" % (order, FWHM) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("1/e2 dia: %0.2f" % np.round(e2_radius*2, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )
        print("Clear aperture dia: %0.2f" % np.round(clear_aperture*2, 2) )

    return rms_radius, e2_radius, one_pc_contour, clear_aperture


def superGauss_beam_conversions_e2(order, e2_radius, show_prints:bool = False):

    rms_radius     = (e2_radius)*( (np.log(np.sqrt(np.exp(1)))/np.log(np.exp(1)**2))**(1/(2*order)) )
    FWHM           = (2*e2_radius)*( (np.log(2)/np.log(np.exp(1)**2))**(1/(2*order)) )
    one_pc_contour = (e2_radius)*( (np.log(1/0.01 )/np.log(np.exp(1)**2))**(1/(2*order)) )
    clear_aperture = (e2_radius)*( (np.log(1/0.001)/np.log(np.exp(1)**2))**(1/(2*order)) )

    if show_prints is True:
        print("For order %d and 1/e^2 radius: %d" % (order, e2_radius) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("FWHM: %0.2f" % np.round(FWHM, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )
        print("Clear aperture dia: %0.2f" % np.round(clear_aperture*2, 2) )

    return rms_radius, FWHM, one_pc_contour, clear_aperture


def superGauss_beam_conversions_CA(order, clear_aperture_rad, show_prints:bool = False):
    # Clear aperture is the 0.1% intenstiy contour

    rms_radius     = (clear_aperture_rad)*( (np.log(np.sqrt(np.exp(1)))/np.log(1/0.001))**(1/(2*order)) )
    FWHM           = (2*clear_aperture_rad)*( (np.log(2)/np.log(1/0.001))**(1/(2*order)) )
    e2_radius      = (clear_aperture_rad)*( (np.log(np.exp(1)**2)/np.log(1/0.001))**(1/(2*order)) )
    one_pc_contour = (clear_aperture_rad)*( (np.log(1/0.01)/np.log(1/0.001))**(1/(2*order)) )

    if show_prints is True:
        print("For order %d and clear aperture radius: %d" % (order, clear_aperture_rad) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("FWHM: %0.2f" % np.round(FWHM, 2) )
        print("1/e2 dia: %0.2f" % np.round(e2_radius*2, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )

    return None