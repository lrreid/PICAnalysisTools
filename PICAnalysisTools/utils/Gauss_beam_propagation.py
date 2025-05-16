"""
Functions related to propagation of Gaussian laser beams
"""


import numpy as np
from scipy.constants import pi
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion


def rayleigh_length(w0, lambda0, spot_unit = "micro", wavelength_unit = "nano", rayleigh_unit = "milli"):

    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")

    zr = (pi*(w0_SI**2))/lambda0_SI    

    return magnitude_conversion(zr, "", rayleigh_unit)


def Gaussian_propagation(w0, lambda0, z, spot_unit = "micro", wavelength_unit = "nano", z_unit = "milli"):

    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")
    z_SI        = magnitude_conversion(z, z_unit, "")

    w = w0_SI*np.sqrt(1+((lambda0_SI*z_SI)/(pi*w0_SI**2))**2)

    return magnitude_conversion(w, "", spot_unit)


def prop_for_beam_size(w0, lambda0, beam_radius, spot_unit = "micro", wavelength_unit = "nano", beam_unit = "micro", z_unit = "milli"):

    w0_SI          = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI     = magnitude_conversion(lambda0, wavelength_unit, "")
    beam_radius_SI = magnitude_conversion(beam_radius, beam_unit, "")

    if w0 == beam_radius:
        Z = 0
    else:
        Z = ((pi*w0_SI**2)/lambda0_SI) * np.sqrt((beam_radius_SI/w0_SI)**2-1)
    
    return magnitude_conversion(Z, "", z_unit)


def waist_from_focal_length(focal_length, lambda0, beam_radius, f_unit = "milli", wavelength_unit = "nano", beam_unit = "micro", spot_unit = "micro", rayleigh_unit = "milli"):

    focal_length_SI = magnitude_conversion(focal_length, f_unit, "")
    lambda0_SI      = magnitude_conversion(lambda0, wavelength_unit, "")
    beam_radius_SI  = magnitude_conversion(beam_radius, beam_unit, "")

    w0 = (focal_length_SI*lambda0_SI)/(pi*beam_radius_SI)
    zr = (pi*(w0**2))/lambda0_SI

    return magnitude_conversion(w0, "", spot_unit), magnitude_conversion(zr, "", rayleigh_unit)


def focal_length_from_waist(lambda0, beam_radius, w0, wavelength_unit = "nano", beam_unit = "micro", spot_unit = "micro", f_unit = "milli"):

    lambda0_SI     = magnitude_conversion(lambda0, wavelength_unit, "")
    beam_radius_SI = magnitude_conversion(beam_radius, beam_unit, "")
    w0_SI          = magnitude_conversion(w0, spot_unit, "")

    focal_length = (pi*beam_radius_SI*w0_SI)/lambda0_SI

    return magnitude_conversion(focal_length, "", f_unit)