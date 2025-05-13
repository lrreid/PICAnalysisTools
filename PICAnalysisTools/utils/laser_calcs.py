"""

Simple laser calculations and conversions

To Do:
    - Add functions for:
        propagation of Gaussian beams
        spot size at distance z from waist
        focusing of top hat beams
    - Tidy up all the functions
    - Add doc strings


"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e, h
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_area


def a0_from_intensity(Int, lambda0=800, int_unit="centi", wavelength_unit="nano"):
    """
    
    Parameters
    ----------
    Int : float
        Laser peak intensity (Wcm-2).
    lambda0 : float, optional
        Central wavelength of laser light (nm). The default is 800e-9.

    Returns
    -------
    a0 : float
        Normalised laser vector potential.

    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")
    Int_SI     = magnitude_conversion_area(Int, int_unit, "", reciprocal_units = True)
    a0         = ((e*lambda0_SI)/(pi*m_e)) * np.sqrt((Int_SI)/(2*epsilon_0*c**5))           # laser a0
    
    return a0


def intensity_from_a0(a0, lambda0=800, wavelength_unit="nano", int_unit="centi"):
    """
    
    Parameters
    ----------
    a0 : float
        Normalised laser vector potential.
    lambda0 : float, optional
        Central wavelength of laser light (nm). The default is 800e-9.

    Returns
    -------
    Int : float
        Laser peak intensity (Wcm-2).
    """
    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")

    Int = (2*a0**2*pi**2*epsilon_0*m_e**2*c**5)/(e**2*lambda0_SI**2)
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True)


def Gaussian_laser_intensitiy(Energy, tau_FWHM, w0, lambda0=800, energy_unit="milli", time_unit = "femto", spot_unit="micro", wavelength_unit="nano", int_unit="centi"):
    """
    Parameters
    ----------
    Energy : float
        Laser energy (J).
    tau_FWHM : float
        FWHM laser pulse duration (s).
    w0 : float
        laser spot size (1/e2 radius) (m).
    lambda0 : float
        laser central wavelength (m).

    Returns
    -------
    a0 : float
        Normalised laser vector potential.
    Int : float
        Laser peak intensity (Wcm^-2).

    """

    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")

    Int = (2*Energy_SI)/(tau_FWHM_SI*pi*w0_SI**2)                                                # Intensity (W/m^2). Gaussian profile assumed.
    a0  = ((e*lambda0_SI)/(pi*m_e*c**2*w0_SI)) * np.sqrt( Energy_SI/(pi*epsilon_0*c*tau_FWHM_SI) )  # Laser amplitude. Gaussian profile assumed.

    return a0, magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True)


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

    # Add options for units.
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


class laser_wavelength_conversions():

    def __init__(self, lambda0, wavelength_unit="nano", freq_unit="Tera", wavenumber_unit = "centi", time_unit = "femto", eV_unit = "", Joule_unit = "" ):

        self.lambda0          = lambda0
        self.wavelength_unit  = wavelength_unit
        self.freq_unit        = freq_unit
        self.wavenumber_unit  = wavenumber_unit
        self.time_unit        = time_unit
        self.eV_unit          = eV_unit
        self.Joule_unit       = Joule_unit
        self.lambda0_SI       = magnitude_conversion(self.lambda0, self.wavelength_unit, "", reciprocal_units = False)
        self.freq, self.omega = self.get_frequency()
        self.wavenumber       = self.get_wavenumber()
        self.period           = self.get_period()
        self.Eph              = self.get_energy()
    

    def get_frequency(self):

        freq  = c/self.lambda0_SI
        omega = (2*pi*c)/self.lambda0_SI

        return magnitude_conversion(freq, "", self.freq_unit), magnitude_conversion(omega, "", self.time_unit, reciprocal_units = True)
    
    def get_wavenumber(self):

        k = 1/self.lambda0_SI

        return magnitude_conversion(k, "", self.wavenumber_unit, reciprocal_units = True)
    

    def get_period(self):

        T = self.lambda0_SI/c

        return magnitude_conversion(T, "", self.time_unit)
    
    def get_energy(self):

        Eph_eV = (h*c)/(self.lambda0*e)
        Eph_J  = (h*c)/(self.lambda0)

        return magnitude_conversion(Eph_eV, "", self.eV_unit), magnitude_conversion(Eph_J, "", self.Joule_unit)

def ponderomotive_energy(a0, energy_unit="kilo"):

    E_pond_max = (m_e*c**2)*(np.sqrt(1+(a0**2))-1)
    E_pond_avg = (m_e*c**2)*(np.sqrt(1+((a0/2)**2))-1)

    return magnitude_conversion(E_pond_max, "", energy_unit), magnitude_conversion(E_pond_avg, "", energy_unit)