"""

Simple laser calculations and conversions

To Do:
    - Add function for focusing of top hat beams
    - Tidy up all the functions
    - Add doc strings


"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e, h
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_area

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

#%% Focusing of top hat beams

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


#%% Focusing of top hat beams

def top_hat_laser_intensitiy(Energy, tau_FWHM, beam_rad, lambda0, focal_length, energy_unit="milli", time_unit = "femto", beam_unit = "milli", f_unit="milli", wavelength_unit="nano", int_unit="centi", spot_unit="micro"):
    #Eqn from Kristjan Poder PhD thesis (2016), page 248. See also: Alex Picksley Thesis page 53.
    #Eqn for Airy disk (Bessel J1 from Wikipedia & Hecht Optics)
    # https://en.wikipedia.org/wiki/Airy_disk

    Energy_SI       = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI     = magnitude_conversion(tau_FWHM, time_unit, "")
    beam_rad_SI     = magnitude_conversion(beam_rad, beam_unit, "")
    lambda0_SI      = magnitude_conversion(lambda0, wavelength_unit, "")
    focal_length_SI = magnitude_conversion(focal_length, f_unit, "")
    
    omega_0 = 0.823*lambda0_SI*(focal_length_SI/(2*beam_rad_SI))
    Int     = (Energy_SI*pi*beam_rad_SI**2)/(tau_FWHM_SI*(lambda0_SI**2)*(focal_length_SI**2))    
    a0      = a0_from_intensity(Int, lambda0=lambda0_SI, int_unit="", wavelength_unit="")
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0, magnitude_conversion(omega_0, "", spot_unit)

def top_hat_laser_intensity_from_spot(Energy, tau_FWHM, lambda0, omega_0, energy_unit="milli", time_unit = "femto", wavelength_unit="nano", spot_unit = "micro", int_unit="centi"):
    #Eqn from Kristjan Poder PhD thesis (2016), page 248. See also: Alex Picksley Thesis page 53.
    #Eqn for Airy disk (Bessel J1 from Wikipedia & Hecht Optics)

    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")
    omega_0_SI  = magnitude_conversion(omega_0, spot_unit, "")
    
    Int = ( (0.823**2)*Energy_SI*pi ) / ( 4*tau_FWHM_SI*(omega_0_SI)**2  )
    a0  = a0_from_intensity(Int, lambda0=lambda0_SI, int_unit="", wavelength_unit="")
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0

#%% Other calculations


def wavelength_from_photon_energy():

    return None

def ponderomotive_energy(a0, energy_unit="kilo"):

    E_pond_max = (m_e*c**2)*(np.sqrt(1+(a0**2))-1)
    E_pond_avg = (m_e*c**2)*(np.sqrt(1+((a0/2)**2))-1)

    return magnitude_conversion(E_pond_max, "", energy_unit), magnitude_conversion(E_pond_avg, "", energy_unit)

def photons_per_pulse():

    return None


def ND_filter_transmission(OD):

    Transmission = 10**(-OD)

    return Transmission

def ND_filter_Energy_transmission(OD, Energy_in, energy_in_unit : str = "milli", energy_out_unit : str = "micro"):

    Energy_in_SI        = magnitude_conversion(Energy_in, energy_in_unit, "")
    filter_transmission = ND_filter_transmission(OD)
    Energy_out          = Energy_in_SI*filter_transmission

    return magnitude_conversion(Energy_out, "", energy_out_unit)