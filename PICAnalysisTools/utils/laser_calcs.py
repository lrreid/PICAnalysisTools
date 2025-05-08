"""

Simple laser calculations and conversions

"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.unit_conversion import magnitude_conversion, magnitude_conversion_area


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
    Int_SI     = magnitude_conversion_area(Int, int_unit, "")
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
    
    return magnitude_conversion_area(Int, "", int_unit)


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

    return a0, magnitude_conversion_area(Int, "", int_unit)

def ponderomotive_energy(a0, energy_unit="kilo"):

    E_pond_max = (m_e*c**2)*(np.sqrt(1+(a0**2))-1)
    E_pond_avg = (m_e*c**2)*(np.sqrt(1+((a0/2)**2))-1)

    return magnitude_conversion(E_pond_max, "", energy_unit), magnitude_conversion(E_pond_avg, "", energy_unit)