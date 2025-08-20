"""
Functions for calculating quantities and matching conditions for laser driven wakefield acceleration based on expressions form the literature.

TO DO:
    - Add descriptions of the assumptions/limitations of the models to docstrings.
"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e, hbar
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol
from PICAnalysisTools.utils.laser_calcs import laser_wavelength_conversions
from PICAnalysisTools.utils.plasma_calcs import PlasmaDen_Conversions, plasma_density_from_wavevector, critical_density
re       = (e**2)/(4*pi*epsilon_0*m_e*c**2)                     # classial electron radius
PC_const = (8*pi*epsilon_0*m_e**2*c**5)/(e**2)                  # critical power constant (W)
P_rel    = 0.5*PC_const                                         # Relativistic power unit (W)

def resonant_density(tau_FWHM, time_unit: str = "femto", den_unit: str = "centi"):
    """
    Calculate the resonant plasma density for laser wakefield acceleration for a given laser pulse duration.

    Parameters
    ----------
    tau_FWHM : float
        Laser pulse FWHM pulse duration. Default unit: fs
    time_unit : str, optional
        Order of magnitude of pulse duration unit, by default "femto"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e_sqrt2: float
        Resonant plasma density according to k_p*L_RMS = sqrt(2). Linear wakefield. Default unit: cm^-3
    n_e_1: float
        Resonant plasma density according to k_p*L_RMS = 1. Non-Linear wakefield. Default unit: cm^-3
    """

    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "" )

    n_e_sqrt2 = (16*np.log(2)*epsilon_0*m_e)/(e**2*tau_FWHM_SI**2)   # Plasma density (electrons.meters^-3). Density according to k_p*L_RMS = sqrt(2). Linear wakefield.
    n_e_1     = (8*np.log(2)*epsilon_0*m_e)/(e**2*tau_FWHM_SI**2)    # Plasma density (electrons.meters^-3). Density according to k_p*L_RMS = 1.       Non-linear wakefield.


    return magnitude_conversion_vol(n_e_sqrt2, "", den_unit, reciprocal_units = True), magnitude_conversion_vol(n_e_1, "", den_unit, reciprocal_units = True)


def pump_depletion_length(n_e, lambda0, a0, den_unit: str = "centi", wavelength_unit: str = "nano", length_unit: str = "milli"):
	"""
     Calculate the pump depletion length of a laser wakefield accelerator.
     
     https://doi.org/10.1103/RevModPhys.81.1229
     https://doi.org/10.1063/1.3124185


	Parameters
	----------
    n_e : float
        Plasma denstiy. Default unit: cm^-3
    lambda0 : float
        Laser central wavelength. Default unit: nm
    a0 : float
        Normalised laser vector potential.
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    length_unit : str, optional
        Order of magnitude of depletion length unit, by default "milli"

	Returns
	-------
	L_pd: float
		Pump depletion length. Default unit: mm
	"""
    
	lambda_p = PlasmaDen_Conversions(n_e=n_e, den_unit=den_unit, wavelength_unit="").lambda_p
    
	if a0**2 <= 1:
		L_pd = ((lambda_p**3)/(magnitude_conversion(lambda0, wavelength_unit, "")**2))*(2/a0**2)
	else:
		L_pd = ((lambda_p**3)/(magnitude_conversion(lambda0, wavelength_unit, "")**2))*((np.sqrt(2)*a0)/pi)
	
	return magnitude_conversion(L_pd, "", length_unit)

def dephasing_length(n_e, lambda0, a0, den_unit: str = "centi", wavelength_unit: str = "nano", length_unit: str = "milli"):
    """
    Calculate the dephasing length of a laser wakefield accelerator.
    
    https://doi.org/10.1103/RevModPhys.81.1229
    https://doi.org/10.1063/1.1842594


    Parameters
    ----------
    n_e : float
        Plasma denstiy. Default unit: cm^-3
    lambda0 : float
        Laser central wavelength. Default unit: nm
    a0 : float
        Normalised laser vector potential.
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    length_unit : str, optional
        Order of magnitude of dephasing length unit, by default "milli"

    Returns
    -------
    L_d: float
        Dephasing length of laser wakefield accelerator. Default unit: mm
    """

    lambda_p = PlasmaDen_Conversions(n_e=n_e, den_unit=den_unit, wavelength_unit="").lambda_p

    if a0**2 <= 1:
        L_d = (lambda_p**3)/(2*magnitude_conversion(lambda0, wavelength_unit, "")**2)
    else:
        L_d = ((lambda_p**3)/(2*magnitude_conversion(lambda0, wavelength_unit, "")**2))*((np.sqrt(2)*a0)/pi)

    return magnitude_conversion(L_d, "", length_unit)

def self_injection_Threshold(n_e, lambda0, den_unit: str = "centi", wavelength_unit: str = "nano"):
    """
    Threshold laser intensity for self-injection into a laser wakefield accelerator.

    Alexander Thomas model proposes an analytical approach assuming a0 >> 1 and that the plasma wave longitudinal and radial electric fields are non-evolving.
    https://doi.org/10.1063/1.3368678


    Benedetti et al obtained an expression based on particle-in-cell simulation assuming a non-evolving laser pulse which is matched to the plasma (ω0 = 2√a0 and kpLRMS ~ 1).
    https://doi.org/10.1063/1.4824811

    Parameters
    ----------
    n_e : float
        Plasma denstiy. Default unit: cm^-3
    lambda0 : float
        Laser central wavelength. Default unit: nm
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"

    Returns
    -------
    Thomas_a0: float
        Threshold laser intensity for self-injection according to Thomas model
    Benedetti_a0: float
        Threshold laser intensity for self-injection according to Benedetti model
    """

    w_p = PlasmaDen_Conversions(n_e=n_e, den_unit=den_unit).w_p
    w_0 = laser_wavelength_conversions(lambda0=lambda0, wavelength_unit=wavelength_unit).omega
    gamma_wake = w_0/w_p

    Thomas_a0    = np.log(2*gamma_wake**2)-1
    Benedetti_a0 = 2.75*np.sqrt(1 + (gamma_wake/22)**2 )

    return Thomas_a0, Benedetti_a0

def density_for_matched_spot(n_e, den_unit: str = "centi", spot_unit: str = "micro"):
    """
    Calculate the matched spot size for laser wakefield acceleration in the ?linear? regime.
    ω0 = λp/2=kp/π

    Parameters
    ----------
    n_e : float
        Plasma density. Default unit: cm^-3
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    spot_unit : str, optional
        Order of magnitude of laser spot size unit, by default "micro"

    Returns
    -------
    w0: float
        Laser spot size (to 1/e^2 radius) matched to plasma. Default unit: microns.
    """

    k_p = PlasmaDen_Conversions(n_e = n_e, den_unit = den_unit).k_p
    w0 = k_p/pi

    return magnitude_conversion(w0, "", spot_unit)


def matched_spot_from_density(w0, spot_unit: str = "micro", den_unit: str = "centi" ):
    """
    Calculate the matched plasma density for laser wakefield acceleration in the ?linear? regime from the laser spot size.
    ω0 = λp/2=kp/π

    Parameters
    ----------
    w0: float
        Laser spot size (to 1/e^2 radius) matched to plasma. Default unit: microns.
    spot_unit : str, optional
        Order of magnitude of laser spot size unit, by default "micro"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e : float
        Plasma density. Default unit: cm^-3
    """

    k_p = magnitude_conversion(w0, spot_unit, "")*pi
    n_e = plasma_density_from_wavevector(k_p = k_p, den_unit = den_unit)

    return n_e

def max_energy_scaling(a0, lambda0, tau, n_e, wavelength_unit: str = "nano", time_unit: str = "femto", den_unit: str = "centi"):
    """
    Calculate the maximum electron energy for self-injected laser wakefield acceleration from theoretical and numerical models.
    Eqns from table 1 of Lu et al. (2007). https://doi.org/10.1103/PhysRevSTAB.10.061301.
    And Gordienko & Pukhov (2005). https://doi.org/10.1063/1.1884126 
    Note these have different assumptions about the laser transverse and longitudinal size and different intensity regimes where they are valid. 

    Parameters
    ----------
    a0 : float
        Laser normalised vector potential
    lambda0 : float
        Laser central wavelength. Default unit: nm
    tau : float
        Laser FWHM pulse duration. Default unit: fs
    n_e : float
        Plasma electron density. Default unit: cm^-3
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    time_unit : str, optional
        Order of magnitude of pulse duration unit, by default "femto"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    E_max_lin: float 
        Maximum electron energy from linear theory.
    E_max_1D: float 
        Maximum electron energy from 1D nonlinear theory.
    E_max_3D: float 
        Maximum electron energy from 3D nonlinear theory.
    E_max_Pukhov: float 
        Maximum electron energy from Gordienko & Pukhov paper. 
    """

    omega_0 = laser_wavelength_conversions(lambda0, wavelength_unit = wavelength_unit).omega   # Angular frequency of laser (rad/s)
    w_p     = PlasmaDen_Conversions(n_e, den_unit=den_unit).w_p
    tau_SI  = magnitude_conversion(tau, time_unit, "")
    
    E_max_lin    = a0**2 * ((omega_0/w_p)**2) * ((m_e*c**2)/e)              # All equations from Lu paper table 1, final column 
    E_max_1D     = 4 * a0**2 * ((omega_0/w_p)**2) * ((m_e*c**2)/e)
    E_max_3D     = (2/3) * a0 * ((omega_0/w_p)**2) * ((m_e*c**2)/e)
    E_max_Pukhov = a0**(3/2) * tau_SI * ((omega_0/w_p)**2) * ((m_e*c**2)/e)       # As in Lu paper (final row in table)
    
    return E_max_lin, E_max_1D, E_max_3D, E_max_Pukhov


def max_energy_scaling_engineering(n_e, lambda0, Energy, tau_FWHM, w0, den_unit: str = "centi", wavelength_unit: str = "nano", laser_energy_unit: str = "milli", time_unit: str = "femto", spot_unit: str = "micro", elec_energy_unit: str = "Mega", power_unit: str = "tera"):
    """
    Calculations of the maximum electron energy from a laser wakefield accelerator using "engineering" models from the literature.
    Eqns from Lu et al. (2007). https://doi.org/10.1103/PhysRevSTAB.10.061301.
    And Gordienko & Pukhov (2005). https://doi.org/10.1063/1.1884126

    Parameters
    ----------
    n_e : float
        Plasma electron density. Default unit: cm^-3
    lambda0 : float
        Laser central wavelength. Default unit: nm
    Energy : float
        Energy of laser. Default unit: mJ
    tau_FWHM : float
        Pulse duration of laser to FWHM. Default unit: fs
    w0 : float
        Laser spot size to 1/e^2 radius. Default unit: microns
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    laser_energy_unit : str, optional
        Order of magnitude of laser energy unit, by default "milli"
    time_unit : str, optional
        Order of magnitude of pulse duration unit, by default "femto"
    spot_unit : str, optional
        Order of magnitude of laser spot size unit, by default "micro"
    elec_energy_unit : str, optional
        Order of magnitude of electron energy unit, by default "Mega"
    power_unit : str, optional
        Order of magnitude of laser power unit, by default "Tera"

    Returns
    -------
    E_Lu: float
        Maximum electron energy accoring to Lu model. Eqn 6. Default unit: MeV
    E_Pukhov: float
        Maximum electron energy accoring to Pukhov model. Eqn 1. Default unit: MeV
    E_PAlt: float
        Maximum electron energy accoring to Pukhov model re-rwitten by Lu. Eqn 13. Default unit: MeV
    P_laser: float
        Laser power. Default unit: TW
    """

    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    P_SI        = magnitude_conversion(Energy, laser_energy_unit, "")/tau_FWHM_SI
    P_TW        = magnitude_conversion(P_SI, "", "Tera")
    n_e_SI      = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    n_c         = critical_density(lambda0, wavelength_unit = wavelength_unit, den_unit = "" )                            # Critical density (m^-3)
    lambda_SI   = magnitude_conversion(lambda0, wavelength_unit, "")
    w0_SI       = magnitude_conversion(w0, spot_unit, "")

    E_Lu     = 1.7 * ((P_TW/100)**(1/3)) * ((1e24/n_e_SI)**(1/3)) * ((0.8e-6/lambda_SI)**(4/3))                         # Max energy according to Lu (GeV) - Lu eqn 6
    E_Pukhov = 0.65 * ((m_e*c**2)/e) * np.sqrt(P_SI/P_rel) * (c*tau_FWHM_SI/lambda_SI)                                     # Max energy according to Pukhov (eV) - Pukhov eqn 1 & Lu eqn 12
    E_PAlt   = 0.16 * ((m_e*c**2)/e) * ((c*tau_FWHM_SI)/w0_SI) * ((P_SI/P_rel)**(2/3)) * ((n_c/n_e_SI)**(1/3))             # Max energy according to Pukhov (eV) - re-rwitten by Lu (eqn 13)
    
    return magnitude_conversion(E_Lu, "Giga", elec_energy_unit), magnitude_conversion(E_Pukhov, "", elec_energy_unit), magnitude_conversion(E_PAlt, "", elec_energy_unit), magnitude_conversion(P_SI, "", power_unit)

def betatron_frequency(n_e, Ek, energy_unit: str = "Mega", den_unit: str = "centi"):
    """
    Calculate the betatron frequency of electrons oscillating in an ion channel.

    Parameters
    ----------
    n_e : float
        Plasma electron density. Default unit: cm^-3
    Ek_converted: float
        Electron energy. Default unit: MeV.
    energy_unit: string
        Order of magnitude of particle energy units, by default "mega"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    w_b: float
        Betatron (angular) frequency
    """

    elec_gamma = (magnitude_conversion(Ek, energy_unit, "")/((m_e*c**2)/e))+1
    w_p        = PlasmaDen_Conversions(n_e, den_unit=den_unit).w_p
    w_b        = w_p/np.sqrt(2*elec_gamma)

    return w_b


def betatron_wavelength(n_e, Ek, den_unit: str = "centi", energy_unit: str = "Mega", wavelength_unit: str = "micro"):
    """
    Calculate the wavelength of betatron oscaillations of electrons in an ion channel.

    Parameters
    ----------
    n_e : float
        Plasma electron density. Default unit: cm^-3
    Ek_converted: float
        Electron energy. Default unit: MeV.
    energy_unit: string
        Order of magnitude of particle energy units, by default "mega"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "micro"

    Returns
    -------
    lambda_b: float
        Wavelength of betatron oscillations
    """

    lambda_b = (2*pi*c)/betatron_frequency(n_e, Ek, energy_unit = energy_unit, den_unit = den_unit)

    return magnitude_conversion(lambda_b, "", wavelength_unit)

def betatron_critical_energy(n_e, Ek, r_beta, den_unit: str = "centi", energy_unit: str = "Mega", r_unit: str = "micro", photon_unit: str = "kilo"):
    """
    Calculate the critical energy of the betatron spectrum.

    Parameters
    ----------
    n_e : float
        Plasma electron density. Default unit: cm^-3
    Ek_converted: float
        Electron energy. Default unit: MeV.
    r_beta : float
        Radius of betatron oscillations. Default unit: microns
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    energy_unit: string
        Order of magnitude of particle energy units, by default "mega"
    r_unit : str, optional
        Order of magnitude of betatron oscillation radius unit, by default "micro"
    photon_unit : str, optional
        Order of magnitude of photon energy unit, by default "kilo"

    Returns
    -------
    E_crit: float
        Critical energy of betatron spectrum. Default unit: keV
    """

    elec_gamma = (magnitude_conversion(Ek, energy_unit, "")/((m_e*c**2)/e))+1
    w_p        = PlasmaDen_Conversions(n_e, den_unit=den_unit).w_p
    E_crit     = ((3*hbar)/(4*c*e)) * elec_gamma**2 * w_p**2 * magnitude_conversion(r_beta, r_unit, "")

    return magnitude_conversion(E_crit, "", photon_unit)


def get_electron_gamma(Ek, energy_unit: str = "Mega"):
    """
    Get the relativistic factor (gamma) of an electron from its energy

    Parameters
    ----------
    Ek_converted: float
        Electron energy. Default unit: MeV.
    energy_unit: string
        Order of magnitude of particle energy units, by default "Mega"

    Returns
    -------
    elec_gamma: float
        Relativistic factor of electrons.
    """

    elec_gamma = (magnitude_conversion(Ek, energy_unit, "")/((m_e*c**2)/e))+1

    return elec_gamma

def capillary_channel_depth(w0, spot_unit: str = "micro", den_unit: str = "centi"):
    """
    Calculate the channel depth of a parabolic radial plasma density profile for laser guiding.

    Parameters
    ----------
    w0 : float
        Laser spot radius to 1/e^2 of peak intensity. Default unit: microns
    spot_unit : str, optional
        Order of magnitude of laser spot unit, by default "micro"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    dne: float
        Channel depth of parabolic profile to guide laser. Default unit: cm^-3
    """

    w0_SI = magnitude_conversion(w0, spot_unit, "" )
    dne   = 1/(pi*re*w0_SI**2)                  # depth of plasma channel (m^-3)

    return magnitude_conversion_vol(dne, "", den_unit, reciprocal_units = True)


def capillary_matched_spot(dne, den_unit: str = "centi", spot_unit: str = "micro"):
    """
    Calculate the laser spot size to match into a parabolic radial plasma channel with given channel depth.

    Parameters
    ----------
    dne : float
        Channel depth of parabolic profile to guide laser. Default unit: cm^-3
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    spot_unit : str, optional
        Order of magnitude of laser spot unit, by default "micro"

    Returns
    -------
    w0 : float
        Laser spot radius to 1/e^2 of peak intensity for guiding through the channel. Default unit: microns
    """

    dne_SI = magnitude_conversion_vol(dne, den_unit, "", reciprocal_units = True)
    w0     = 1/np.sqrt(pi*re*dne_SI)

    return magnitude_conversion(w0, "", spot_unit)
