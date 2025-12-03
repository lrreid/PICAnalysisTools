"""
Functions for calculating quantities and matching conditions for beam driven wakefield acceleration.
"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol

def resonant_ne(condition: str, sigmaZ, sigma_unit: str = "micro", den_unit: str = "centi"):
    """
    Calculate the resonant plasma density for beam driven wakefield acceleration

    Parameters
    ----------
    condition : str
        Condition for choosing plasma density. Choose between
        k_p * sigma_z = sqrt(2)
        k_p * sigma_z = pi/4
    sigmaZ : float
        RMS longitudial beam length. Default unit: microns
    sigma_unit : str, optional
        Order of magnitude of beam length unit, by default "micro"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Resonant plasma density for given beam length. Default unit: cm^-3
    """
    
    sigmaZ_SI = magnitude_conversion(sigmaZ, sigma_unit, "" )

    if condition == "sqrt(2)":
        n_e = 2*epsilon_0*m_e*(c/(sigmaZ_SI*e))**2
    elif condition == "pi/4":
        n_e = epsilon_0*m_e*((c*pi)/(4*sigmaZ_SI*e))**2    
    else:
        n_e = 0
        print("Resonance condition incorrectly defined!")
    
    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)

def sigma_matched(Ek, n_e, eta_N, energy_unit: str = "Mega", den_unit: str = "centi", emit_unit: str ="micro", sigma_unit: str ="micro"):
    """
    Calculate matched transverse beam size for given plasma density.

    Parameters
    ----------
    Ek : float
        Particle species central energy. Default unit: MeV
    n_e : float
        Plasma denstiy. Default unit: cm^-3
    eta_N : float
        Particle species normalised emittance. Default unit:mm mrad.
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "Mega"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    emit_unit : str, optional
        Order of magnitude of normalised emittance unit, by default "micro"
    sigma_unit : str, optional
        Order of magnitude of transverse beam size unit, by default "micro"

    Returns
    -------
    sigma_m: float
        Matched transverse beam size (RMS) to plasma. Default unit: microns
    """

    Ek_SI    = magnitude_conversion(Ek, energy_unit, "" )
    n_e_SI   = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    eta_N_SI = magnitude_conversion(eta_N, emit_unit, "")
    
    gamma0  = ((Ek_SI*e)/(m_e*c**2)) + 1
    sigma_m = ((2*c**2*epsilon_0*(eta_N_SI*1e-6)**2*m_e)/(gamma0*n_e_SI*e**2))**0.25
    
    return magnitude_conversion(sigma_m, "", sigma_unit)


def charge_density(sig_r, sig_z, Q, sigma_unit: str ="micro", charge_unit: str = "pico", den_unit: str = "centi"):
    """
    Calculate the charge density of a particle beam with Gaussian distribution.

    Parameters
    ----------
    sig_r : float
        Transverse beam size (RMS). Default unit: microns
    sig_z : float
        Longitudinal beam size (RMS). Default unit: microns
    Q : float
        Total particle beam charge. Default unit: pC
    sigma_unit : str, optional
        Order of magnitude of transverse beam size unit, by default "micro"
    charge_unit : str, optional
         Order of magnitude of beam charge unit, by default "pico"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    nb: float
        Charge density of particle beam. Default unit: cm^-3
    """

    sig_r_SI = magnitude_conversion(sig_r, sigma_unit, "")
    sig_z_SI = magnitude_conversion(sig_z, sigma_unit, "")
    Q_SI     = magnitude_conversion(Q, charge_unit, "")
        
    nb = (Q_SI) / (e * ((2*pi)**(3/2)) * (sig_r_SI)**2  * (sig_z_SI) ) # Beam charge density (m-3)
    
    return magnitude_conversion_vol(nb, "", den_unit, reciprocal_units = True)

def sigr_for_nonlinear_thres(sig_z, Q, condition, sigma_unit: str ="micro", charge_unit: str = "pico"):
    """
    Calculate longitudinal beam size for beam density to match plasma density (nb/ne = 1)
    i.e. Threshold for non-linear regime of pwfa for given
    bunch length and charge.

    Parameters
    ----------
    sig_z : float
        Longitudinal beam size (RMS). Default unit: microns
    Q : float
        Total particle beam charge. Default unit: pC
    condition : str
        Condition for choosing plasma density. Choose between
        k_p * sigma_z = sqrt(2)
        k_p * sigma_z = pi/4
    sigma_unit : str, optional
        Order of magnitude of beam length unit, by default "micro"
    charge_unit : str, optional
         Order of magnitude of beam charge unit, by default "pico"

    Returns
    -------
    sig_r: float
        Transverse beam size which gives nb = ne. Default unit: um
    """

    sig_z_SI = magnitude_conversion(sig_z, sigma_unit, "")
    Q_SI     = magnitude_conversion(Q, charge_unit, "")
    
    if condition == "sqrt(2)":
        sig_r = (1/(2*2**(0.25)*(pi**(3/4))*c))*np.sqrt((Q_SI*sig_z_SI*e)/(epsilon_0*m_e))
    elif condition == "pi/4":
        sig_r = ((2*2**(0.25))/(c*pi**(7/4))) * np.sqrt( (Q_SI*e*sig_z_SI)/(epsilon_0*m_e) )
    else:
        sig_r = 0
    
    return magnitude_conversion(sig_r, "", sigma_unit)



def sigz_for_nonlinear_thres(sig_r, Q, condition, sigma_unit: str ="micro", charge_unit: str = "pico"):
    """
    Calculate longitudinal beam size for beam density to match plasma density (nb/ne = 1)
    i.e. Threshold for non-linear regime of pwfa for given
    bunch length and charge.

    Parameters
    ----------
    sig_r : float
        Transverse beam size (RMS). Default unit: microns
    Q : float
        Total particle beam charge. Default unit: pC
    condition : str
        Condition for choosing plasma density. Choose between
        k_p * sigma_z = sqrt(2)
        k_p * sigma_z = pi/4
    sigma_unit : str, optional
        Order of magnitude of beam length unit, by default "micro"
    charge_unit : str, optional
         Order of magnitude of beam charge unit, by default "pico"

    Returns
    -------
    sig_z: float
        Longitudinal beam size which gives nb = ne. Default unit: um
    """

    sig_r_SI = magnitude_conversion(sig_r, sigma_unit, "")
    Q_SI     = magnitude_conversion(Q, charge_unit, "")
    
    if condition == "sqrt(2)":
        sig_z = (4*np.sqrt(2)*(pi**(3/2))*sig_r_SI**2*epsilon_0*m_e*c**2)/(e*Q_SI)
    elif condition == "pi/4":
        sig_z = (((pi**(7/2))*(sig_r_SI)**2*epsilon_0*m_e*c**2)/(4*np.sqrt(2)*e*Q_SI))
    else:
        sig_z = 0
        
    return magnitude_conversion(sig_z, "", sigma_unit)