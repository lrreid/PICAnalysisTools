"""
Functions for calculating quantities and matching conditions for laser driven wakefield acceleration.


All of these need re-written to add options for units and to be compatible with other functions within package.

"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol
re = (e**2)/(4*pi*epsilon_0*m_e*c**2)                   # classial electron radius


def resonant_density(tau_FWHM, time_unit : str = "femto", den_unit : str = "centi"):

    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "" )

    n_e_sqrt2 = (16*np.log(2)*epsilon_0*m_e)/(e**2*tau_FWHM_SI**2)   # Plasma density (electrons.meters^-3). Density according to k_p*L_RMS = sqrt(2). Linear wakefield.
    n_e_1     = (8*np.log(2)*epsilon_0*m_e)/(e**2*tau_FWHM_SI**2)    # Plasma density (electrons.meters^-3). Density according to k_p*L_RMS = 1.       Non-linear wakefield.


    return magnitude_conversion_vol(n_e_sqrt2, "", den_unit, reciprocal_units = True), magnitude_conversion_vol(n_e_1, "", den_unit, reciprocal_units = True)


def pump_depletion_length():

    return None

def dephasing_length():

    return None

def self_injection_Threshold():
    # Add Thomas and Benedetti models.

    return None


def density_for_matched_spot(): # ω0 = λp/2=kp/π

    return None


def matched_spot_from_density(): # ω0 = λp/2=kp/π

    return None


def max_energy_scaling(E, w0, lambda0, tau, ne):

    # This needs re-written to add options for units and to be compatible with other functions within package

    a0       = ((e*lambda0)/(pi*m_e*c**2*w0)) * np.sqrt( E/(pi*epsilon_0*c*tau) )       # Laser amplitude
    omega_0  = (2*pi*c)/lambda0                                                         # Angular frequency of laser (rad/s)
    w_p      = np.sqrt( (ne*e**2) / (epsilon_0*m_e) )                                   # Plasma frequency
    
    E_max_lin    = a0**2*(omega_0/w_p)**2 * ((m_e*c**2)/e)              # All equations from Lu paper table 1, final column 
    E_max_1D     = 4*a0**2*(omega_0/w_p)**2 * ((m_e*c**2)/e)
    E_max_3D     = (2/3) * a0 * ((omega_0/w_p)**2) * ((m_e*c**2)/e)
    E_max_Pukhov = ((omega_0**2*tau*a0**(3/2))/w_p) * ((m_e*c**2)/e)        # As in Lu paper (final row in table)
    
    return E_max_lin, E_max_1D, E_max_3D, E_max_Pukhov


def max_energy_scaling_engineering(E, w0, lambda0, tau, ne):

    # This needs re-written to add options for units and to be compatible with other functions within package

    omega_0   = (2*pi*c)/lambda0
    P_GW      = (E/tau)*1e-9 
    P_TW      = (E/tau)*1e-12
    ne_cm     = ne*1e-6
    nc        = (epsilon_0*m_e*omega_0**2)/(e**2)                                                       # Critical density (m^-3)
    lambda_um = lambda0*1e6

    PC_const = (1e-9*8*pi*epsilon_0*m_e**2*c**5)/(e**2)                                                 # (half) critical power constant (GW)
    P_rel    = 0.5*PC_const                                                                             # Relativistic power unit (GW)

    E_Lu     = 1e3 * 1.7 * ((P_TW/100)**(1/3)) * ((1e18/ne_cm)**(1/3)) * ((0.8/lambda_um)**(4/3))       # Max energy according to Lu (MeV) - Lu eqn 6
    E_Pukhov = 1e-6 *0.65 * ((m_e*c**2)/e) * np.sqrt(P_GW/P_rel) * (c*tau/lambda0)                      # Max energy according to Pukhov (MeV) - Lu eqn 12
    E_PAlt   = 1e-6 * 0.16 * ((m_e*c**2)/e) * ((c*tau)/w0) * ((P_GW/P_rel)**(2/3)) * ((nc/ne)**(1/3))   # Max energy according to Pukhov (MeV) - re-rwitten by Lu (eqn 13)
    
    return E_Lu, E_Pukhov, E_PAlt


def betatron_wavelength():

    return None

def betatron_critical_energy():

    return None


def capillary_channel_depth(w0, spot_unit : str = "micro", den_unit : str = "centi"):


    w0_SI = magnitude_conversion(w0, spot_unit, "" )
    dne   = 1/(pi*re*w0_SI**2)                  # depth of plasma channel (m^-3)

    return magnitude_conversion_vol(dne, "", den_unit, reciprocal_units = True)


def capillary_matched_spot(dne, den_unit : str = "centi", spot_unit : str = "micro"):

    dne_SI = magnitude_conversion_vol(dne, den_unit, "", reciprocal_units = True)
    w0     = 1/np.sqrt(pi*re*dne_SI)

    return magnitude_conversion(w0, "", spot_unit)
