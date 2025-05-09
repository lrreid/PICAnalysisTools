"""

Simple plasma calculations and conversions

TO DO:
    - All functions need tested!
    - Add doc strings
    - Add lwfa scalings? separate script?

    
Calculations for plasma wavelength, skin depth, frequency and angular frequency are all wrong. 
Numbers are correct but order of magnitude is not. Think this is an error in the unit conversion function!

"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol

class PlasmaDen_Conversions():

    def __init__(self, n_e, den_unit="centi", wavelength_unit="micro", freq_unit="Tera", skin_unit = "micro", time_unit = "femto"):

        self.n_e             = n_e
        self.den_unit        = den_unit
        self.wavelength_unit = wavelength_unit
        self.freq_unit       = freq_unit
        self.skin_unit       = skin_unit
        self.time_unit       = time_unit
        self.n_e_SI          = magnitude_conversion_vol(self.n_e, self.den_unit, "", reciprocal_units = True)        # convert plasma density to si units (m^-3)
        self.w_p, self.f_p   = self.plasma_frequency_from_density()                         # Plasma (angular) frequency
        self.lambda_p        = self.plasma_wavelength_from_density()                        # Plasma wavelength
        self.k_p, self.skin_depth = self.plasma_wavevector_from_density()                   # plasma wavevector and skin depth
        self.T_p             = self.plasma_period_from_density()                            # Plasma period

    def plasma_wavelength_from_density(self):

        lambda_p  = (2*pi*c)*np.sqrt( (epsilon_0*m_e)/(e**2*self.n_e_SI) )
        
        return magnitude_conversion(lambda_p, "", self.wavelength_unit)


    def plasma_frequency_from_density(self):

        w_p = np.sqrt( (self.n_e_SI*e**2) / (epsilon_0*m_e) )                # Plasma (angular) frequency (rad/s)
        f_p = (1/(2*pi)) * np.sqrt( (self.n_e_SI*e**2) / (epsilon_0*m_e) )   # Plasma ("normal"?) frequency (Hz)

        return w_p, magnitude_conversion(f_p, "", self.freq_unit)


    def plasma_period_from_density(self):

        T_p = ((2*pi)*np.sqrt((epsilon_0*m_e)/(self.n_e_SI*e**2)) )          # Plasma period (s)

        return magnitude_conversion(T_p, "", self.time_unit)


    def plasma_wavevector_from_density(self):

        k_p = (e/c) * np.sqrt(self.n_e_SI/(epsilon_0*m_e))
        skin_depth = 1/k_p

        return k_p, magnitude_conversion(skin_depth, "", self.skin_unit)

#### end of class ####

def plasma_density_from_frequency(w_p, den_unit="centi"):

    n_e = (w_p**2*epsilon_0*m_e)/(e**2)

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def plasma_density_from_wavelength(lambda_p, wavelength_unit = "", den_unit="centi"):

    lambda_p_SI = magnitude_conversion(lambda_p, wavelength_unit, "" ) # convert plasma period to time in SI units (s).
    
    n_e = epsilon_0*m_e*((2*pi*c)/(lambda_p_SI*e))**2

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def plasma_density_from_period(T_p, time_unit = "femto", den_unit="centi"):

    T_p_SI = magnitude_conversion(T_p, time_unit, "" ) # convert plasma period to time in SI units (s).

    n_e = epsilon_0*m_e* ((2*pi)/(e*T_p_SI))**2

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)

def plasma_density_from_wavevector(k_p, den_unit="centi"):

    n_e = ((c*k_p)/e)**2 * epsilon_0*m_e

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def critical_power(n_e, lambda0 = 800e-9, den_unit = "centi", wavelength_unit = "nano", power_unit="Tera"):
    # Calculate the critical power for laser self-focusing in a plasma.

    n_e_SI     = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                  # Get laser wavelength in SI unit (m)

    Pcrit = (32*pi**3*epsilon_0**2*m_e**3*c**7)/(n_e_SI*e**4*lambda0_SI**2)  # critical power in W

    return magnitude_conversion(Pcrit, "", power_unit)


def plasma_density_from_critical_power(Pcrit, lambda0, power_unit="Tera", wavelength_unit = "nano", den_unit = "centi" ):

    Pcrit_SI   = magnitude_conversion(Pcrit, power_unit, "")                         # Get critical power in SI unit (Watt)
    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                  # Get laser wavelength in SI unit (m)
    n_e        = (32*pi**3*epsilon_0**2*m_e**3*c**7)/(Pcrit_SI*e**4*lambda0_SI**2)      # Plasma density at critical power (m^-3)

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def critical_density(lambda0, wavelength_unit = "nano", den_unit = "centi" ):

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")         # Get laser wavelength in SI unit (m)
    n_c        = m_e*epsilon_0*((2*pi*c)/(lambda0_SI*e))**2                 # Critical density (m^-3)

    return magnitude_conversion_vol(n_c, "", den_unit, reciprocal_units = True)


def critical_density_from_frequency(omega_l, den_unit = "centi"):

    n_c = (m_e*epsilon_0*omega_l**2)/(e**2)                                 # Critical density (m^-3)

    return magnitude_conversion_vol(n_c, "", den_unit, reciprocal_units = True)


def plasma_refractive_index(lambda0, n_e, wavelength_unit = "nano", den_unit = "centi"):

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                  # Get laser wavelength in SI unit (m)
    n_e_SI     = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    eta        = np.sqrt(1 - ( (n_e_SI/(4*epsilon_0*m_e)) * ((e*lambda0_SI)/(pi*c))**2 ) )

    return eta

def laser_group_velocity(lambda0, n_e, wavelength_unit = "nano", den_unit = "centi"):
    
    v_group = c * plasma_refractive_index(lambda0, n_e, wavelength_unit, den_unit)

    return v_group


def wavebreaking_field(n_e, den_unit = "centi", field_unit = "Giga"):

    #w_p  = plasma_frequency_from_density(n_e)[0]
    #E_wb = (c*m_e*w_p)/e                           # Wave breaking field (V/m)

    n_e_SI = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    E_wb   = c * np.sqrt((m_e*n_e_SI)/epsilon_0)

    return magnitude_conversion(E_wb, "", field_unit) 
 