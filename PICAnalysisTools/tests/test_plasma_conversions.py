"""

Test PlasmaDen_Conversions class.
Check that the functions give the expected answers with the correct units. Compare against known good solutions.

"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.plasma_calcs import PlasmaDen_Conversions
from PICAnalysisTools.utils.plasma_calcs import *

#%% Script variables

n_e          = 1e17                             # plasma density (cm^-3)
lambda_laser = 800                              # laser central wavelength (nm)
freq_laser   = (2*pi*c)/(lambda_laser*1e-9)     # laser angular frequency (rad/s)

#%% Run instance of plasma density calculations class

plasma_info = PlasmaDen_Conversions(n_e, den_unit="centi", wavelength_unit="micro", freq_unit="Tera", skin_unit = "micro", time_unit = "femto")

#%% Manually preform calculations

n_e_SI = n_e * 1e6        # plasma density in SI units (m^-3)

w_p = np.sqrt((n_e_SI*e**2)/(epsilon_0*m_e))                            # Plasma frequency (rad/s)
f_p = 1e-12 * (1/(2*pi)) * np.sqrt( (n_e_SI*e**2) / (epsilon_0*m_e) )   # Plasma frequency (THz)
l_p = 1e6*2*pi*c*np.sqrt( (epsilon_0*m_e)/(e**2*n_e_SI) )               # Plasma wavelength (um)
k_p = (e/c) * np.sqrt(n_e_SI/(epsilon_0*m_e))                           # Plasma wavevector
skin_depth = (c/e) * np.sqrt((epsilon_0*m_e)/n_e_SI) * 1e6              # plasma skin depth (um)
T_p = 1e15 * ((2*pi)*np.sqrt((epsilon_0*m_e)/(n_e_SI*e**2)) )           # Plasma period (fs)

print("From direct calculation:")
print("For a plasma density of %0.2e cm^-3" % n_e)
print("Plasma frequency:\t\t%0.3e rad/s" % w_p )
print("Plasma wavelength:\t\t%0.3f um" % np.round(l_p,2) )
print("Plasma wavevector:\t\t%0.3e rad/s" % k_p )
print("Plasma period:\t\t\t%0.3f fs" % np.round(T_p,2) )
print("Plasma frequency:\t\t%0.3f THz" % np.round(f_p,3) )
print("Plasma skin depth:\t\t%0.2f um" % np.round(skin_depth,2) )

print("\n\n")

print("From PlasmaDen_Conversions class:")
print("For a plasma density of %0.2e cm^-3" % n_e)
print("Plasma frequency:\t\t%0.3e rad/s" % plasma_info.w_p )
print("Plasma wavelength:\t\t%0.3f um" % np.round(plasma_info.lambda_p,2) )
print("Plasma wavevector:\t\t%0.3e rad/s" % plasma_info.k_p )
print("Plasma period:\t\t\t%0.3f fs" % np.round(plasma_info.T_p,2) )
print("Plasma frequency:\t\t%0.3f THz" % np.round(plasma_info.f_p,3) )
print("Plasma skin depth:\t\t%0.2f um" % np.round(plasma_info.skin_depth,2) )



print("\n\n")

if np.round(plasma_info.w_p, 6) == np.round(w_p, 6):
    print("plasma angular frequency correct")
else:
    print("################ plasma angular frequency calculations disagree ################")
    
if np.round(plasma_info.lambda_p, 6) == np.round(l_p, 6):
    print("plasma wavelength correct")
else:
    print("################ plasma wavelength calculations disagree ################")
    
if np.round(plasma_info.k_p, 6) == np.round(k_p, 6):
    print("plasma wavevector correct")  
else:
    print("################ plasma wavevector calculations disagree ################")
    
if np.round(plasma_info.T_p, 6) == np.round(T_p, 6):
    print("plasma period correct")
else:
    print("################ plasma period calculations disagree ################")
    
if np.round(plasma_info.f_p, 6) == np.round(f_p, 6):
    print("plasma frequency correct")
else:
    print("################ plasma frequency calculations disagree ################")
    
if np.round(plasma_info.skin_depth, 6) == np.round(skin_depth, 6):
    print("plasma skin depth correct")
else:
    print("################ plasma skin depth calculations disagree ################")

print("\n\n")    

#%% Test other plasma functions

print("Each of these tests should return the original plasma density")

ne_freq = plasma_density_from_frequency(plasma_info.w_p, den_unit="centi")
ne_lambda = plasma_density_from_wavelength(plasma_info.lambda_p, wavelength_unit = "micro", den_unit="centi")
ne_period = plasma_density_from_period(plasma_info.T_p, time_unit = "femto", den_unit="centi")
ne_kp     = plasma_density_from_wavevector(plasma_info.k_p, den_unit="centi")

print("expected answer:\t\t\t%0.2e cm^-3" % n_e )
print("Plasma density from angular frequency:\t%0.2e cm^-3" % ne_freq )
print("Plasma density from wavelength:\t\t%0.2e cm^-3" % ne_lambda )
print("Plasma density from period:\t\t%0.2e cm^-3" % ne_period )
print("Plasma density from wavevector:\t\t%0.2e cm^-3" % ne_kp )

print("\n\n")  
#%% Test critical power for self-focusing functions

print("Testing functions for critical power for self-focusing")

P_Crit = critical_power(n_e, lambda0 = lambda_laser, den_unit = "centi", wavelength_unit = "nano", power_unit="Tera")

w0    = (2*pi*c)/(lambda_laser*1e-9)                                # Laser angular frequency
const = 8*pi*epsilon_0*((m_e**2*c**5)/(e**2))       # Constant for critical power (W)
Pcrit = const*((w0/w_p)**2)*1e-12                   # Critical power (TW)

print("Critical power:\t\t\t%0.3f TW" % np.round(Pcrit,3) )
print("Critical power from function:\t%0.3f TW" % np.round(P_Crit,3) )
if np.round(P_Crit,6) == np.round(Pcrit,6):
    print("Critical power calculations agree")
else:
    print("################ Critical power calculations disagree ################")

ne_critical_power = plasma_density_from_critical_power(Pcrit, lambda0=lambda_laser, power_unit="Tera", wavelength_unit = "nano", den_unit = "centi" )

print("\n\n")  

print("Reverse critical power calculation")
print("Plasma density from calculation:\t\t%0.2e cm^-3" % n_e )
print("Plasma density from reverse calculation:\t%0.2e cm^-3" % ne_critical_power )

if np.round(n_e/ne_critical_power,6) == 1:
    print("Reverse critical power calculations agree")
else:
    print("################ Reverse critical power calculations disagree ################")

print("\n\n") 


#%% Test critical density functions

ncrit      = (1e-6*epsilon_0*m_e*freq_laser**2)/(e**2)      # plasma critical density (cm^-3) 

n_c_lambda = critical_density(lambda0 = lambda_laser, wavelength_unit = "nano", den_unit = "centi" )
n_c_freq   = critical_density_from_frequency(omega_l=freq_laser, den_unit = "centi")


print("Critical plasma density calculation")
print("Critical plasma density:\t\t\t%0.2e cm^-3" % ncrit )
print("Critical plasma density from wavelength:\t%0.2e cm^-3" % n_c_lambda )
print("Critical plasma density from frequency:\t\t%0.2e cm^-3" % n_c_freq )
if np.round(ncrit/n_c_lambda,6) == 1:
    print("Critical plasma density from wavelength calculations agree")
else:
    print("################ Critical plasma density from wavelength calculations disagree ################")

if np.round(ncrit/n_c_freq,6) == 1:
    print("Critical plasma density from frequency calculations agree")
else:
    print("################ Critical plasma density from frequency calculations disagree ################")

print("\n\n") 

#%% Test group velocity function

v_group = laser_group_velocity(lambda0 = lambda_laser, n_e = n_e, wavelength_unit = "nano", den_unit = "centi")
u_group = c * np.sqrt(1-(((lambda_laser*1e-9)/(l_p*1e-6))**2))              # laser group velocity (m/s)

print("laser group velocity calculation")
print("laser group velocity:\t\t\t%0.4e m/s" % u_group )
print("laser group velocity from function:\t%0.4e m/s" % v_group )
if np.round(u_group/v_group,6) == 1:
    print("laser group velocity calculations agree")
else:
    print("################ laser group velocity calculations disagree ################")
print("\n\n") 


#%% Plasma wavebreaking field

E_WB = wavebreaking_field(n_e, den_unit = "centi", field_unit = "Giga")
E_wb = ((c*m_e*w_p)/e)*1e-9                                                 # Wave breaking field (GV/m)

print("Plasma wavebreaking field calculation")
print("Plasma wavebreaking field:\t\t\t%0.2f GV/m" % np.round(E_wb,3) )
print("Plasma wavebreaking field from function:\t%0.2f GV/m" % np.round(E_WB,3) )
if np.round(E_wb/E_WB,6) == 1:
    print("Plasma wavebreaking field calculations agree")
else:
    print("################ Plasma wavebreaking field calculations disagree ################")
print("\n\n")

