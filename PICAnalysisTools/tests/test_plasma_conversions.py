"""

Test PlasmaDen_Conversions class.
Check that the functions give the expected answers with the correct units. Compare against known good solutions.

"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.plasma_calcs import PlasmaDen_Conversions


#%% Run function

n_e = 1e18              # plasma density (cm^-3)

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
    