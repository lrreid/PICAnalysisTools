
"""
Check that all unit order of magnitde conversions are working correctly by testing them all in both directions.

"""

import numpy as np
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_area, magnitude_conversion_vol


#%% Length check

k = 0

print("########### 1D conversion tests ###########")

lambda_p    = 106e-6
lambda_p_um = magnitude_conversion(lambda_p, "", "micro")
print("plasma wavelength in SI units: %0.3e m" % (lambda_p) )
print("plasma wavelength: %d um" % (lambda_p_um) )
print("Expected answer: %d m" % (lambda_p*1e6) )
if np.round(lambda_p_um/(lambda_p*1e6),6) == 1:
    print("1D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1

print("\n")
lambda_0_nm    = 800
lambda_0       = magnitude_conversion(lambda_0_nm, "nano", "")
print("laser wavelength: %d nm" % (lambda_0_nm) )
print("laser wavelength in SI units: %0.2e m" % (lambda_0) )
print("Expected answer: %0.2e m" % (lambda_0_nm*1e-9) )
if np.round(lambda_0/(lambda_0_nm*1e-9),6) == 1:
    print("1D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1
    
print("\n")
print("\n")

#%% Area check
print("########### Area conversion tests: ###########")

Intensity = 1.15e19                             # peak laser intensity (Wcm^-2)
Int_SI    = magnitude_conversion_area(Intensity, "centi", "")
print("Peak intensity: %0.2e Wcm^-2" % Intensity)
print("Peak intensity in SI units: %0.2e Wm^-2" % Int_SI)
print("Expected answer: %0.2e m" % (Intensity*1e4) )
if np.round(Int_SI/(Intensity*1e4),6) == 1:
    print("2D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1

print("\n")
    
Intensity_SI = 5.67e22                          # peak laser intensity (Wcm^-2)
Int    = magnitude_conversion_area(Intensity_SI, "", "centi")
print("Peak intensity in SI units: %0.2e Wm^-2" % Intensity_SI)
print("Peak intensity in converted units: %0.2e Wcm^-2" % Int)
print("Expected answer: %0.2e m" % (Intensity_SI*1e-4) )
if np.round(Int/(Intensity_SI*1e-4),6) == 1:
    print("2D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1


print("\n")
print("\n")


#%% Volume check
print("########### Volume conversion tests: ###########")

n_e    = 1e17                       # plasma density (cm^-3)
n_e_SI = magnitude_conversion_vol(n_e, "centi", "")
print("Plasma density: %0.2e cm^-3" % n_e)
print("Plasma density in SI units: %0.2e m^-3" % n_e_SI)
print("Expected answer: %0.2e m^-3" % (n_e*1e6) )
if np.round(n_e_SI/(n_e*1e6),6) == 1:
    print("3D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1
    
print("\n")

N_E_SI = 1e25                       # plasma density (cm^-3)
N_E    = magnitude_conversion_vol(N_E_SI, "", "centi")
print("Plasma density in SI units: %0.2e m^-3" % N_E_SI)
print("Plasma density in cm^-3: %0.2e cm^-3" % N_E)
print("Expected answer: %0.2e m^-3" % (N_E_SI*1e-6) )
if np.round(N_E/(N_E_SI*1e-6),6) == 1:
    print("3D unit conversion test passed!")
else:
    print("test failed")
    k = k + 1

print("\n")


if k > 0:
    print("\n######################### Warning #########################\nOne or more tests failed\n\n")
else:
    print("\n######################### All tests passed #########################")