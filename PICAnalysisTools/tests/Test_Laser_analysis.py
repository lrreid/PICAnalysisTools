"""
Test that LaserProperties class functions as intended.

To use, you must set your own FolderPath and Simulation to simulaton data which includes a laser pulse.
"""

import numpy as np
from openpmd_viewer.addons import LpaDiagnostics
from PICAnalysisTools.Laser_Properties import LaserProperties
from PICAnalysisTools.utils.sim_path import set_sim_path
from PICAnalysisTools.Laser_Properties import get_central_wavelength
fsize = 12

#%% Read data from files

FolderPath = r'/home/lewis'
# FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools'
Simulation = 'example_data'
FilePath, SimPath   = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = LpaDiagnostics(FilePath, check_all_files=False, backend='h5py')
LP = LaserProperties(ts, Coord = "x", Method = "fit", Profile_Method = "projection", centre_unit = "micro", spot_unit = "micro", time_unit = "femto")


#%% Single instance

snapshot = 0

centroid, a0, waist, tau = LP.get_laser_properties_single(Snapshot=snapshot)

print("\n")
print("For snapshot: %d" % snapshot)
print("Centroid position: %0.2f um" % np.round(centroid,2))
print("Normalised intensity: %0.2f" % np.round(a0,2))
print("Laser waist: %0.2f um" % np.round(waist,2))
print("Pulse duration: %0.2f fs" % np.round(tau,2) )
print("\n")

#%% Loop through all snapshots

output = LP.get_laser_properties_loop(save_txt=True, Sim_Path=SimPath)

#%% Find the central wavelength of the spectrum

avg, width, peak = get_central_wavelength(ts, snapshot)

print("\n")
print("Average wavelength: %0.2f nm" % np.round(avg,2))
print("RMS width: %0.2f nm" % np.round(width,2))
print("Wavelength at peak intensity: %0.2f nm" % np.round(peak,2))
print("\n")