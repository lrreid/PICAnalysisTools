"""
Test that ParticleEnergy and ParticleTransverseProperties classes functions as intended.
To use, you must set your own FolderPath, Simulation and Species to simulaton data which includes a particle species.
"""

import matplotlib.pyplot as plt
from scipy.constants import c
from openpmd_viewer import OpenPMDTimeSeries        # This should be used for pmd viewer version 1.x.
from PICAnalysisTools.Particle_Properties import ParticleEnergy, ParticleTransverseProperties, get_normalised_momentum
from PICAnalysisTools.utils.sim_path import set_sim_path
from PICAnalysisTools.utils.plot_limits import plt_limits_log

# from matplotlib import use
# use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

# FolderPath = r'/home/lewis'
FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
Species    = "electrons"

FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts   = OpenPMDTimeSeries(FilePath)
K    = 2                                       # snapshot to analyse
ctau = round(ts.t[K]*c*1e3,2)

bunch_thresh_MeV  = 600                                                     # Threshold for including electrons in particle diagnostic (MeV)
bunch_thresh_norm = get_normalised_momentum(bunch_thresh_MeV, "mega")       # Threshold for including electrons in particle diagnostic (normalised units)

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species=Species,
                                            iteration=ts.iterations[K], plot=False, select = {'uz':[bunch_thresh_norm, None] } )

PE = ParticleEnergy(ux, uy, uz, q, w, energy_unit = "Mega", charge_unit = "pico") # Create instance of Particle energy class

#%% Test beam energy properties

print("\n")
print("Mean energy: %0.2f MeV" % round(PE.Ek_mean,2) )
print("RMS energy spread: %0.2f MeV" % round(PE.de_rms,2) )
print("Percentage energy spread: %0.2f %%" % round(PE.PC_spread,2) )
print("Total charge: %0.2f pC" % round(PE.charge,2) )

#%% Plot beam energy spectrum

e_spec, [Spec_Min, Spec_Max] = PE.get_energy_spectrum(Spec_Res =  0.2, E_Round =  10) # Units match energy_unit of the instance of the Particle energy class
plt_y_min, plt_y_max         = plt_limits_log(e_spec[1], max_offset = -1)

fig, ax = plt.subplots()
plt.plot(e_spec[0], e_spec[1], '-', color='firebrick', linewidth=1.0)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.axis([Spec_Min, Spec_Max, plt_y_min, plt_y_max ])
ax.set_title(('c$\\tau$ = %0.2f mm' % ctau ), fontsize=fsize)
ax.set_xlabel('Beam energy (MeV)', fontsize=fsize)
ax.set_ylabel('Number of electrons', fontsize=fsize)
plt.grid(True)
plt.show()

#%% Create class for transverse particle properties in x and y

TX = ParticleTransverseProperties(x, ux, z, uz, w, r_unit = "micro", div_unit = "milli", emit_unit = "micro")
TY = ParticleTransverseProperties(y, uy, z, uz, w, r_unit = "micro", div_unit = "milli", emit_unit = "micro")

print("\n")
print("Beam properties in x")
print("Beam size: %0.2f um" % round(TX.beam_r,2) )
print("Divergence: %0.2f mrad" % round(TX.div_r,2) )
print("Normalised emittance: %0.2f mm mrad" % round(TX.eta_tr_norm_r,2) )
print("Twiss alpha: %0.2f" % round(TX.Twiss_alpha,2) )
print("Twiss beta: %0.2f m" % round(TX.Twiss_beta,2) )
print("Twiss gamma: %0.2f m^-1" % round(TX.Twiss_gamma,2) )

print("\n")
print("Beam properties in y")
print("Beam size: %0.2f um" % round(TY.beam_r,2) )
print("Divergence: %0.2f mrad" % round(TY.div_r,2) )
print("Normalised emittance: %0.2f mm mrad" % round(TY.eta_tr_norm_r,2) )
print("Twiss alpha: %0.2f" % round(TY.Twiss_alpha,2) )
print("Twiss beta: %0.2f m" % round(TY.Twiss_beta,2) )
print("Twiss gamma: %0.2f m^-1" % round(TY.Twiss_gamma,2) )

print("\n")
print("Beam length: %0.2f um" % round(TX.beam_z,2) )