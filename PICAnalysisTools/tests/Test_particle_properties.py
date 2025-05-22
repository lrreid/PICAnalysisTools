"""
Test that ParticleEnergy and ParticleTransverseProperties classes functions as intended.

To use, you must set your own FolderPath, Simulation and Species to simulaton data which includes a particle species.

TO DO:
    - Add plotting for phase space. Make sure these functions work. 

"""

import matplotlib.pyplot as plt
from scipy.constants import c
from openpmd_viewer import OpenPMDTimeSeries        # This should be used for pmd viewer version 1.x.
from PICAnalysisTools.Particle_Properties import ParticleEnergy, ParticleTransverseProperties
from PICAnalysisTools.utils.sim_path import set_sim_path
# from PICAnalysisTools.utils.rounding import roundup
from PICAnalysisTools.utils.plot_limits import plt_limits_log

# from matplotlib import use
# use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

FolderPath = r'D:\Lewis\fbpic\20250501_pwfa_at_CLARA'
Simulation = '20250501_sigz_2p2um_sigr_5um_kpLrms_pi4'
Species = "beam"

FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts   = OpenPMDTimeSeries(FilePath)
K    = 12                                       # snapshot to analyse
ctau = round(ts.t[K]*c*1e3,2)

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species=Species, iteration=ts.iterations[K], plot=False)

PE = ParticleEnergy(ux, uy, uz, q, w, energy_unit = "Mega", charge_unit = "pico") # Create instance of class

#%% Test beam energy properties

print("\n")
print("Mean energy: %0.2f MeV" % round(PE.Ek_mean,2) )
print("RMS energy spread: %0.2f MeV" % round(PE.de_rms,2) )
print("Percentage energy spread: %0.2f %%" % round(PE.PC_spread,2) )
print("Total charge: %0.2f pC" % round(PE.charge,2) )

#%% Plot beam energy spectrum

e_spec, Spec_Min, Spec_Max, E_Bins = PE.get_energy_spectrum(Spec_Res =  0.2, E_Round =  10)
plt_y_min, plt_y_max               = plt_limits_log(e_spec[0], max_offset = -1)

fig, ax = plt.subplots()
plt.plot(E_Bins[:-1], e_spec[0], '-', color='firebrick', linewidth=1.0)
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

div_x, eta_tr_norm_x, beam_x, Twiss_alpha_x, Twiss_beta_x, Twiss_gamma_x = TX.transverse_beam_properties()
div_y, eta_tr_norm_y, beam_y, Twiss_alpha_y, Twiss_beta_y, Twiss_gamma_y = TY.transverse_beam_properties()
beam_z = TX.beam_length()


print("\n")
print("Beam properties in x")
print("Divergence: %0.2f mrad" % round(div_x,2) )
print("Normalised emittance: %0.2f mm mrad" % round(eta_tr_norm_x,2) )
print("Beam size: %0.2f um" % round(beam_x,2) )
print("Twiss alpha: %0.2f" % round(Twiss_alpha_x,2) )
print("Twiss beta: %0.2f m" % round(Twiss_beta_x,2) )
print("Twiss gamma: %0.2f m^-1" % round(Twiss_gamma_x,2) )

print("\n")
print("Beam properties in y")
print("Divergence: %0.2f mrad" % round(div_y,2) )
print("Normalised emittance: %0.2f mm mrad" % round(eta_tr_norm_y,2) )
print("Beam size: %0.2f um" % round(beam_y,2) )
print("Twiss alpha: %0.2f" % round(Twiss_alpha_y,2) )
print("Twiss beta: %0.2f m" % round(Twiss_beta_y,2) )
print("Twiss gamma: %0.2f m^-1" % round(Twiss_gamma_y,2) )

print("\n")
print("Beam length: %0.2f um" % round(beam_z,2) )