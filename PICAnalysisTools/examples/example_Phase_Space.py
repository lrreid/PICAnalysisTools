"""
Test class which gets the phase space of a particle beam

To use, you must set your own FolderPath, Simulation and Species to simulaton data which includes a particle species.

TO DO:
    - 

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from openpmd_viewer import OpenPMDTimeSeries        # This should be used for pmd viewer version 1.x.
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter

from PICAnalysisTools.utils.white_background_colormap import cmap_white
from PICAnalysisTools.Particle_Properties import PhaseSpace, BeamProjection, get_normalised_momentum
from PICAnalysisTools.utils.particle_selection import radial_selection
from PICAnalysisTools.utils.sim_path import set_sim_path
from PICAnalysisTools.utils.plot_limits import plt_limits_log
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion

# from matplotlib import use
# use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

# FolderPath = r'/home/lewis'
FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
Species    = "electrons"

FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')

K = 2               # snapshot to analyse

bunch_thresh_MeV  = 600                                                     # Threshold for including electrons in particle diagnostic (MeV)
bunch_thresh_norm = get_normalised_momentum(bunch_thresh_MeV, "mega")       # Threshold for including electrons in particle diagnostic (normalised units)

_, info_Ex = ts.get_field( iteration=ts.iterations[K], field = 'E', m=1, coord = 'x')

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species=Species,
                                                iteration=ts.iterations[K], plot=False, select = {'uz':[bunch_thresh_norm, None] } )

# Take radial selection from full dataset
x, y, z, ux, uy, uz, w = radial_selection(0, 6e-6, x, y, z, ux, uy, uz, w)

ctau = round(ts.t[K]*c*1e3,2)

#%% Energy Z space

Phase_Space = PhaseSpace(x, y, z, ux, uy, uz, w, r_unit = "micro", energy_unit= "mega", div_unit= "milli", time_unit= "femto")      # Create instance of transverse phase space class

E_Phase, plt_limits, E_Bins, R_Bins, Energy_line, z_line = Phase_Space.Energy_z_space(z_res=0.1, z_round=5, Spec_Res=0.2, E_Round=50, Centre_z=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(E_Phase, cmap=cmap_white(plt.cm.inferno), extent=(plt_limits), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons', size=fsize)
plt.clim(plt_limits_log(E_Phase, min_offset = 0, max_offset = 1))

plt.plot(Energy_line, E_Bins, color = 'blue', linewidth = 0.5)
plt.plot(R_Bins, z_line, color = 'blue', linewidth = 0.5)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('E (MeV)', fontsize=fsize)
plt.axis(plt_limits)
plt.show()

#%% Energy time space

E_Phase_t, Spec_Min, Spec_Max, E_Bins, t_min, t_max, T_Bins, Row_sum, T_line, Col_sum = Phase_Space.Energy_time_space(Spec_Res=0.2, E_Round=50, t_res=0.1, t_round=5, Find_Lineouts=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(E_Phase_t), cmap=cmap_white(plt.cm.inferno), extent=([t_min, t_max, Spec_Min, Spec_Max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(E_Phase_t, min_offset = 0, max_offset = 1))

plt.plot(Row_sum, E_Bins[:-1], color = 'forestgreen', linewidth = 1)
plt.plot(T_line[:-1], Col_sum, color = 'forestgreen', linewidth = 1)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('t (fs)', fontsize=fsize)
ax1.set_ylabel('E (MeV)', fontsize=fsize)
plt.axis([t_min, t_max, Spec_Min, Spec_Max])
plt.show()


#%% Test transverse phase space and create plot

Trans_phase = Phase_Space.Div_r_space('x', r_res=0.20, r_round=10, div_round=10, div_res=0.1,
                              Find_Lineouts=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(Trans_phase[0]), cmap=cmap_white(plt.cm.viridis), extent=([Trans_phase[3], Trans_phase[4], -Trans_phase[1], Trans_phase[1]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(Trans_phase[0], min_offset = 0, max_offset = 1))

plt.plot(Trans_phase[6], Trans_phase[2][:-1], color = 'red', linewidth = 1)
plt.plot(Trans_phase[7], Trans_phase[8], color = 'red', linewidth = 1)

ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x\' (mrad)', fontsize=fsize)
plt.axis([Trans_phase[3], Trans_phase[4], -Trans_phase[1], Trans_phase[1]])
plt.show()

#%% Test beam projection class - YAG screen projection

BP_xy = BeamProjection(x, y, w)

YAG, Hist_Limits, _, _ = BP_xy.beam_projection(r_res=0.2, r_round=1, independant_bins=False, equal_bins=True)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(YAG), cmap=cmap_white(plt.cm.winter), extent=([Hist_Limits[0], Hist_Limits[1], Hist_Limits[0], Hist_Limits[1]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(YAG, min_offset = 0, max_offset = 0))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nYAG screen projection" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('y ($\\mu$m)', fontsize=fsize)
plt.axis([Hist_Limits[0], Hist_Limits[1], Hist_Limits[0], Hist_Limits[1]])
plt.show()


#%% Test beam projection class - xz space projection

BP_xz = BeamProjection(x, z, w)

xz_project, hist_limits, X_Bins, Z_Bins = BP_xz.beam_projection(r_res=0.2, r_round=10, independant_bins=True, equal_bins=False)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(xz_project), cmap=cmap_white(plt.cm.spring), extent=([hist_limits[2], hist_limits[3], hist_limits[0], hist_limits[1]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(xz_project, min_offset = 0, max_offset = 1))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nUsing class" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x ($\\mu$m)', fontsize=fsize)
plt.axis([hist_limits[2], hist_limits[3], hist_limits[0], hist_limits[1]])
plt.show()


#%% Test beam projection class - xz space projection with predetermined histogram limits

window_limits = magnitude_conversion(np.array([info_Ex.rmin, info_Ex.rmax, info_Ex.zmin, info_Ex.zmax]), "", "micro")

print("Simulation window limits:")
print("zmin = %0.2f um" % np.round(info_Ex.zmin*1e6,2) )
print("zmax = %0.2f um" % np.round(info_Ex.zmax*1e6,2) )
print("rmin = %0.2f um" % np.round(info_Ex.rmin*1e6,2) )
print("rmax = %0.2f um" % np.round(info_Ex.rmax*1e6,2) )


xz_project_full, _, _ = BP_xz.beam_projection_fixed_window(window_extent = window_limits, r_res=0.3)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(xz_project_full), cmap=cmap_white(plt.cm.spring), extent=([window_limits[2], window_limits[3], window_limits[0], window_limits[1]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(xz_project_full, min_offset = 0, max_offset = 1))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nUsing class. Full simulation window" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x ($\\mu$m)', fontsize=fsize)
plt.axis([window_limits[2], window_limits[3], window_limits[0], window_limits[1]])
plt.show()


#%% Do I need the separate whole window function? Can I achieve the same thing by separating the extent of the imshow and plot limits?
# This looks like a better way of doing this!

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(xz_project), cmap=cmap_white(plt.cm.plasma), extent=([hist_limits[2], hist_limits[3], hist_limits[0], hist_limits[1]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(plt_limits_log(xz_project, min_offset = 0, max_offset = 1))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nWhole window plot using standard function" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x ($\\mu$m)', fontsize=fsize)
plt.axis([window_limits[2], window_limits[3], window_limits[0], window_limits[1]])
plt.show()