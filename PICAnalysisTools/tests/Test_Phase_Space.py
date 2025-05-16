"""
Test class which gets the phase space of a particle beam


TO DO:
    - Work through changes in classes to these examples
    - Remove return of order and replace with plot_limits functions.
    - Collect all bin mins and maxes into single variable

"""

##################################################################################################################
##################################################################################################################
####################### SCRIPT NOT COMPATIBLE WITH CHANGES MADE TO FUNCTIONS AND CLASSES #########################
##################################################################################################################
##################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from openpmd_viewer import OpenPMDTimeSeries        # This should be used for pmd viewer version 1.x.
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter

from PICAnalysisTools.utils.white_background_colormap import cmap_white
from PICAnalysisTools.Particle_Properties import PhaseSpace, BeamProjection
from PICAnalysisTools.utils.particle_selection import radial_selection
from PICAnalysisTools.utils.sim_path import set_sim_path

# from matplotlib import use
# use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

FolderPath = r'D:\Lewis\fbpic\20250501_pwfa_at_CLARA'
Simulation = '20250501_sigz_2p2um_sigr_5um_kpLrms_pi4'
Species    = "beam"

FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')

K = -3               # snapshot to analyse

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species=Species,
                                                iteration=ts.iterations[K], plot=False)

# Take radial selection from full dataset
#x, y, z, ux, uy, uz, w = radial_selection(0, 5e-6, x, y, z, ux, uy, uz, w)

ctau = round(ts.t[K]*c*1e3,2)

Phase_Space = PhaseSpace(x, y, z, ux, uy, uz, w, r_unit = "micro", energy_unit= "mega", div_unit= "milli", time_unit= "femto")      # Create instance of transverse phase space class

E_Phase, Order, Spec_Min, Spec_Max, E_Bins, r_min, r_max, R_Bins, Row_sum, Z_line, Col_sum = Phase_Space.Energy_z_space(z_res=0.1, z_round=5, Spec_Res=0.2, E_Round=50,
                                                                                                                        Centre_z=True, Find_Lineouts=True, lineout_height=0.2)


#%% Plotting

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(E_Phase), cmap=cmap_white(plt.cm.inferno), extent=([r_min, r_max, Spec_Min, Spec_Max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Order+1))

plt.plot(Row_sum, E_Bins[:-1], color = 'blue', linewidth = 1)
plt.plot(Z_line[:-1], Col_sum, color = 'blue', linewidth = 1)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('E (MeV)', fontsize=fsize)
plt.axis([r_min, r_max, Spec_Min, Spec_Max])
plt.show()


#%% Energy time space

E_Phase_t, Order_t, Spec_Min, Spec_Max, E_Bins, t_min, t_max, T_Bins, Row_sum, T_line, Col_sum = Phase_Space.Energy_time_space(Spec_Res=0.2, E_Round=50, t_res=0.1, t_round=5, Find_Lineouts=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(E_Phase_t), cmap=cmap_white(plt.cm.inferno), extent=([t_min, t_max, Spec_Min, Spec_Max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Order_t+1))

plt.plot(Row_sum, E_Bins[:-1], color = 'forestgreen', linewidth = 1)
plt.plot(T_line[:-1], Col_sum, color = 'forestgreen', linewidth = 1)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('t (fs)', fontsize=fsize)
ax1.set_ylabel('E (MeV)', fontsize=fsize)
plt.axis([t_min, t_max, Spec_Min, Spec_Max])
plt.show()


#%% Test transverse phase space and create plot

Trans_phase = Phase_Space.Div_r_space('x', r_res=0.30, r_round=10, div_round=10, div_res=0.1,
                              Find_Lineouts=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(Trans_phase[0]), cmap=cmap_white(plt.cm.viridis), extent=([Trans_phase[4], Trans_phase[5], -Trans_phase[2], Trans_phase[2]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Trans_phase[1]+1))

plt.plot(Trans_phase[7], Trans_phase[3][:-1], color = 'red', linewidth = 1)
plt.plot(Trans_phase[8], Trans_phase[9], color = 'red', linewidth = 1)

ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x\' (mrad)', fontsize=fsize)
plt.axis([Trans_phase[4], Trans_phase[5], -Trans_phase[2], Trans_phase[2]])
plt.show()

#%% Test beam projection class - YAG screen projection

BP_xy = BeamProjection(x, y, w)

YAG, order, x_min, x_max, x_Bins = BP_xy.beam_projection(r_res=0.3, r_round=1, independant_bins=False, equal_bins=True)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(YAG), cmap=cmap_white(plt.cm.winter), extent=([x_min, x_max, x_min, x_max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(order+0))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('y ($\\mu$m)', fontsize=fsize)
plt.axis([x_min, x_max, x_min, x_max])
plt.show()


#%% Test beam projection class - xz space projection

BP_xz = BeamProjection(x, z, w)

xz_project, Order, Z_min, Z_max, Z_Bins, X_min, X_max, X_Bins = BP_xz.beam_projection(r_res=0.3, r_round=10, independant_bins=True, equal_bins=False)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(xz_project), cmap=cmap_white(plt.cm.spring), extent=([X_min, X_max, Z_min, Z_max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Order+1))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nUsing class" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax1.set_ylabel('x ($\\mu$m)', fontsize=fsize)
plt.axis([X_min, X_max, Z_min, Z_max])
plt.show()

