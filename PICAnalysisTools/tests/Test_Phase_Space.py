"""
Test class which gets the phase space of a particle beam
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
from white_background_colormap import cmap_white

from Particle_Properties import PhaseSpace, BeamProjection
import Useful_Functions as uf

from matplotlib import use
use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

#FilePath = r'D:\Lewis\fbpic\20231206_FEBE_Plasma_Channel_110TW\20231208_FEBE_Matched_100mm_0p53Lp_Lab\diags\hdf5'
FilePath = r'D:\Lewis\fbpic\20231206_FEBE_Plasma_Channel_110TW\20231208_100mm_Timing_scan\20231208_FEBE_Matched_100mm_0p55Lp\lab_diags\hdf5'
#FilePath = r'C:\Users\lewis\Desktop\plasma_channel\20240102_FEBE_Plasma_Channel_2170mJ_23fs_w0_30um_20um_10fs_5pC_1E17\lab_diags\hdf5'

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')

K = 12               # snapshot to analyse

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species='beam',
                                                iteration=ts.iterations[K], plot=False)

# Take radial selection from full dataset
x, y, z, ux, uy, uz, w = uf.radial_selection(0, 50e-6, x, y, z, ux, uy, uz, w)

ctau = round(ts.t[K]*c*1e3,2)

LPS = PhaseSpace(z*1e6, ux, uy, uz, w)      # Create instance of transverse phase space class
TPS = PhaseSpace(x*1e6, ux, uy, uz, w)      # Create instance of longitudinal phase space class

E_Phase, Order, Spec_Min, Spec_Max, E_Bins, r_min, r_max, R_Bins, Row_sum, Z_line, Col_sum = LPS.Energy_r_space(Spec_Res=0.1, E_Round=50, r_res=0.1, r_round=5,
                                                                                                                Normalise_r=True, Find_Lineouts=True, lineout_height=0.2, equal_r_bins=False)


#%% Plotting

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(E_Phase), cmap=cmap_white(plt.cm.inferno), extent=([r_min, r_max, Spec_Min, Spec_Max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Order+1))

plt.plot(Row_sum, E_Bins[:-1], color = 'blue', linewidth = 1)
plt.plot(Z_line[:-1], Col_sum, color = 'blue', linewidth = 1)

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\mu$m)', fontsize=fsize)
ax1.set_ylabel('E (MeV)', fontsize=fsize)
plt.axis([r_min, r_max, Spec_Min, Spec_Max])
plt.show()


#%% Energy time space

E_Phase_t, Order_t, Spec_Min, Spec_Max, E_Bins, t_min, t_max, T_Bins, Row_sum, T_line, Col_sum = LPS.Energy_time_space(Spec_Res=0.1, E_Round=50, t_res=0.1, t_round=5)

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

Trans_phase = TPS.Div_r_space('x', div_round=1, div_res=0.1, r_res=0.30, r_round=10,
                              Find_Lineouts=True, lineout_height=0.2)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(Trans_phase[0]), cmap=cmap_white(plt.cm.viridis), extent=([Trans_phase[4], Trans_phase[5], -Trans_phase[2], Trans_phase[2]]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Trans_phase[1]+1))

plt.plot(Trans_phase[7], Trans_phase[3][:-1], color = 'red', linewidth = 1)
plt.plot(Trans_phase[8], Trans_phase[9], color = 'red', linewidth = 1)

ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\mu$m)', fontsize=fsize)
ax1.set_ylabel('x\' (mrad)', fontsize=fsize)
plt.axis([Trans_phase[4], Trans_phase[5], -Trans_phase[2], Trans_phase[2]])
plt.show()

#%% Test beam projection class - YAG screen projection

BP_xy = BeamProjection(x*1e6, y*1e6, w)

YAG, order, x_min, x_max, x_Bins = BP_xy.beam_projection(r_res=0.3, r_round=1, independant_bins=False, equal_bins=True)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(YAG), cmap=cmap_white(plt.cm.winter), extent=([x_min, x_max, x_min, x_max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(order+0))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('x ($\mu$m)', fontsize=fsize)
ax1.set_ylabel('y ($\mu$m)', fontsize=fsize)
plt.axis([x_min, x_max, x_min, x_max])
plt.show()


#%% Test beam projection class - xz space projection

BP_xz = BeamProjection(x*1e6, z*1e6, w)

xz_project, Order, Z_min, Z_max, Z_Bins, X_min, X_max, X_Bins = BP_xz.beam_projection(r_res=0.3, r_round=10, independant_bins=True, equal_bins=False)

fig1, ax1 = plt.subplots()
plt.imshow(np.flipud(xz_project), cmap=cmap_white(plt.cm.spring), extent=([X_min, X_max, Z_min, Z_max]), aspect='auto', norm=LogNorm())
plt.colorbar().set_label(label='Number of electrons',size=fsize)
plt.clim(1, 10**(Order+1))

ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_title("c$\\tau$ = %0.2f mm\nUsing class" % round(ctau,2), fontsize=fsize)
ax1.set_xlabel('z ($\mu$m)', fontsize=fsize)
ax1.set_ylabel('x ($\mu$m)', fontsize=fsize)
plt.axis([X_min, X_max, Z_min, Z_max])
plt.show()

