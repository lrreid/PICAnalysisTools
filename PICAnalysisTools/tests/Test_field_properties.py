"""
Script to test functions related to analysing field data.

To use, you must set your own FolderPath and Simulation to simulaton data which includes a laser pulse.

TO DO:
    - Get plot limits from function

"""

import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries
from PICAnalysisTools.Field_Properties import FieldProperites, PlasmaField
from PICAnalysisTools.Laser_Properties import get_laser_cenroid
from PICAnalysisTools.Particle_Properties import BeamProjection
from PICAnalysisTools.utils.sim_path import set_sim_path, set_analysis_path
from PICAnalysisTools.Laser_Properties import get_a0_field_map
from PICAnalysisTools.utils.plot_limits import plt_limits_log, plt_limits_log_absolute
from PICAnalysisTools.utils.white_background_colormap import cmap_white
from PICAnalysisTools.utils.plot_limits import roundup
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion

fsize = 12

#%% Read data from files

Save_Plots          = False
show_laser_countour = True
show_beam           = True

plasma_species = "rho_plasma_elec"
beam_species   = 'electrons'

FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
Ana_name   = 'field_plotting_tests'
FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

if Save_Plots is True:
    Ana_dir = set_analysis_path(SimPath, Ana_name)
    DPI = 300
else:
    DPI = 150

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')

#%% Extract data from file

K = 2

Ex, info_Ex   = ts.get_field( iteration=ts.iterations[K], field = 'E', m=1, coord = 'x')
den, info_den = PlasmaField(ts, Species_name = plasma_species, den_unit = "centi").get_plasma_density_map(K)

ctau = np.round(ts.t[K]*c*1e6,0)

if show_laser_countour is True:
    a0_field, a0_max = get_a0_field_map(ts, K)

if show_beam is True:

    window_limits = magnitude_conversion(np.array([info_Ex.rmin, info_Ex.rmax, info_Ex.zmin, info_Ex.zmax]), "", "micro")
    x, y, z, w    = ts.get_particle( ['x', 'y', 'z', 'w'], species=beam_species, iteration=ts.iterations[K], plot=False)
    beam_dist     = BeamProjection(x, z, w).beam_projection_fixed_window(window_extent = window_limits, r_res=0.2)


#%% Take lineouts of field and plot the data

FP = FieldProperites(Ex, info_Ex, z_unit = "micro", r_unit = "micro")               # Create instance of FieldProperties class for laser Ex field

centroid_z, centroid_r, _, _ = FP.find_field_centroid()                             # Find centroid of laser field
peak_z, peak_r, _, _         = FP.find_field_max(use_absolute = True)               # Find location of peak of laser field 
laser_centroid, laser_peak   = get_laser_cenroid(ts, K, centroid_unit = "micro")    # Find centroid of laser in Z axis.

Ex_long  = FP.get_longitudinal_lineout(0)                                           # Extract longitudinal lineout through r = 0
Ex_trans = FP.get_transverse_lineout(laser_centroid, position_unit = "micro")       # Extract transverse lineout through laser centroid

#%% Plotting

Ex_plt_max = plt_limits_log_absolute(Ex, offset=0)                                  # calculate plot limits for Ex field plotting

fig, ax = plt.subplots(dpi=150)
plt.imshow(Ex,cmap=plt.cm.PRGn, extent=(np.round(info_Ex.imshow_extent*1e6,0)), aspect='auto')
plt.plot(centroid_z, centroid_r, '+', color='red', markersize=4, label="centroid")
plt.plot(peak_z, peak_r, '+', color='blue', markersize=4, label="peak")
plt.colorbar().set_label(label=('E$_{x}$ ($Vm^{-1}$)'), size=fsize)
plt.clim(-Ex_plt_max, Ex_plt_max)
plt.axis(np.round(info_Ex.imshow_extent*1e6,0))
ax.set_title("Laser $E_{x} field$\nc$\\tau$ = %d $\\mu$m" % ctau, fontsize=fsize)
ax.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax.set_ylabel('r ($\\mu$m)', fontsize=fsize)
plt.legend(loc='lower left', prop={'size': 10})

if Save_Plots is True:
    plt.savefig("%s/Ex_Field_ctau_%dum.png" % (Ana_dir, ctau),dpi=DPI,format="png")
    plt.close()
else:
    plt.show()


#%% Plot longitudinal lineout

fig, ax = plt.subplots(dpi=150)             # Plot line out
plt.plot(info_Ex.z*1e6, Ex_long, color='seagreen', linewidth=1, label=("$E_{x}$"))
plt.axis([np.round(info_Ex.zmin*1e6,0), np.round(info_Ex.zmax*1e6,0), -Ex_plt_max, Ex_plt_max])          # set axes [xmin, xmax, ymin, ymax]
ax.set_title("Longitudional lineout along\nr = 0 $\\mu$m", fontsize=fsize)
ax.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax.set_ylabel('E$_{x}$ ($Vm^{-1}$)', fontsize=fsize)
plt.grid(True)
# plt.legend(loc='lower left', prop={'size': 10})
if Save_Plots is True:
    plt.savefig("%s/Longitudinal_lineout_ctau_%dum.png" % (Ana_dir, ctau), dpi=DPI, format="png")
    plt.close()
else:
    plt.show()
    
#%% Plot longitudinal lineout

plt_trans_min, plt_trans_max = plt_limits_log(Ex_trans)

fig, ax = plt.subplots(dpi=150)             # Plot line out
plt.plot(info_Ex.r*1e6, Ex_trans, color='firebrick', linewidth=1, label=("$E_{x}$"))
plt.axis([np.round(info_Ex.rmin*1e6,0), np.round(info_Ex.rmax*1e6,0), plt_trans_min, plt_trans_max])          # set axes [xmin, xmax, ymin, ymax]
ax.set_title("Transverse lineout through\nlaser centroid", fontsize=fsize)
ax.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax.set_ylabel('E$_{x}$ ($Vm^{-1}$)', fontsize=fsize)
plt.grid(True)
# plt.legend(loc='lower left', prop={'size': 10})s
if Save_Plots is True:
    plt.savefig("%s/Transverse_lineout_ctau_%dum.png" % (Ana_dir, ctau), dpi=DPI, format="png")
    plt.close()
else:
    plt.show()

#%% Plot plasma density

fig, ax = plt.subplots(dpi=150)
plt.imshow(den, cmap=plt.cm.viridis_r, extent=(np.round(info_den.imshow_extent*1e6,0)), aspect='auto') #, norm=LogNorm())
plt.colorbar().set_label(label='$n_{e}$ (cm$^{-3}$)', size=fsize)
plt.clim(0, 3e19)
plt.axis(np.round(info_den.imshow_extent*1e6,0))          # set axes [xmin, xmax, ymin, ymax]
ax.set_title("c$\\tau$ = %d $\\mu$m" % ctau, fontsize=fsize)
ax.set_xlabel('z ($\\mu$m)', fontsize=fsize)
ax.set_ylabel('r ($\\mu$m)', fontsize=fsize)

if show_laser_countour is True:
    # Add countour plot of laser a0
    CS = ax.contour(np.array(info_den.z)*1e6, np.array(info_den.r)*1e6, a0_field, [0.135*a0_max, 0.2*a0_max, 0.4*a0_max, 0.6*a0_max, 0.8*a0_max, 0.98*a0_max], cmap='cool')
    ax.clabel(CS, inline=0, fontsize=0)

if show_beam is True:
    beam_proj = plt.imshow(beam_dist[0], cmap=cmap_white(plt.cm.plasma), extent=(np.round(info_den.imshow_extent*1e6,0)), aspect='auto')
    beam_proj.set_clim(0, roundup(0.80*np.max(beam_dist[0]), 1))

fig.tight_layout()

if Save_Plots is True:
    plt.savefig("%s/Plasma_density_ctau_%dum.png" % (Ana_dir, ctau), dpi=DPI, format="png")
    plt.close()
else:
    plt.show()