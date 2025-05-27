"""
Functions for plotting of plasma fields.

TO DO:
    - option for Accelerating field lineout
    - Option for plasma density lineout
    - Option for normalisation of the z axis (Centre on laser centroid or normalise to plasma density) - This should be in its own function
    - option for normalisation of the r axis (to the laser spot or skin depth)
    - Option for cropping in r
    - option for cropping in z
    - option for defining plasma density limits
    - option for chaning units in z axis and r axis. (mm would be really useful for z axis)
    - option for choosing units of ctau in file name and plot title
    - Modify set_Analysis_path to have optional inputs.
    - options for particle species selection in space and energy

    - Need to think about how units are used. Particularly for r_max

    - Create function to loop through all simulation snapshots. It should call plt_plasma_density

"""
from os import getcwd
import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from openpmd_viewer import OpenPMDTimeSeries
from PICAnalysisTools.Field_Properties import PlasmaField
from PICAnalysisTools.Particle_Properties import BeamProjection
from PICAnalysisTools.utils.sim_path import set_sim_path, set_analysis_path
from PICAnalysisTools.Laser_Properties import get_a0_field_map
from PICAnalysisTools.utils.plot_limits import plt_limits_log, plt_limits_log_absolute
from PICAnalysisTools.utils.white_background_colormap import cmap_white
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion

plasma_species = "rho_plasma_elec"
beam_species   = 'electrons'

FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')
K  = 2

#%% Extract data from file

def plt_plasma_denstiy(ts, snapshot, plasma_species, den_unit = "centi", ne_lims = None, log_scale = False, r_unit = "micro", z_unit = "micro", crop_r = False, r_max = 40, show_laser_countour = False, show_beam = False, beam_species = None, save_plots = False, SimPath = None, Ana_name = "ne_plots", fsize = 12):

    if save_plots is True:
        # If an analysis path is not defined, the current working directory is chosen as a default.
        if SimPath is None:
            SimPath = getcwd()
            print("No directory specified for saved plots.\nFiles automatically saved to: %s" % SimPath )
        Ana_dir = set_analysis_path(SimPath, Ana_name)
        DPI = 300
    else:
        DPI = 150

    den, info_den = PlasmaField(ts, Species_name = plasma_species, den_unit = den_unit).get_plasma_density_map(snapshot)

    ctau = magnitude_conversion(ts.t[snapshot]*c, "", z_unit)

    if ne_lims is None:
        if log_scale is True:
            ne_lims = plt_limits_log(den, min_offset = 0, max_offset = -1)
            # print("Min ne: %0.2e" % np.min(den) )
            # print("Found limits: min: %0.2e \t max: %0.2e" % (ne_lims[0], ne_lims[1]) )
        else:
            ne_lims = np.array(plt_limits_log(den, min_offset = 1, max_offset = -1))*0.01

    if show_laser_countour is True:
        a0_field, a0_max = get_a0_field_map(ts, snapshot)

    # If the name of a particle species is not specified, try to find one automatically and choose the zeroth. 
    # If none available, switch off show beam.
    if beam_species is None:
        try:
            beam_species = ts.avail_species[0]
        except:
            print("No particle species available")
            show_beam = False

    if show_beam is True:

        window_limits = magnitude_conversion(np.array([info_den.rmin, info_den.rmax, info_den.zmin, info_den.zmax]), "", "micro")
        x, z, w       = ts.get_particle( ['x', 'z', 'w'], species=beam_species, iteration=ts.iterations[snapshot], plot=False)
        beam_dist     = BeamProjection(x, z, w).beam_projection_fixed_window(window_extent = window_limits, r_res=0.2)


    if crop_r is False or r_max is None:
        r_max = magnitude_conversion(info_den.rmax, "", r_unit)

    # imshow_extent = 

    #%% Plot plasma density

    fig, ax = plt.subplots(dpi=150)

    if log_scale is True:
        plt.imshow(den, cmap=plt.cm.viridis_r, extent=(np.round(info_den.imshow_extent*1e6,0)), aspect='auto', norm=LogNorm())
    else:
        plt.imshow(den, cmap=plt.cm.viridis_r, extent=(np.round(info_den.imshow_extent*1e6,0)), aspect='auto') #, norm=LogNorm())
    plt.colorbar().set_label(label='n$_{e}$ (cm$^{-3}$)', size=fsize)
    plt.clim(ne_lims)
    # plt.axis(np.round(info_den.imshow_extent*1e6,0))          # set axes [xmin, xmax, ymin, ymax]
    plt.axis([info_den.zmin*1e6, info_den.zmax*1e6, -r_max, r_max])
    ax.set_title("c$\\tau$ = %d $\\mu$m" % ctau, fontsize=fsize)
    ax.set_xlabel('z ($\\mu$m)', fontsize=fsize)
    ax.set_ylabel('r ($\\mu$m)', fontsize=fsize)

    if show_laser_countour is True:
        # Add countour plot of laser a0
        CS = ax.contour(np.array(info_den.z)*1e6, np.array(info_den.r)*1e6, a0_field, [0.135*a0_max, 0.2*a0_max, 0.4*a0_max, 0.6*a0_max, 0.8*a0_max, 0.98*a0_max], cmap='cool')
        ax.clabel(CS, inline=0, fontsize=0)

    if show_beam is True:
        beam_proj = plt.imshow(beam_dist[0], cmap=cmap_white(plt.cm.plasma), extent=(np.round(info_den.imshow_extent*1e6,0)), aspect='auto')
        beam_proj.set_clim(0, plt_limits_log_absolute(beam_dist[0], offset = -1)*0.2)

    fig.tight_layout()

    if save_plots is True:
        plt.savefig("%s/Plasma_density_ctau_%dum.png" % (Ana_dir, ctau), dpi=DPI, format="png")
        plt.close()
    else:
        plt.show()

    return None


#%% Run instance of function under test

plt_plasma_denstiy(ts=ts, snapshot=K, plasma_species=plasma_species, den_unit = "centi", log_scale = False, crop_r = True, r_max = 60, show_laser_countour = True, show_beam = True, beam_species = None, save_plots = False, SimPath = SimPath, Ana_name = "ne_plots", fsize = 12)
