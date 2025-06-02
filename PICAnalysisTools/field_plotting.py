"""
Warning! Function still under construction and testing.

Functions for plotting of plasma fields.

TO DO:
    - option for Accelerating field lineout
    - Option for plasma density lineout
    - option for normalisation of the r axis (to the laser spot or skin depth)
    - option for cropping in z
    - options for particle species selection in space and energy
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
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, get_order_letter
from PICAnalysisTools.Laser_Properties import get_laser_cenroid
from PICAnalysisTools.utils.plasma_calcs import PlasmaDen_Conversions
from PICAnalysisTools.Particle_Properties import get_normalised_momentum

plasma_species = "rho_plasma_elec"
beam_species   = 'electrons'

FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='h5py')
K  = 2

#%% Extract data from file

def plt_plasma_denstiy(ts, snapshot, plasma_species, den_unit = "centi", ne0 = None, ne_lims = None, log_scale = False, r_unit = "micro", z_unit = "micro", crop_r = False, r_max = 40, z_norm = False, Z_norm_type = "laser",  show_laser_countour = False, show_beam = False, beam_species = None, save_plots = False, SimPath = None, Ana_name = "ne_plots", fsize = 12):

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


    #%% Set plot limits and

    # Get r limits
    r_data = magnitude_conversion(info_den.r, "", r_unit)
    if crop_r is False or r_max is None:
        r_max = magnitude_conversion(info_den.rmax, "", r_unit)

    # Get z limits
    if ne0 is None:
        ne0 = den[int(len(info_den.r)//2),-1] # Assume plasma density is the value at the far right of the simualtion box at r = 0.

    if z_norm is True:
        z_data = normalise_z(ts, snapshot, norm_type = Z_norm_type, n_e = ne0, z_unit = z_unit)
    else:
        z_data = magnitude_conversion(info_den.z, "", z_unit)

    z_min = np.min(z_data)
    z_max = np.max(z_data)

    extent_converted   = np.array([z_min, z_max, -r_max, r_max])
    extent_full_window = np.array([z_min, z_max, magnitude_conversion(info_den.rmin, "", r_unit), magnitude_conversion(info_den.rmax, "", r_unit)])


    #%% Plot plasma density

    title   = "c$\\tau$ = %0.2f %sm" % (ctau, "$\\mu$" if get_order_letter(z_unit) == "u" else get_order_letter(z_unit))
    if z_norm is True and Z_norm_type == "l_p":
        z_label = "z ($\\times\\lambda_{p}$)"
    else:
        z_label = "z (%sm)" % ("$\\mu$" if get_order_letter(z_unit) == "u" else get_order_letter(z_unit))
    r_label = "r (%sm)" % ("$\\mu$" if get_order_letter(r_unit) == "u" else get_order_letter(r_unit))


    fig, ax = plt.subplots(dpi=150)

    if log_scale is True:
        plt.imshow(den, cmap=plt.cm.viridis_r, extent=(extent_full_window), aspect='auto', norm=LogNorm())
    else:
        plt.imshow(den, cmap=plt.cm.viridis_r, extent=(extent_full_window), aspect='auto')
    plt.colorbar().set_label(label='n$_{e}$ (cm$^{-3}$)', size=fsize)
    plt.clim(ne_lims)   
    plt.axis(extent_converted)                                      # set axes [xmin, xmax, ymin, ymax]
    ax.set_title("%s" % title, fontsize=fsize)
    ax.set_xlabel('%s' % z_label, fontsize=fsize)
    ax.set_ylabel('%s' % r_label, fontsize=fsize)

    if show_laser_countour is True:
        # Add countour plot of laser a0
        CS = ax.contour(z_data, r_data, a0_field, [0.135*a0_max, 0.2*a0_max, 0.4*a0_max, 0.6*a0_max, 0.8*a0_max, 0.98*a0_max], cmap='cool')
        ax.clabel(CS, inline=0, fontsize=0)

    if show_beam is True:
        beam_proj = plt.imshow(beam_dist[0], cmap=cmap_white(plt.cm.plasma), extent=(extent_full_window), aspect='auto')
        beam_proj.set_clim(0, plt_limits_log_absolute(beam_dist[0], offset = -1)*0.2)

    fig.tight_layout()

    if save_plots is True:
        plt.savefig("%s/Plasma_density_ctau_%dum.png" % (Ana_dir, ctau), dpi=DPI, format="png")
        plt.close()
    else:
        plt.show()

    return None


# move to a file in utils??
def normalise_z(ts, snapshot, norm_type: str = "laser", n_e = None, z_unit: str = "micro", den_unit: str = "centi"):
    # Is there a way to make this more robust? What if there is no E field data in the file or no x axis?
    # Can I use something like ts.avail_species[0]?
    # Have z data as an input rather than re-reading the data file?

    _, info_Ex     = ts.get_field( iteration=ts.iterations[snapshot], field='E', m=1, coord='x')
    z_data         = magnitude_conversion(info_Ex.z, "", z_unit)
    laser_centroid = get_laser_cenroid(ts, snapshot, centroid_unit = z_unit)[0]

    # Is there any more interesting normalisation types to have?
    if norm_type == 'l_p':
        lambda_p = PlasmaDen_Conversions(n_e, den_unit=den_unit, wavelength_unit=z_unit).lambda_p
        Z_norm = ( z_data - laser_centroid ) / lambda_p
    elif norm_type == 'laser':
        Z_norm = z_data - laser_centroid
    else:
        print("norm type incorrectly defined. Original simulation axis chosen as default.")
        Z_norm = z_data

    return Z_norm


#%% Run instance of function under test

plt_plasma_denstiy(ts=ts, snapshot=K, plasma_species=plasma_species, den_unit = "centi", log_scale = False, r_unit = "micro", z_unit = "milli", crop_r = True, r_max = 60, z_norm = True, Z_norm_type = "l_p", show_laser_countour = True, show_beam = True, beam_species = None, save_plots = False, SimPath = SimPath, Ana_name = "ne_plots", fsize = 12)
