"""
Warning! Function still under construction and testing.

Functions for plotting of plasma fields.

TO DO:
    - option for Accelerating field lineout
    - Option for plasma density lineout
    - option for normalisation of the r axis (to the laser spot or skin depth)
    - option for cropping in z
    - Create function to loop through all simulation snapshots. It should call plt_plasma_density

    
What if you want to save or view multiple fields? Most of this code, file reads etc would have to be done multiple times.
Should I split the function into lots of separate ones, move the plotting out the function and lego brick all the pieces together in a plotting script?

"""
from os import getcwd
import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from openpmd_viewer import OpenPMDTimeSeries
from PICAnalysisTools.Field_Properties import FieldProperites, PlasmaField, get_focusing_field_map
from PICAnalysisTools.Particle_Properties import BeamProjection, get_normalised_momentum
from PICAnalysisTools.Laser_Properties import get_a0_field_map, get_laser_cenroid
from PICAnalysisTools.utils.sim_path import set_sim_path, set_analysis_path
from PICAnalysisTools.utils.plot_limits import plt_limits_log, plt_limits_log_absolute, plt_limits_absolute
from PICAnalysisTools.utils.white_background_colormap import cmap_white
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, get_order_letter
from PICAnalysisTools.utils.plasma_calcs import PlasmaDen_Conversions
from PICAnalysisTools.utils.rounding import normalise

plasma_species = "rho_plasma_elec"
beam_species   = 'electrons'

FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='openpmd-api')
K  = 2

#%% Extract data from file

def plt_plasma_field(ts, snapshot, field, mode, coord, plasma_species = "rho", field_unit = "centi", ne0 = None, colour_lims = None, log_scale = False, r_unit = "micro", z_unit = "micro",
                     crop_r = False, r_max = 40, z_norm = False, Z_norm_type = "laser",  show_laser_countour = False, show_beam = False, beam_species = None, particle_selection = False, selection_limits = [None, None, None, None, None, None],
                     energy_unit = "Mega", show_accel_lineout = False, show_ne_lineout = False, save_plots = False, SimPath = None, Ana_name = "ne_plots", fsize = 12):

    if save_plots is True:
        # If an analysis path is not defined, the current working directory is chosen as a default.
        if SimPath is None:
            SimPath = getcwd()
            print("No directory specified for saved plots.\nFiles automatically saved to: %s" % SimPath )
        Ana_dir = set_analysis_path(SimPath, Ana_name)
        DPI = 300
    else:
        DPI = 150

    ctau = magnitude_conversion(ts.t[snapshot]*c, "", z_unit)

    if field == "plasma":
        field_map, info_field = PlasmaField(ts, Species_name = plasma_species, den_unit = field_unit).get_plasma_density_map(snapshot)
        cbar_label = 'n$_{e}$ (%sm$^{-3}$)' % get_order_letter(field_unit)
        colour_map = plt.cm.viridis_r
    elif field == "focus":
        field_map, info_field = get_focusing_field_map(ts = ts, Snapshot = snapshot, field_unit = field_unit)
        cbar_label = '-($E_{x}$-$cB_{y})$ (%sVm$^{-1}$)' % get_order_letter(field_unit)
        colour_map = plt.cm.bwr
    else:
        field_map, info_field = ts.get_field( iteration=ts.iterations[snapshot], field=field, m=mode, coord=coord)             # Extract data for field
        field_map  = magnitude_conversion(field_map, "", field_unit)
        cbar_label = '$%s_{%s}$ (%sVm$^{-1}$)' % (field, coord, get_order_letter(field_unit))
        colour_map = plt.cm.PRGn_r

    if colour_lims is None:
        if field == "plasma":
            if log_scale is True:
                colour_lims = plt_limits_log(field_map, min_offset = 0, max_offset = -1)
                # print("Min ne: %0.2e" % np.min(field_map) )
                # print("Found limits: min: %0.2e \t max: %0.2e" % (colour_lims[0], colour_lims[1]) )
            else:
                colour_lims = np.array(plt_limits_log(field_map, min_offset = 1, max_offset = -1))*0.01
        else:
            if log_scale is True:
                colour_lims = plt_limits_log(field_map, min_offset = 0, max_offset = 0)
            else:
                colour_lims = plt_limits_absolute(field_map, 1)

    # If the name of a particle species is not specified, try to find one automatically and choose the zeroth. 
    # If none available, switch off show beam.
    if beam_species is None:
        try:
            beam_species = ts.avail_species[0]
        except:
            print("No particle species available")
            show_beam = False

    #%% Set plot limits and

    # Get r limits
    r_data = magnitude_conversion(info_field.r, "", r_unit)
    if crop_r is False or r_max is None:
        r_max = magnitude_conversion(info_field.rmax, "", r_unit)

    # Get z limits
    if ne0 is None and Z_norm_type == "l_p":
        if field == "plasma":
            den_unit = field_unit
        else:
            den_unit = "centi"

        if field == "plasma":
            ne0 = field_map[int(len(info_field.r)//2),-1] # Assume plasma density is the value at the far right of the simualtion box at r = 0.
        else:
            ne0 = get_background_density(ts, snapshot, Species_name = plasma_species, den_unit = den_unit)

    if z_norm is True:
        if Z_norm_type == "l_p":
            z_data = normalise_z(ts, snapshot, norm_type = Z_norm_type, n_e = ne0, z_unit = z_unit, den_unit = den_unit) # Need to think about units for plasma density. I can't use field_unit here. 
        else:
            z_data = normalise_z(ts, snapshot, norm_type = Z_norm_type, z_unit = z_unit)
    else:
        z_data = magnitude_conversion(info_field.z, "", z_unit)

    z_min = np.min(z_data)
    z_max = np.max(z_data)

    extent_converted   = np.array([z_min, z_max, -r_max, r_max])
    extent_full_window = np.array([z_min, z_max, magnitude_conversion(info_field.rmin, "", r_unit), magnitude_conversion(info_field.rmax, "", r_unit)])

    #%% Plot field

    title   = "c$\\tau$ = %0.2f %sm" % (ctau, "$\\mu$" if get_order_letter(z_unit) == "u" else get_order_letter(z_unit))
    if z_norm is True and Z_norm_type == "l_p":
        x_label = "z ($\\times\\lambda_{p}$)"
    else:
        x_label = "z (%sm)" % ("$\\mu$" if get_order_letter(z_unit) == "u" else get_order_letter(z_unit))
    y_label = "r (%sm)" % ("$\\mu$" if get_order_letter(r_unit) == "u" else get_order_letter(r_unit))

    fig, ax = plt.subplots(dpi=150)

    if log_scale is True:
        plt.imshow(field_map, cmap=colour_map, extent=(extent_full_window), aspect='auto', norm=LogNorm())
    else:
        plt.imshow(field_map, cmap=colour_map, extent=(extent_full_window), aspect='auto')
    plt.colorbar().set_label(label=cbar_label, size=fsize)
    plt.clim(colour_lims)   
    plt.axis(extent_converted)                                      # set axes [xmin, xmax, ymin, ymax]
    ax.set_title("%s" % title, fontsize=fsize)
    ax.set_xlabel('%s' % x_label, fontsize=fsize)
    ax.set_ylabel('%s' % y_label, fontsize=fsize)

    if show_accel_lineout is True:
        if field == "focus":
            LO_color   = "seagreen" 
        else:
            LO_color   = "red"
        E_Accel_LO = get_Accel_lineout(ts, snapshot)
        ax.plot(z_data, normalise(E_Accel_LO)*abs(extent_converted[2]), color=LO_color, linewidth=1) 

    if show_ne_lineout is True:
        if field == "focus":
            LO_color   = "magenta" 
        else:
            LO_color   = "blue"
        ne_lineout = get_ne_lineout(ts, snapshot, plasma_species)
        ax.plot(z_data, normalise(ne_lineout)*abs(extent_converted[2]), color=LO_color, linewidth=1) 

    if show_laser_countour is True:
        a0_field, a0_max = get_a0_field_map(ts, snapshot)
        a0_contour       = ax.contour(z_data, r_data, a0_field, [0.135*a0_max, 0.2*a0_max, 0.4*a0_max, 0.6*a0_max, 0.8*a0_max, 0.98*a0_max], cmap='cool')
        ax.clabel(a0_contour, inline=0, fontsize=0)

    if show_beam is True:
        beam_dist = get_beam_distribution(ts, info_field, snapshot, particle_selection, selection_limits, z_unit, r_unit, energy_unit)
        beam_proj = plt.imshow(beam_dist, cmap=cmap_white(plt.cm.plasma), extent=(extent_full_window), aspect='auto')
        beam_proj.set_clim(0, plt_limits_log_absolute(beam_dist, offset = -1)*0.2)

    fig.tight_layout()

    if save_plots is True:
        if field == "focus" or "plasma":
            plt.savefig("%s/%s_field_ctau_%dum.png" % (Ana_dir, field, ctau), dpi=DPI, format="png")
        else:
            plt.savefig("%s/%s_%s%s_field_ctau_%dum.png" % (Ana_dir, field, coord, mode, ctau), dpi=DPI, format="png")
        plt.close()
    else:
        plt.show()

    return None


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


def get_background_density(ts, snapshot, Species_name, den_unit: str = "centi"):
    # This method is not robust against simulations using ADK ionisation
    # Background density is taken from a single pixel at the far right of the simulation window at r = 0.

    field_map, info_field = PlasmaField(ts, Species_name = Species_name, den_unit = den_unit).get_plasma_density_map(snapshot)
    ne0 = field_map[int(len(info_field.r)//2),-1]

    return ne0

def get_Accel_lineout(ts, snapshot, unit = "Giga"):
    E_Accel, info_E_Accel = ts.get_field( iteration=ts.iterations[snapshot], field="E", m=0, coord="z") 
    E_Accel_LO = FieldProperites(E_Accel, info_E_Accel).get_longitudinal_lineout(0)

    return magnitude_conversion(E_Accel_LO, "", unit)

def get_ne_lineout(ts, snapshot, plasma_species, unit = "centi"):
    ne_map, info_ne = PlasmaField(ts, Species_name = plasma_species, den_unit = unit).get_plasma_density_map(snapshot)
    ne_lineout      = FieldProperites(ne_map, info_ne).get_longitudinal_lineout(0)

    return ne_lineout

def get_beam_distribution(ts, info_field, snapshot, particle_selection, selection_limits, z_unit, r_unit, energy_unit):
    if particle_selection is True:
        if selection_limits[0] is None:
            selection_limits[0] = info_field.zmin
        else:
            selection_limits[0] = magnitude_conversion(selection_limits[0], z_unit, "")
        if selection_limits[1] is None:
            selection_limits[1] = info_field.zmax
        else:
            selection_limits[1] = magnitude_conversion(selection_limits[1], z_unit, "")
        if selection_limits[2] is None:
            selection_limits[2] = info_field.rmin
        else:
            selection_limits[2] = magnitude_conversion(selection_limits[2], r_unit, "")
        if selection_limits[3] is None:
            selection_limits[3] = info_field.rmax
        else:
            selection_limits[3] = magnitude_conversion(selection_limits[3], r_unit, "")

        if selection_limits[4] is not None:
            selection_limits[4] = get_normalised_momentum(selection_limits[4], energy_unit)
        if selection_limits[5] is not None:
            selection_limits[5] = get_normalised_momentum(selection_limits[5], energy_unit)

        x, z, w = ts.get_particle( ['x', 'z', 'w'], species=beam_species, iteration=ts.iterations[snapshot], plot=False, select={'x':[selection_limits[2],selection_limits[3]], 'y':[selection_limits[2],selection_limits[3]], 'z':[selection_limits[0],selection_limits[1]], 'uz':[selection_limits[4],selection_limits[5]]})
    else:
        x, z, w = ts.get_particle( ['x', 'z', 'w'], species=beam_species, iteration=ts.iterations[snapshot], plot=False)

    window_limits = magnitude_conversion(np.array([info_field.rmin, info_field.rmax, info_field.zmin, info_field.zmax]), "", "micro")
    beam_dist     = np.flipud(BeamProjection(x, z, w).beam_projection_fixed_window(window_extent = window_limits, r_res=magnitude_conversion(info_field.dr, "", "micro"))[0])

    return beam_dist

#%% Run instance of function under test

particle_select = [None, None, None, None, 50, None] # Selection limits for showing particles on plot. [z_min, z_max, r_min, r_max, E_min, E_max]. Units are assumed to match: z_unit, r_unit, energy_unit

plt_plasma_field(ts=ts, snapshot=K, field="E", mode=0, coord="z", plasma_species=plasma_species, field_unit = "giga", log_scale = False,
                 r_unit = "micro", z_unit = "micro", crop_r = True, r_max = 60, z_norm = True, Z_norm_type = "laser", show_laser_countour = True,
                 show_beam = True, beam_species = None, particle_selection = True, selection_limits = particle_select, energy_unit = "Mega",
                 show_accel_lineout = False, show_ne_lineout = False, save_plots = False, SimPath = SimPath, Ana_name = "field_plots", fsize = 12)
