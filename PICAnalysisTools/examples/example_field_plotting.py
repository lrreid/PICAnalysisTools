

from openpmd_viewer import OpenPMDTimeSeries
from PICAnalysisTools.field_plotting import plt_plasma_field
from PICAnalysisTools.utils.sim_path import set_sim_path

plasma_species = "rho_plasma_elec"
beam_species   = 'electrons'

FolderPath = r'C:\Users\ryi76833\OneDrive - Science and Technology Facilities Council\Documents\Python_Programs\PICAnalysisTools\PICAnalysisTools'
Simulation = 'example_data'
FilePath, SimPath = set_sim_path(FolderPath, Simulation, boosted_frame=False)

ts = OpenPMDTimeSeries(FilePath, check_all_files=False, backend='openpmd-api')
K  = 2

#%% Run instance of function under test

particle_select = [None, None, None, None, 50, None] # Selection limits for showing particles on plot. [z_min, z_max, r_min, r_max, E_min, E_max]. Units are assumed to match: z_unit, r_unit, energy_unit

plt_plasma_field(ts=ts, snapshot=K, field="E", mode=0, coord="z", plasma_species=plasma_species, field_unit = "giga", log_scale = False,
                 r_unit = "micro", z_unit = "micro", crop_r = True, r_max = 60, z_norm = True, Z_norm_type = "laser", show_laser_countour = True,
                 show_beam = True, beam_species = None, particle_selection = True, selection_limits = particle_select, energy_unit = "Mega",
                 show_accel_lineout = False, show_ne_lineout = False, save_plots = False, SimPath = SimPath, Ana_name = "field_plots", fsize = 12)