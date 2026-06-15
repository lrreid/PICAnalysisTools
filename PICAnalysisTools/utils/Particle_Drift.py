
import numpy as np
from PICAnalysisTools.core.unit_conversions import magnitude_conversion
from PICAnalysisTools.utils.particle_selection import radial_selection

def simple_drift(x, y, ux, uy, uz, drift_dist, drift_unit: str = "milli"):

    drift_SI = magnitude_conversion(drift_dist, drift_unit, "")

    #%% Calculate drift of the electron beam in vacuum.
    px_pz = np.arctan2(ux, uz)
    py_pz = np.arctan2(uy, uz)

    # Calculate drift of electron beam in x and y
    x_drift = x+(px_pz*drift_SI)
    y_drift = y+(py_pz*drift_SI)

    return x_drift, y_drift



def drift_with_aperture(x, y, z, ux, uy, uz, w, drift_dist, ap_min, ap_max, drift_unit: str = "milli", ap_unit: str = "milli"):

    # Do the simple drift
    x_drift, y_drift = simple_drift(x=x, y=y, ux=ux, uy=uy, uz=uz, drift_dist=drift_dist, drift_unit = drift_unit)

    # Take radial selection based on aperture
    x_radial_select, y_radial_select, z_radial_select, ux_radial_select, uy_radial_select, uz_radial_select, w_radial_select = radial_selection(R_Select_min=ap_min, R_Select_max=ap_max, x=x_drift, y=y_drift, z=z, ux=ux, uy=uy, uz=uz, w=w, select_unit=ap_unit)

    return x_radial_select, y_radial_select, z_radial_select, ux_radial_select, uy_radial_select, uz_radial_select, w_radial_select
