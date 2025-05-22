"""
This file contains functions relevant for analysing the fields in PIC simualtions.


Date created: 10/01/2024 @ 22:26
Authors: Lewis R Reid

Inspiration for this class is "PIC_View_save_ne.py" it is the most
sophisticated I have in terms of options.

TO DO:
    - Check that the functions I have here can provide the correct information for scripts like "PIC_View_save_Field_plasma_inc_beam_loop.py" which I use a lot.
    - Add options for different normalisations of the z (propagation) axis. To laser centroid, plasma wavelength and edge of simulation box. Other options?

"""
import numpy as np
from scipy.constants import c, e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol
from PICAnalysisTools.utils.statistics import D4S_centroid, find_nearest

class FieldProperites():

    def __init__(self, field, info_field, z_unit: str = "micro", r_unit: str = "micro"):
        self.field      = field
        self.info_field = info_field
        self.z_unit     = z_unit
        self.r_unit     = r_unit


    def find_pixel_number(self, array, position, position_unit: str = "micro"):

        position_SI = magnitude_conversion(position, position_unit, "")
        pixel_no, _ = find_nearest(array, position_SI)

        return pixel_no

    def get_longitudinal_lineout(self, position, position_unit: str = "micro"):

        if position == 0:
            lineout  = self.field[int(len(self.info_field.r))//2,:]
        else:
            pixel_no = self.find_pixel_number(self.info_field.r, position, position_unit)
            lineout  = self.field[pixel_no,:]

        return lineout
    
    def get_transverse_lineout(self, position, position_unit: str = "micro"):
        
        pixel_no = self.find_pixel_number(self.info_field.z, position, position_unit)
        lineout  = self.field[:,pixel_no]

        return lineout

    def find_field_max(self, use_absolute: bool = False, centroid_unit: str = "micro"):

        if use_absolute is True:
            peak          = np.where(abs(self.field) == np.max(abs(self.field)))
        else:
            peak          = np.where(self.field == np.max(self.field))

        peak_z        = self.info_field.z[peak[1][0]]
        peak_r        = self.info_field.r[peak[0][0]]

        return magnitude_conversion(peak_z, "", centroid_unit), magnitude_conversion(peak_r, "", centroid_unit), peak[1][0], peak[0][0]
    
    def find_field_centroid(self, centroid_unit: str = "micro"):

        centroid_z_px, centroid_r_px = D4S_centroid(abs(self.field), rtn_int = True)
        centroid_z = self.info_field.z[centroid_z_px]
        centroid_r = self.info_field.r[centroid_r_px]

        return magnitude_conversion(centroid_z, "", centroid_unit), magnitude_conversion(centroid_r, "", centroid_unit), centroid_z_px, centroid_r_px
    


class PlasmaField():

    def __init__(self, ts, Species_name: str, den_unit: str = "centi"):
        self.ts           = ts
        self.Species_name = Species_name
        self.den_unit     = den_unit


    def get_plasma_density_map(self, Snapshot):

        rho, info_rho = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='%s' % self.Species_name, m='all')
        den           = np.abs( (-1/e)*rho )         # Plasma density (m^-3)

        return magnitude_conversion_vol(den, "", self.den_unit, reciprocal_units = True), info_rho
    

def get_focusing_field_map(ts, Snapshot):
    # get focusing field of plasma wave

    Ex0, info_field = ts.get_field( iteration=ts.iterations[Snapshot], field='E', m=0, coord='x')                     # Extract Ex field of plasma wave
    By, _           = ts.get_field( iteration=ts.iterations[Snapshot], field='B', m=0, coord='y')                     # Extract By field of plasma wave
    focusing_field = -1*(Ex0-(c*By))

    return focusing_field, info_field