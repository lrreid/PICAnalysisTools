"""
This file contains functions relevant for analysing the fields in PIC simualtions.


Date created: 10/01/2024 @ 22:26
Authors: Lewis R Reid

Inspiration for this class is "PIC_View_save_ne.py" it is the most
sophisticated I have in terms of options.

TO DO:
    - Check that the functions I have here can provide the correct information for scripts like "PIC_View_save_Field_plasma_inc_beam_loop.py" which I use a lot.
    - Add a better way of finding laser centroid rather than highest intensity point.
    

"""
import numpy as np
from scipy.constants import c, m_e, e, pi
from scipy.signal import hilbert
from PICAnalysisTools.utils.unit_conversion import magnitude_conversion

class plasma_fields():

    def __init__(self, ts):
        self.ts = ts


    def find_pixel_number(self, array, position, array_unit="", position_unit="micro"):

        array_in_pos_units = magnitude_conversion(array, array_unit, position_unit)
        pixel_no           = np.where( np.round(array_in_pos_units,2) == position )

        return pixel_no

    def get_longitudinal_lineout(self, field, info_field, position, position_unit="micro"):

        if position == 0:
            lineout  = field[int(len(info_field.r))//2,:]
        else:
            pixel_no = self.find_pixel_number(info_field.r, position, "", position_unit)
            lineout  = field[pixel_no,:]

        return lineout
    
    def get_transverse_lineout(self, field, info_field, position, position_unit="micro"):
        
        pixel_no = self.find_pixel_number(info_field.z, position, "", position_unit)
        lineout  = field[:,pixel_no]

        return lineout
    

    def get_a0_field(self, Snapshot, central_wavelength):
        # This should be updated to get the central wavelength automatically from the laser diagnostics.
        Ex, _    = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='E', m=1, coord='x')
        a0_field = np.array((np.absolute(hilbert(np.array(Ex)))*e*central_wavelength)/(m_e*c*c*2*pi))    # Laser a0 field
        a0_max   = np.max(a0_field)

        return a0_field, a0_max
    

    def get_plasma_density(self, Snapshot, Species_name):

        rho, info_rho = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='%s' % Species_name, m='all')
        den           = np.abs( ((-1/e)*1e-6)*rho )         # Plasma density (cm^-3)

        return den, info_rho    


    def find_field_max(self, Snapshot, Field, M, Coord, centroid_unit="micro"):

        field, info_field = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field=Field, m=M, coord=Coord)
        centroid          = np.where(field == np.max(abs(field)))
        centroid_z        = info_field.z[centroid[1][0]]
        centroid_r        = info_field.r[centroid[0][0]]

        return magnitude_conversion(centroid_z, "", centroid_unit), magnitude_conversion(centroid_r, "", centroid_unit), centroid[1][0], centroid[0][0]
    
    def find_laser_centroid(self, Snapshot, centroid_unit="micro"):

        centroid_z, centroid_r, centroid_z_px, centroid_r_px = self.find_field_max(self, Snapshot, Field='E', M=1, Coord='x', centroid_unit=centroid_unit)

        return centroid_z, centroid_r, centroid_z_px, centroid_r_px
    
    def get_focusing_field(self, Snapshot):
        # get focusing field of plasma wave

        Ex0, _   = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='E', m=0, coord='x')                     # Extract Ex field of plasma wave
        By, _    = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='B', m=0, coord='y')                     # Extract By field of plasma wave
        focusing_field = -1*(Ex0-(c*By))

        return focusing_field