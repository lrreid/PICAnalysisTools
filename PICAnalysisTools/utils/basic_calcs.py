"""

Basic calculations


TO DO:
    - Add docstrings

"""

import numpy as np
from scipy.constants import pi
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_area


def area_circle(diameter, dia_unit: str = "milli", area_unit: str = "centi"):

    diameter_SI = magnitude_conversion(diameter, dia_unit, "", reciprocal_units = False)

    area = pi * (diameter_SI*0.5)**2
    return magnitude_conversion_area(area, "", area_unit, reciprocal_units = False)



def area_circle_projection(diameter, angle_deg, dia_unit: str = "milli", area_unit: str = "centi"):

    diameter_SI = magnitude_conversion(diameter, dia_unit, "", reciprocal_units = False)
    
    angle_rad   = angle_deg * (pi/180)
    dia_project = diameter_SI / np.cos(angle_rad)
    area        = pi * (diameter_SI*0.5) * (dia_project*0.5)
    
    return magnitude_conversion_area(area, "", area_unit, reciprocal_units = False)