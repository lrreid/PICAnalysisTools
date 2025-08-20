"""
This file contains functions relevant for analysing the fields in PIC simualtions.

TO DO:
    - Add options for different normalisations of the z (propagation) axis. To laser centroid, plasma wavelength and edge of simulation box. Other options?

"""
import numpy as np
from scipy.constants import c, e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol
from PICAnalysisTools.utils.statistics import D4S_centroid
from PICAnalysisTools.utils.rounding import find_nearest

class FieldProperites():

    def __init__(self, field, info_field, z_unit: str = "micro", r_unit: str = "micro"):
        """
        Class for extracting properties from laser and plasma fields.

        Parameters
        ----------
        field : array_like
            2D array containing field data
        info_field : openPMD FieldMetaInformation object
            object containing meta-information about the grid. Typically obtained via ts.get_field()
        z_unit : str, optional
            Order of magnitude of propagation longitudianal axis unit, by default from z_unit "milli", by default "micro"
        r_unit : str, optional
            Order of magnitude of propagation transverse axis unit, by default from z_unit "milli", by default "micro"
        """
        self.field      = field
        self.info_field = info_field
        self.z_unit     = z_unit
        self.r_unit     = r_unit


    def find_pixel_number(self, array, position, position_unit: str = "micro"):
        """
        Find the index in the field closest to a defined value

        Parameters
        ----------
        array : array
            array containing data
        position : float
            Target value that you want to find an array. Default unit: microns
        position_unit : str, optional
            Order of magnitude of target position unit, by default "micro"

        Returns
        -------
        pixel_no: int
            Index in array where value is closest to the target position.
        """

        position_SI = magnitude_conversion(position, position_unit, "")
        pixel_no, _ = find_nearest(array, position_SI)

        return pixel_no

    def get_longitudinal_lineout(self, position, position_unit: str = "micro", field_unit: str = ""):
        """
        Get longitudianal lineout of field.
        Function currently assumes cylindrical geometery.

        Parameters
        ----------
        position : float
            Radial position in field where the lineout should be taken
        position_unit : str, optional
            Order of magnitude of target position unit, by default "micro"

        Returns
        -------
        lineout: array
            Numpy array containing lineout of of field along chosen radial position
        """

        if position == 0:
            lineout  = self.field[int(len(self.info_field.r))//2,:]
        else:
            pixel_no = self.find_pixel_number(self.info_field.r, position, position_unit)
            lineout  = self.field[pixel_no,:]

        return magnitude_conversion(self.info_field.z, "", self.z_unit), magnitude_conversion(lineout, "", field_unit)
    
    def get_transverse_lineout(self, position, position_unit: str = "micro", field_unit: str = ""):
        """
        Get transverse lineout of field.

        Parameters
        ----------
        position : float
            Longitudinal position in field where the lineout should be taken
        position_unit : str, optional
            Order of magnitude of target position unit, by default "micro"

        Returns
        -------
        lineout: array
            Numpy array containing lineout of of field along chosen longitudinal position
        """
        
        pixel_no = self.find_pixel_number(self.info_field.z, position, position_unit)
        lineout  = self.field[:,pixel_no]

        return magnitude_conversion(self.info_field.r, "", self.r_unit), magnitude_conversion(lineout, "", field_unit)

    def find_field_max(self, use_absolute: bool = False, centroid_unit: str = "micro"):
        """
        Find the location of the maximum amplitude of the field

        Parameters
        ----------
        use_absolute : bool, optional
            Option to use absolute value of field to locate maximum. abs(field), by default False
        centroid_unit : str, optional
            Order of magnitude of centroid units, by default "micro"

        Returns
        -------
        peak_z: float
            Longitudinal component of field maximum. Default unit: microns
        peak_r: float
            Transverse component of field maximum. Default unit: microns
        centroid_z_px: int
            Index in field array at the field maximum in the longitudinal direction.
        centroid_r_px: int
            Index in field array at the field maximum in the transverse direction.
        """

        if use_absolute is True:
            peak = np.where(abs(self.field) == np.max(abs(self.field)))
        else:
            peak = np.where(self.field == np.max(self.field))

        peak_z = self.info_field.z[peak[1][0]]
        peak_r = self.info_field.r[peak[0][0]]

        return magnitude_conversion(peak_z, "", centroid_unit), magnitude_conversion(peak_r, "", centroid_unit), peak[1][0], peak[0][0]
    
    def find_field_centroid(self, use_absolute: bool = True, centroid_unit: str = "micro"):
        """
        Find the location of the centroid of the field

        Parameters
        ----------
        use_absolute : bool, optional
            Option to use absolute value of field to locate maximum. abs(field), by default False
        centroid_unit : str, optional
            Order of magnitude of centroid units, by default "micro"

        Returns
        -------
        centroid_z: float
            Longitudinal component of field centroid. Default unit: microns
        centroid_r: float
            Transverse component of field centroid. Default unit: microns
        centroid_z_px: int
            Index in field array at the field centroid in the longitudinal direction.
        centroid_r_px: int
            Index in field array at the field centroid in the transverse direction.
        """

        if use_absolute is True:
            centroid_z_px, centroid_r_px = D4S_centroid(abs(self.field), rtn_int = True)
        else:
            centroid_z_px, centroid_r_px = D4S_centroid(self.field, rtn_int = True)
        
        centroid_z = self.info_field.z[centroid_z_px]
        centroid_r = self.info_field.r[centroid_r_px]

        return magnitude_conversion(centroid_z, "", centroid_unit), magnitude_conversion(centroid_r, "", centroid_unit), centroid_z_px, centroid_r_px
    


class PlasmaField():

    def __init__(self, ts, Species_name: str, den_unit: str = "centi"):
        """
        Class to obtain the plasma (electron) density from a PIC simulation snapshot

        Parameters
        ----------
        ts : Object
            LpaDiagnostics "time series" object.
        Species_name : str
            Name of particle species which has charge density data.
        den_unit : str, optional
           Order of magnitude of plasma density unit, by default "centi"
        """

        self.ts           = ts
        self.Species_name = Species_name
        self.den_unit     = den_unit


    def get_plasma_density_map(self, Snapshot):
        """
        Get the plasma density map for a PIC simulation snapshot

        Parameters
        ----------
        Snapshot : int
            Snapshot from simulation to be analysed.

        Returns
        -------
        den: array
            2D array containing plasma density data. Default unit: cm^-3
        info_field : openPMD FieldMetaInformation object
            object containing meta-information about the grid.
        """

        rho, info_rho = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='%s' % self.Species_name, m='all')
        den           = np.abs( (-1/e)*rho )         # Plasma density (m^-3)

        return magnitude_conversion_vol(den, "", self.den_unit, reciprocal_units = True), info_rho
    

def get_focusing_field_map(ts, Snapshot, field_unit: str = "Giga"):
    """
    get focusing field of plasma wave

    Parameters
    ----------
    ts : Object
        LpaDiagnostics "time series" object.
    Snapshot : int
        Snapshot from simulation to be analysed.
    field_unit : str, optional
        Order of magnitude of focusing field unit, by default "Giga"

    Returns
    -------
    focusing_field: array
        2D array containing plasma wave focusing field data. Default unit: GV/m
    info_field : openPMD FieldMetaInformation object
        object containing meta-information about the grid.
    """

    Ex0, info_field = ts.get_field( iteration=ts.iterations[Snapshot], field='E', m=0, coord='x')                     # Extract Ex field of plasma wave
    By, _           = ts.get_field( iteration=ts.iterations[Snapshot], field='B', m=0, coord='y')                     # Extract By field of plasma wave
    focusing_field = -1*(Ex0-(c*By))

    return magnitude_conversion(focusing_field, "", field_unit), info_field