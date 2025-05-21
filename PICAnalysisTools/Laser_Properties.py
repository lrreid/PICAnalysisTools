"""
This file contains functions for collecting together analysis of laser pulse properties


Date created: 05/01/2024 @ 16:11
Authors: Lewis R Reid

This is only compatible with OpenPMD viewer version 1.9.0 or greater.


TO DO:
    - 


Possible change:
    - Add someting like self.laser_tracking = get_laser_properties_loop(SimPath) to the init
    - Then create two new functions for save_txt_files() and save_plts() which call self.laser_tracking

"""

import numpy as np
from scipy.constants import c, m_e, e, pi
from scipy.signal import hilbert
from os.path import exists, join
from os import makedirs, getcwd
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, get_order_letter
from PICAnalysisTools.utils.statistics import w_std


class LaserProperties():

    def __init__(self, ts, Coord = "x", Method: str = "fit", Profile_Method: str = "projection", centre_unit: str = "micro", spot_unit: str = "micro", time_unit: str = "femto"):
        self.ts             = ts
        self.centre_unit    = centre_unit
        self.spot_unit      = spot_unit
        self.time_unit      = time_unit
        self.Coord          = Coord
        self.Method         = Method
        self.Profile_Method = Profile_Method

    def get_laser_properties_single(self, Snapshot):
        '''
        Function to extract properties of laser properties from simualtion snapshot.
    
        Parameters
        ----------
        ts: Object
            LpaDiagnostics "time series" object.
        Snapshot: int
            Snapshot from simulation to be analysed.
        Coord: string, default 'x'
            Polarization of the field. Options are 'x', 'y'.
        method: str, optional, default 'fit'
           The method which is used to compute the waist
           'fit': Gaussian fit of the transverse profile
           'rms': RMS radius, weighted by the transverse profile
           ('rms' tends to give more weight to the "wings" of the pulse)

        profile_Method: str, optional, default 'projection' ---> note this is different to the openPMD viewer default.
            Method used to obtain the transverse profile. Options are: 'peak','projection'
    
        Returns
        -------
        Waist: float
            laser waist (m).
        t_FWHM: float
            laser pulse duration FWHM (s).
        a0: float
            Normalised laser vector potential.
        Centroid: float
            Position of centre of laser pusle (m)
    
        '''

        Waist        = self.ts.get_laser_waist(iteration = self.ts.iterations[Snapshot], pol=self.Coord, method=self.Method, profile_method=self.Profile_Method) # OpenPMD viewwer defaults are fit and peak. New to version 1.9.0.
        Pulse_Length = self.ts.get_ctau(iteration        = self.ts.iterations[Snapshot], pol=self.Coord, method=self.Method)
        t_FWHM       = np.sqrt(2*np.log(2))*(Pulse_Length/c)                 # Calculate pulse duration (s)
        a0           = self.ts.get_a0(iteration          = self.ts.iterations[Snapshot], pol=self.Coord)

        Centroid, _ = get_laser_cenroid(self.ts, Snapshot, centroid_unit = self.centre_unit)

        return Centroid, a0, magnitude_conversion(Waist, "", self.spot_unit), magnitude_conversion(t_FWHM, "", self.time_unit)
    

    def get_laser_properties_loop(self, save_txt: bool = False, Sim_Path=getcwd()):

        output = np.zeros([len(self.ts.iterations), 4])                                  # Create array of zeros to store analysed data

        for k in range(len(self.ts.iterations)):
            laser_properties = self.get_laser_properties_single(k)

            output[k, 0] = np.round(laser_properties[0], 4)       # Add laser centroid to output matrix
            output[k, 1] = np.round(laser_properties[1], 4)       # Add laser normalsied intensity to output matrix
            output[k, 2] = np.round(laser_properties[2], 4)       # Add laser waist to output matrix
            output[k, 3] = np.round(laser_properties[3], 4)       # Add laser pulse duration to output matrix

        if save_txt is True:
            
            Ana_dir = join(Sim_Path, 'Analysed', 'Laser_properties_tracking')

            if exists(Ana_dir) is False:
                makedirs(Ana_dir)

            title_centroid = "Centroid (%sm)" % get_order_letter(self.centre_unit)
            title_a0       = "a0"
            title_spot     = "w0 (%sm)" % get_order_letter(self.spot_unit)
            title_duration = "tau_FWHM (%ss)" % get_order_letter(self.time_unit)

            titles   = ["%s" % title_centroid, "%s" % title_a0, "%s" % title_spot, "%s" % title_duration]
            np.savetxt("%s/Laser_Tracking_summary.txt" % (Ana_dir), output, fmt=('%.4f'), delimiter = '\t', header  = '\t'.join(titles), comments = '' )

        return output


def get_a0_field_map(ts, Snapshot):

    central_wavelength, _, _ = get_central_wavelength(ts, Snapshot, m='all', coord='x', wavelength_unit = "" )

    Ex, _    = ts.get_field( iteration=ts.iterations[Snapshot], field='E', m=1, coord='x')
    a0_field = np.array((np.absolute(hilbert(np.array(Ex)))*e*central_wavelength)/(m_e*c*c*2*pi))    # Laser a0 field
    a0_max   = np.max(a0_field)

    return a0_field, a0_max


def get_spectrum(ts, snapshot, m: str = 'all', coord: str = 'x', wavelength_unit: str = "nano" ):
    # Adapted from openpmd_viewer function of the same name

    E_field, info_field = ts.get_field( iteration=ts.iterations[snapshot], field='E', m=m, coord=coord)

    # FFT of 1d data
    dt        = info_field.dz/c                                     #  Integration step for the FFT
    fft_field = np.fft.fft(E_field[int(len(info_field.r)//2),:]) * dt
    spectrum  = abs( fft_field[ : int( len(fft_field) / 2 ) ] )     # Take half of the data (positive frequencies only)

    dO         = (2*pi*c)/(info_field.zmax - info_field.zmin)       # Grid spacing of angular frequency 
    Omega_max  = dO*len(spectrum)                                   # Peak angular frequency in spectrum
    Omega      = np.arange(dO,Omega_max,dO)                         # Array of angular frequency points
    wavelength = ((2*pi*c)/Omega)                                   # Convert omega to wavelength (m)

    return magnitude_conversion(wavelength, "", wavelength_unit), spectrum[1:]


def get_central_wavelength(ts, snapshot, m: str = 'all', coord: str = 'x', wavelength_unit: str = "nano" ):

    wavelength, spectrum = get_spectrum(ts=ts, snapshot=snapshot, m=m, coord=coord, wavelength_unit = wavelength_unit)
    average, width      = w_std(wavelength, spectrum)
    peak                = wavelength[np.where(spectrum == np.max(spectrum))[0][0]]

    return magnitude_conversion(average, "", wavelength_unit), magnitude_conversion(width, "", wavelength_unit), magnitude_conversion(peak, "", wavelength_unit)


def get_laser_cenroid(ts, snapshot, centroid_unit: str = "micro"):

    Ex, info_Ex = ts.get_field( iteration=ts.iterations[snapshot], field='E', m=1, coord='x')
    Ex_env      = np.abs(hilbert(Ex[int(len(info_Ex.r)//2),:]))         # Get the laser envelope
    centroid, _ = w_std(info_Ex.z, Ex_env)                              # Find the laser centroid
    Peak        = info_Ex.z[np.where(Ex_env == np.max(Ex_env))[0][0]]

    return magnitude_conversion(centroid, "", centroid_unit), magnitude_conversion(Peak, "", centroid_unit)