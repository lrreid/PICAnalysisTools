"""
This file contains functions for collecting together analysis of laser pulse properties

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
from PICAnalysisTools.utils.plot_limits import plt_limits

import matplotlib.pyplot as plt
from matplotlib import use
use("Agg")  

class LaserProperties():

    def __init__(self, ts, Coord = "x", Method: str = "fit", Profile_Method: str = "projection", centroid_unit: str = "micro", spot_unit: str = "micro", time_unit: str = "femto", Sim_Path: str = getcwd() ):
        """
        Class containing functions to extract properties of the laser field from PIC simulations

        Parameters
        ----------
        ts : Object
            LpaDiagnostics "time series" object.
        Coord: string, default 'x'
            Polarization of the field. Options are 'x', 'y'.
        Method : str, optional
           The method which is used to compute the waist
           'fit': Gaussian fit of the transverse profile
           'rms': RMS radius, weighted by the transverse profile
           ('rms' tends to give more weight to the "wings" of the pulse)
           by default "fit"
        Profile_Method : str, optional
            Method used to obtain the transverse profile. Options are: 'peak','projection', by default "projection"
            note this is different to the openPMD viewer default.
        centroid_unit : str, optional
            Order of magnitude of units for laser centroid, by default "micro"
        spot_unit : str, optional
            Order of magnitude of the laser spot size units, by default "micro"
        time_unit : str, optional
            Order of magnitude of pulse duration units, by default "femto"
        Sim_Path : str, optional
            Full path to simulation directory, by default getcwd()
        """

        self.ts             = ts
        self.centroid_unit  = centroid_unit
        self.spot_unit      = spot_unit
        self.time_unit      = time_unit
        self.Coord          = Coord
        self.Method         = Method
        self.Profile_Method = Profile_Method
        self.Sim_Path       = Sim_Path
        self.Ana_dir        = join(Sim_Path, 'Analysed', 'Laser_properties_tracking')

    def get_laser_properties_single(self, Snapshot):
        '''
        Function to extract properties of laser properties from simualtion snapshot.
    
        Parameters
        ----------
        Snapshot: int
            Snapshot from simulation to be analysed.

    
        Returns
        -------
        Centroid: float
            Position of laser centroid. Default unit: microns
        Waist: float
            laser waist. Default unit: microns
        t_FWHM: float
            laser pulse duration FWHM. Default unit: fs
        a0: float
            Normalised laser vector potential.    
        '''

        Waist        = self.ts.get_laser_waist(iteration = self.ts.iterations[Snapshot], pol=self.Coord, method=self.Method, profile_method=self.Profile_Method) # OpenPMD viewwer defaults are fit and peak. New to version 1.9.0.
        Pulse_Length = self.ts.get_ctau(iteration        = self.ts.iterations[Snapshot], pol=self.Coord, method=self.Method)
        t_FWHM       = np.sqrt(2*np.log(2))*(Pulse_Length/c)                 # Calculate pulse duration (s)
        a0           = self.ts.get_a0(iteration          = self.ts.iterations[Snapshot], pol=self.Coord)

        Centroid, _ = get_laser_cenroid(self.ts, Snapshot, centroid_unit = self.centroid_unit)

        return Centroid, a0, magnitude_conversion(Waist, "", self.spot_unit), magnitude_conversion(t_FWHM, "", self.time_unit)
    

    def get_laser_properties_loop(self, save_txt: bool = False):
        """
        Loop through all available h5 files and extract laser properties

        Parameters
        ----------
        save_txt : bool, optional
            If true, text file containing will be saved, by default False

        Returns
        -------
        output: array_like
            Array containing analysed data.
            Col 0 - centroid
            Col 1 - a0
            Col 2 - waist
            Col 3 - FWHM pulse duration
        """

        output = np.zeros([len(self.ts.iterations), 4])                                  # Create array of zeros to store analysed data
        print("\nAnalysing laser properties\nIteration:")

        for k in range(len(self.ts.iterations)):
            print("%d / %d" % (k+1, len(self.ts.iterations)) )

            laser_properties = self.get_laser_properties_single(k)

            output[k, 0] = np.round(laser_properties[0], 4)       # Add laser centroid to output matrix
            output[k, 1] = np.round(laser_properties[1], 4)       # Add laser normalsied intensity to output matrix
            output[k, 2] = np.round(laser_properties[2], 4)       # Add laser waist to output matrix
            output[k, 3] = np.round(laser_properties[3], 4)       # Add laser pulse duration to output matrix

        if save_txt is True:

            if exists(self.Ana_dir) is False:
                makedirs(self.Ana_dir)

            title_centroid = "Centroid (%sm)" % get_order_letter(self.centroid_unit)
            title_a0       = "a0"
            title_spot     = "w0 (%sm)" % get_order_letter(self.spot_unit)
            title_duration = "tau_FWHM (%ss)" % get_order_letter(self.time_unit)

            titles   = ["%s" % title_centroid, "%s" % title_a0, "%s" % title_spot, "%s" % title_duration]
            np.savetxt("%s/Laser_Tracking_summary.txt" % (self.Ana_dir), output, fmt=('%.4f'), delimiter = '\t', header  = '\t'.join(titles), comments = '' )

        print("\n")

        return output
    
    def plot_laser_properties(self, output, centroid_round, a0_round, w0_round, time_round):
        fsize = 12

        if exists(self.Ana_dir) is False:
            makedirs(self.Ana_dir)

        centroid_min, centroid_max = plt_limits(output[:, 0], centroid_round)
        a0_min, a0_max             = plt_limits(output[:, 1], a0_round)
        w0_min, w0_max             = plt_limits(output[:, 2], w0_round)
        tau_min, tau_max           = plt_limits(output[:, 3], time_round)

        fig, ax = plt.subplots(dpi=300)
        plt.plot(output[:,0], output[:,1],'-o', color='dodgerblue', linewidth=2)
        plt.axis([centroid_min, centroid_max, a0_min, a0_max])         # set axes [xmin, xmax, ymin, ymax]
        ax.set_title("laser $a_{0}$ tracking", fontsize=fsize)
        ax.set_xlabel('z (%sm)' % (get_order_letter(self.centroid_unit, return_mu=True)) , fontsize=fsize)
        ax.set_ylabel('$a_{0}$', fontsize=fsize)
        ax.grid()
        plt.savefig("%s/a0_tracking.png" % (self.Ana_dir), dpi=300, format="png")
        plt.close()

        fig, ax = plt.subplots(dpi=300)
        plt.plot(output[:,0], output[:,2],'-o', color='firebrick', linewidth=2)
        plt.axis([centroid_min, centroid_max, w0_min, w0_max])         # set axes [xmin, xmax, ymin, ymax]
        ax.set_title("laser spot size tracking", fontsize=fsize)
        ax.set_xlabel('z (%sm)' % (get_order_letter(self.spot_unit, return_mu=True)) , fontsize=fsize)
        ax.set_ylabel('$\\omega$ (%sm)' % (get_order_letter(self.centroid_unit, return_mu=True)), fontsize=fsize)
        ax.grid()
        plt.savefig("%s/w0_tracking.png" % (self.Ana_dir), dpi=300, format="png")
        plt.close()

        fig, ax = plt.subplots(dpi=300)
        plt.plot(output[:,0], output[:,3],'-o', color='seagreen', linewidth=2)
        plt.axis([centroid_min, centroid_max, tau_min, tau_max])         # set axes [xmin, xmax, ymin, ymax]
        ax.set_title("laser pulse duration tracking", fontsize=fsize)
        ax.set_xlabel('z (%sm)' % (get_order_letter(self.centroid_unit, return_mu=True)) , fontsize=fsize)
        ax.set_ylabel('$\\tau_{FWHM}$ (%ss)' % (get_order_letter(self.time_unit)), fontsize=fsize)
        ax.grid()
        plt.savefig("%s/tau_tracking.png" % (self.Ana_dir), dpi=300, format="png")
        plt.close()


        return None


def get_a0_field_map(ts, Snapshot):
    """
    Get field map for laser a0

    Parameters
    ----------
    ts : Object
        LpaDiagnostics "time series" object.
    Snapshot: int
        Snapshot from simulation to be analysed.

    Returns
    -------
    a0_field: array_like
        2D array containing a0 field
    a0_max: float
        Maximum a0 value in snapshot
    """

    central_wavelength, _, _ = get_central_wavelength(ts, Snapshot, m='all', coord='x', wavelength_unit = "" )

    Ex, _    = ts.get_field( iteration=ts.iterations[Snapshot], field='E', m=1, coord='x')
    a0_field = np.array((np.absolute(hilbert(np.array(Ex)))*e*central_wavelength)/(m_e*c*c*2*pi))    # Laser a0 field
    a0_max   = np.max(a0_field)

    return a0_field, a0_max


def get_spectrum(ts, snapshot, m: str = 'all', coord: str = 'x', wavelength_unit: str = "nano" ):
    """
    Calculate the laser spectrum from the electric field
    Adapted from openpmd_viewer function of the same name

    Parameters
    ----------
    ts : Object
        LpaDiagnostics "time series" object.
    Snapshot: int
        Snapshot from simulation to be analysed.
    m : str, optional
        Azimuthal mode of 2D radial grid, by default 'all'
        Either 'all' (for the sum of all the modes) or an integer (for the selection of a particular mode)
    coord : str, optional
        Coordinate of electric field to analyse, by default 'x'
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"

    Returns
    -------
    wavelength: array
        Spectrum data wavelength axis
    intensity: array
        Spectrum data in intensitiy axis
    """

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
    """
    Calculate the central wavelength and width of laser spectrum

    Parameters
    ----------
    ts : Object
        LpaDiagnostics "time series" object.
    Snapshot: int
        Snapshot from simulation to be analysed.
    m : str, optional
        Azimuthal mode of 2D radial grid, by default 'all'
        Either 'all' (for the sum of all the modes) or an integer (for the selection of a particular mode)
    coord : str, optional
        Coordinate of electric field to analyse, by default 'x'
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"

    Returns
    -------
    centroid: float
        centroid wavelength of spectrum. Default unit: nm
    width: float
        RMS width of spectrum. Default unit: nm
    peak: float
        Wavelength at peak intensity. Default unit: nm
    """

    wavelength, spectrum = get_spectrum(ts=ts, snapshot=snapshot, m=m, coord=coord, wavelength_unit = wavelength_unit)
    centroid, width      = w_std(wavelength, spectrum)
    peak                 = wavelength[np.where(spectrum == np.max(spectrum))[0][0]]

    return centroid, width, peak


def get_laser_cenroid(ts, snapshot, centroid_unit: str = "micro"):
    """
    Get position of laser centroid in the longitudinal axis

    Parameters
    ----------
    ts : Object
        LpaDiagnostics "time series" object.
    Snapshot: int
        Snapshot from simulation to be analysed.
    centroid_unit : str, optional
        Order of magnitude of units for laser centroid, by default "micro"

    Returns
    -------
    centroid: float
        Longitudinal position of laser centroid. Default unit: microns
    Peak: float
        Longitudinal position of laser peak intensity. Default unit: microns
    """

    Ex, info_Ex = ts.get_field( iteration=ts.iterations[snapshot], field='E', m=1, coord='x')
    Ex_env      = np.abs(hilbert(Ex[int(len(info_Ex.r)//2),:]))         # Get the laser envelope
    centroid, _ = w_std(info_Ex.z, Ex_env)                              # Find the laser centroid
    Peak        = info_Ex.z[np.where(Ex_env == np.max(Ex_env))[0][0]]

    return magnitude_conversion(centroid, "", centroid_unit), magnitude_conversion(Peak, "", centroid_unit)