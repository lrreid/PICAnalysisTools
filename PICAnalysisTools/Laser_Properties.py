"""
This file contains functions for collecting together analysis of laser pulse properties


Date created: 05/01/2024 @ 16:11
Authors: Lewis R Reid

This is only compatible with OpenPMD viewer version 1.9.0 or greater.

TO DO:
    - Allow for different input units in a0_from_intensity, intensity_from_a0
"""

import numpy as np
from scipy.constants import c
from os.path import exists
from os import makedirs, getcwd
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion, get_order_letter


class LaserProperties():

    def __init__(self, ts, centre_unit = "micro", spot_unit = "micro", time_unit = "femto"):
        self.ts          = ts
        self.centre_unit = centre_unit
        self.spot_unit   = spot_unit
        self.time_unit   = time_unit


    def get_laser_properties_single(self, Snapshot, Coord='x', Method='fit', Profile_Method='projection'):
        '''
        Function to extract properties of laser properties from simualtion snapshot.
    
        Parameters
        ----------
        ts : Object
            LpaDiagnostics "time series" object.
        Snapshot : int
            Snapshot from simulation to be analysed.
        Coord : string, default 'x'
            Polarization of the field. Options are 'x', 'y'.
        method : str, optional, default 'fit'
           The method which is used to compute the waist
           'fit': Gaussian fit of the transverse profile
           'rms': RMS radius, weighted by the transverse profile
           ('rms' tends to give more weight to the "wings" of the pulse)

        profile_Method : str, optional, default 'projection' ---> note this is different to the openPMD viewer default.
            Method used to obtain the transverse profile. Options are: 'peak','projection'
    
        Returns
        -------
        Waist : float
            laser waist (m).
        t_FWHM : float
            laser pulse duration FWHM (s).
        a0 : float
            Normalised laser vector potential.
        Centroid : float
            Position of centre of laser pusle (m)
    
        '''

        Waist        = self.ts.get_laser_waist(iteration = self.ts.iterations[Snapshot], pol=Coord, method=Method, profile_method=Profile_Method) # OpenPMD viewwer defaults are fit and peak. New to version 1.9.0.
        Pulse_Length = self.ts.get_ctau(iteration        = self.ts.iterations[Snapshot], pol=Coord, method=Method)
        t_FWHM       = np.sqrt(2*np.log(2))*(Pulse_Length/c)                 # Calculate pulse duration (s)
        a0           = self.ts.get_a0(iteration          = self.ts.iterations[Snapshot], pol=Coord)

        Ex, info_Ex = self.ts.get_field( iteration=self.ts.iterations[Snapshot], field='E', m=1, coord='x')
        Ex_lineout  = np.array(Ex[int(np.size(info_Ex.r))//2,:])
        Ex_max_ind  = np.where(abs(Ex_lineout) == np.amax(abs(Ex_lineout)))[0][0]
        Centroid    = info_Ex.z[Ex_max_ind]                                 # This is the max electric field. Change to centroid from statistical calculation? Might work better for complex laser pulses

        return magnitude_conversion(Centroid, "", self.centre_unit), a0, magnitude_conversion(Waist, "", self.spot_unit), magnitude_conversion(t_FWHM, "", self.time_unit)
    

    def get_laser_properties_loop(self, coord='x', method='fit', profile_Method='projection', save_txt=False, Sim_path=getcwd()):

        output = np.zeros([len(self.ts.iterations), 4])                                  # Create array of zeros to store analysed data

        for k in range(len(self.ts.iterations)):
            laser_properties = self.get_laser_properties_single(k, coord, method, profile_Method)

            output[k, 0] = round(laser_properties[0], 4)       # Add laser centroid to output matrix
            output[k, 1] = round(laser_properties, 4)       # Add laser normalsied intensity to output matrix
            output[k, 2] = round(laser_properties[2], 4)       # Add laser waist to output matrix
            output[k, 3] = round(laser_properties[3], 4)       # Add laser pulse duration to output matrix

        if save_txt == True:
            
            Ana_dir = Sim_path + "\\Analysed\\Laser_properties_tracking" # Is there a way to automatically find simpath rather than having path to data as an input??

            if exists(Ana_dir) == False:
                makedirs(Ana_dir)

            title_centroid = "Centroid (%sm)" % get_order_letter(self.centre_unit)
            title_a0       = "a0"
            title_spot     = "w0 (%sm)" % get_order_letter(self.spot_unit)
            title_duration = "tau_FWHM (%ss)" % get_order_letter(self.time_unit)

            titles   = ["%s" % title_centroid, "%s" % title_a0, "%s" % title_spot, "%s" % title_duration] # Need to automate what these units are!
            np.savetxt("%s/Laser_Tracking_summary.txt" % (Ana_dir), output, fmt=('%.4f'), delimiter='\t', header ='\t'.join(titles), comments='' )

        return output