"""
This file contains functions relevant for analysing the particles in PIC simualtions.

TO DO:
    - Take lineouts of histograms (summed and line) should be their own function/class
    - Shift x, z, time axis to centre of beam should be its own function
    - Adapt transverse particle properties to choose order of magnitude of Twiss beta values.

"""

import numpy as np
from scipy.constants import c, m_e, e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion
from PICAnalysisTools.utils import binning, statistics, rounding
from PICAnalysisTools.utils.statistics import w_std
elec_rest_mass = (m_e*c**2)/e                   # Electron rest mass (eV)

class ParticleEnergy():

    def __init__(self, ux, uy, uz, q, w, energy_unit: str = "mega", charge_unit: str = "pico"):
        """
        Parameters
        ----------
        ux : array_like
            x-coordinates of particle momenta (px/mc)
        uy : array_like
            y-coordinates of particle momenta (py/mc)
        uz : array_like
            z-coordinates of particle momenta (pz/mc)
        q: float
            particle charge (C).
        w : array_like
            Macroparticle weights
        energy_unit: string
            Order of magnitude of particle energy units, by default "mega"
        charge_unit: string
            Order of magnitude of particle charge units, by default "pico"
        """
        
        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.q  = q
        self.w  = w
        self.energy_unit = energy_unit
        self.charge_unit = charge_unit
        self.Ek          = self.get_particle_energy()
        self.charge      = self.get_charge()
        self.Ek_mean, self.de_rms, self.PC_spread = self.beam_energy_properties() 

        
    def get_particle_energy(self):
        """
        Returns
        -------
        Ek_converted: float
            Energy of each beam macroparticles. Default unit: MeV.
        """
        
        gamma  = np.sqrt(1 + self.ux ** 2 + self.uy ** 2 + self.uz ** 2)        # Get weighted average of the normalised momenta of each macroparticle
        Ek = ((m_e*c**2)/e)*(gamma-1)                                           # Convert normalised momenta to energy in eV 

        return magnitude_conversion(Ek, "", self.energy_unit)
    

    def get_charge(self):
        """
        Returns
        -------
        charge: float
            Total beam charge. Default unit: pC.
        """

        charge = -1*np.sum(self.w * self.q)                               # Calculate charge (C) within selection

        return magnitude_conversion(charge, "", self.charge_unit)
    

    def get_energy_spectrum(self, Spec_Res, E_Round):
        """
        Get the energy spectrum of the particle species

        Parameters
        ----------
        Spec_Res : float
            Resolution of histogram. Units in eV with order of magnitude energy_unit.
        E_Round : float
            Rounding value to determine max & min energy histogram bin.

        Returns
        -------
        energy_spec: tuple
            [0] - Energy bins of the histogram. [1] - Number of electrons in each energy bin
        energy_limits: list
            List containing the maximum and minimum particle energy bins. 
        """

        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)
        energy_hist                = np.histogram(self.Ek, bins=E_Bins, weights=self.w)            # Spectrum with resolution Spec_Res
        energy_spec                = (E_Bins, (np.append(energy_hist[0], 0)))                      # Tuple with histogram bins and and number of particles per bin. Extra zero to ensure arrays are the same length

        return energy_spec, [Spec_Min, Spec_Max]
    
    def beam_energy_properties(self):
        """
        Get statistical properties of particle species

        Returns
        -------
        Ek_mean: float
            Average particle energy. Default unit: MeV
        de_rms: float
            RMS energy spread from statistical calculation. Default unit: MeV
        PC_spread: float
            Percentage RMS energy spread.
        """

        Ek_mean, de_rms = statistics.w_std(self.Ek, weights=self.w )        # Get mean energy and RMS energy spread
        PC_spread       = (de_rms/Ek_mean)*100                              # calculate percentage energy spread

        return Ek_mean, de_rms, PC_spread


class ParticleTransverseProperties:

    def __init__(self, r, ur, z, uz, w, r_unit: str = "micro", div_unit: str = "milli", emit_unit: str = "micro", beta_unit: str = "milli", gamma_unit: str = ""):
        """
        Class to calculate transverse properties of a particle species

        Parameters
        ----------
        r : array_like
            Particle positions. Could be x, y or r. (m).
        ur : array_like
            Particle momenta. Could be ux, uy or ur to match positions. (pr/mc)
        z : array_like
            z-coordinates of particle positions (m)
        uz : array_like
            z-coordinates of particle momenta (pz/mc)
        w : array_like
            Macroparticle weights
        r_unit : str, optional
            Order of magnitude of beam size unit, by default "micro"
        div_unit : str, optional
            Order of magnitude of divergence unit, by default "milli"
        emit_unit : str, optional
            Order of magnitude of emittance unit, by default "micro"
        """

        self.r  = r
        self.ur = ur
        self.z  = z
        self.uz = uz
        self.w  = w
        self.r_unit     = r_unit
        self.div_unit   = div_unit
        self.emit_unit  = emit_unit
        self.beta_unit  = beta_unit
        self.gamma_unit = gamma_unit
        self.beam_r, self.div_r, self.eta_tr_norm_r, self.Twiss_alpha, self.Twiss_beta, self.Twiss_gamma = self.transverse_beam_properties()
        self.beam_z    = self.beam_length()
    
    def beam_length(self):
        """
        Calculate the length of the particle beam

        Returns
        -------
        beam_z: float
            RMS longitudinal length of particle beam. Default unit: microns
        """
        z_mean        = np.average( self.z, weights=self.w)
        z2_mean       = np.average( ((self.z-z_mean)**2), weights=self.w)
        beam_z        = (np.sqrt( np.abs(z2_mean)))         # beam size (m)

        return magnitude_conversion(beam_z, "", self.r_unit)

    def transverse_beam_properties(self):
        """
        Calculate transverse properties of particle species

        Returns
        -------
        beam_r: float
            RMS transverse length of particle beam. Default unit: microns
        div_r: float
            Divergence of particle beam. Default unit: mrad
        eta_tr_norm_r: float
            Normalised emittance of particle beam. mm mrad
        Twiss_alpha: float
            Twiss alpha parameter of particle beam. Dimensionless
        Twiss_beta: float
            Twiss beta parameter of particle beam. Unit: m
        Twiss_gamma: float
            Twiss gamma parameter of particle beam. Unit: 1/m
        """
        r_prime       = np.arctan2(self.ur, self.uz)                                    # Get divergence in r (rad)
        _, div_r      = statistics.w_std( np.arctan2(self.ur, self.uz), self.w )        # Get whole beam divergence in r (rad)

        r_mean        = np.average( self.r, weights=self.w)
        r2_mean       = np.average( ((self.r-r_mean)**2), weights=self.w)
        beam_r        = (np.sqrt( np.abs(r2_mean)))         # beam size (m)

        rPrime2_mean  = np.average( (r_prime**2), weights=self.w)
        rrPrime_mean  = np.average( (self.r*r_prime), weights=self.w)
        rrPrime_mean2 = (np.average( ((self.r-r_mean)*r_prime), weights=self.w))**2
        eta_tr_r      = np.sqrt( (r2_mean*rPrime2_mean)  -  rrPrime_mean2 )
        eta_tr_norm_r = (np.average( self.uz, weights=self.w))  * eta_tr_r
        eta_tr        = np.sqrt( (r2_mean*rPrime2_mean)  -  rrPrime_mean2 )

        Twiss_alpha = -1*(rrPrime_mean/eta_tr)
        Twiss_beta  = r2_mean/eta_tr
        Twiss_gamma = rPrime2_mean/eta_tr

        return magnitude_conversion(beam_r, "", self.r_unit), magnitude_conversion(div_r, "", self.div_unit), magnitude_conversion(eta_tr_norm_r, "", self.emit_unit), Twiss_alpha, magnitude_conversion(Twiss_beta, "", self.beta_unit), magnitude_conversion(Twiss_gamma, "", self.gamma_unit, reciprocal_units = True)


class PhaseSpace():
    
    def __init__(self, x, y, z, ux, uy, uz, w, r_unit: str = "micro", energy_unit: str = "mega", div_unit: str = "milli", time_unit: str = "femto"):
        """
        Calculate various phase space distributions of a particle beam species.

        Parameters
        ----------
        x : array_like
            x-coordinates of particle positions (m)
        y : array_like
            y-coordinates of particle positions (m)
        z : array_like
            z-coordinates of particle positions (m)
        ux : array_like
            x-coordinates of particle momenta (px/mc)
        uy : array_like
            y-coordinates of particle momenta (py/mc)
        uz : array_like
            z-coordinates of particle momenta (pz/mc)
        w : array_like
            Macroparticle weights
        r_unit : str, optional
            Order of magnitude of beam size unit, by default "micro"
        energy_unit: string
            Order of magnitude of particle energy units, by default "mega"
        div_unit : str, optional
            Order of magnitude of divergence unit, by default "milli"
        time_unit : str, optional
            Order of magnitude of time unit, by default "femto"
        """

        self.x  = x
        self.y  = y
        self.z  = z
        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.w  = w
        self.r_unit      = r_unit
        self.energy_unit = energy_unit
        self.div_unit    = div_unit
        self.time_unit   = time_unit
        self.Ek          = self.get_particle_energy()
        self.div_x       = self.get_divergence(Coord = 'x')
        self.div_y       = self.get_divergence(Coord = 'y') # This won't work

    def get_particle_energy(self):
        """
        Returns
        -------
        Ek_converted: float
            Energy of each beam macroparticles. Default unit: MeV.
        """
        
        gamma  = np.sqrt(1 + self.ux ** 2 + self.uy ** 2 + self.uz ** 2)        # Get weighted average of the normalised momenta of each macroparticle
        Ek = ((m_e*c**2)/e)*(gamma-1)                                           # Convert normalised momenta to energy in eV 

        return magnitude_conversion(Ek, "", self.energy_unit)


    def get_divergence(self, Coord: str = 'x'):
        """
        Get particle species divergence in a given axis. Calculates divergence of each individial macroparticle, not the mean of all particles.

        Parameters
        ----------
        Coord : str, optional
            Coordinate of the particle species that the divergence will be calculated for, by default 'x'

        Returns
        -------
        Div_r: array_like
            Array of divegence of each beam macroparticle. Default unit: mrad.
        """

        if Coord == 'x':
            Div_r = np.arctan2(self.ux, self.uz)                            # Get divergence in x (rad)
        elif Coord == 'y':
            Div_r = np.arctan2(self.uy, self.uz)                            # Get divergence in y (rad)
        else:
            print("Coordinate incorrectly defined. Choose x or y.\nWARNING: x coordinate chosen as default")
            Div_r = np.arctan2(self.ux, self.uz)                            # Get divergence in x (rad)

        return magnitude_conversion(Div_r, "", self.div_unit)

    def Energy_z_space(self, z_res, z_round, Spec_Res, E_Round, Centre_z: bool = True, lineout_height: float = 0.2):
        """
        Calculate the longitudinal phase / energy space with lineouts.

        Parameters
        ----------
        z_res : float
            Restolution of histogram binning in z-axis. Default unit: microns
        z_round : float
            Rounding value for z-axis to determine max & min z histogram bin. Default unit: microns
        Spec_Res : float
            Restolution of histogram binning in particle energy axis. Default unit: MeV.
        E_Round : _type_
            Rounding value to determine max & min energy histogram bin. Default unit: MeV.
        Centre_z : bool, optional
            Option to centre the beam on z = 0, by default True
        lineout_height : float, optional
            Height of lineout relative to the plot axes limits, by default 0.2

        Returns
        -------
        Ez_Phase: array_like
            The bi-dimensional histogram of samples z and energy. Values in z are histogrammed along the first dimension and values in energy are histogrammed along the second dimension.
        plt_limits: list
            List containing the maximum and minimum bin values in each axis. Useful for setting plot axes limits.
            [z_min, z_max, Spec_Min, Spec_Max]
        Z_Bins
            Array containing histogram bins of z axis
        E_Bins: array_like
            Array containing histogram bins of energy axis
        z_line: array_like
            Lineout of phase space data in z axis
        Energy_line: array_like
            Lineout of phase space data in energy axis
        """

        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)                                               # Get histogram bins for energy axis
        z_converted                = magnitude_conversion(self.z, "", self.r_unit)                                              # Convert z-coordinates of macroparticles from meters to chosen order of magnitude
        z_min, z_max, Z_Bins       = binning.get_bins(z_converted, z_res, z_round)                                              # Get histogram bins for z axis
        Ez_Phase                   = np.flipud(np.histogram2d(self.Ek, z_converted, bins=(E_Bins, Z_Bins), weights=self.w )[0]) # flipped so min in y-axis is on the bottom of the plot in imshow, this is opposite to default imshow.

        if Centre_z is True:
            # Set the centre of the beam to z = 0
            z_centroid = w_std(Z_Bins, np.append(Ez_Phase.sum(axis=0), 0))[0]    # Sum all columns of histogram and find statistical centroid of the beam
            z_at_max   = rounding.round_nearest(z_centroid, z_res)              # Round centroid to the nearest histogram bin
            z_max      = z_max - z_at_max
            z_min      = z_min - z_at_max
            Z_Bins     = Z_Bins - z_at_max                                      # Set centre of beam to z = 0
            
        plt_limits = [z_min, z_max, Spec_Min, Spec_Max]

        # Calculate lineouts for z and energy axes. Height of the line normalised to the plot limits
        Ene_height = (Spec_Max-Spec_Min) * lineout_height       # lineout 20% of plot window in height by default
        R_height   = (z_max - z_min) * lineout_height           # lineout 20% of plot window in height by default

        Energy_line = np.flip(np.append( (R_height * rounding.normalise(np.sum(Ez_Phase, axis=1))) + z_min , 0))
        z_line      = np.append( (Ene_height * rounding.normalise(np.sum(Ez_Phase, axis=0))) + Spec_Min , 0)

        return Ez_Phase, plt_limits, Z_Bins, E_Bins, z_line, Energy_line
        

    def Energy_time_space(self, t_res, t_round, Spec_Res, E_Round, Centre_t: bool = True, lineout_height: float = 0.2):
        """
        Calculate the longitudinal phase (in time) / energy space with lineouts.
        
        TO DO:
        This is not very robust against the location of t = 0.
        Should I add an option for a z position offset?


        Parameters
        ----------
        t_res : float
            Restolution of histogram binning in time-axis. Default unit: fs
        t_round : _type_
            Rounding value for time-axis to determine max & min t histogram bin. Default unit: fs
        Spec_Res : float
            Restolution of histogram binning in particle energy axis. Default unit: MeV.
        E_Round : _type_
            Rounding value to determine max & min energy histogram bin. Default unit: MeV.
        Centre_z : bool, optional
            Option to centre the beam on z = 0, by default True
        lineout_height : float, optional
            Height of lineout relative to the plot axes limits, by default 0.2

        Returns
        -------
        Et_Phase: array_like
            The bi-dimensional histogram of samples time and energy. Values in t are histogrammed along the first dimension and values in energy are histogrammed along the second dimension.
        plt_limits: list
            List containing the maximum and minimum bin values in each axis. Useful for setting plot axes limits.
            [t_min, t_max, Spec_Min, Spec_Max]
        T_Bins
            Array containing histogram bins of time axis
        E_Bins: array_like
            Array containing histogram bins of energy axis
        t_line: array_like
            Lineout of phase space data in time axis
        Energy_line: array_like
            Lineout of phase space data in energy axis
        """

        time_data                  = magnitude_conversion((self.z/c), "", self.time_unit)
        t_min, t_max, T_Bins       = binning.get_bins(time_data, t_res, t_round)
        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)
        Et_Phase                   = np.flipud(np.histogram2d(self.Ek, time_data, bins=(E_Bins, T_Bins), weights=self.w )[0]) # flipped so min in y-axis is on the bottom of the plot in imshow, this is opposite to default imshow.

        if Centre_t is True:
            # Set the centre of the beam to t = 0
            z_centroid = w_std(T_Bins, np.append(Et_Phase.sum(axis=0), 0))[0]       # Sum all columns of histogram and find statistical centroid of the beam
            t_at_max   = rounding.round_nearest(z_centroid, t_res)                  # Round centroid to the nearest histogram bin
            t_max      = t_max - t_at_max
            t_min      = t_min - t_at_max
            T_Bins     = T_Bins - t_at_max                                          # Set centre of beam to t = 0

        plt_limits = [t_min, t_max, Spec_Min, Spec_Max]

        # Calculate lineouts for time and energy axes. Height of the line normalised to the plot limits
        Ene_height = (Spec_Max-Spec_Min) * lineout_height                       # lineout 20% of plot window in height by default
        R_height   = (t_max - t_min) * lineout_height                           # lineout 20% of plot window in height by default

        Energy_line = np.flip(np.append( (R_height * rounding.normalise(np.sum(Et_Phase, axis=1))) + t_min , 0))
        t_line      = np.append( (Ene_height * rounding.normalise(np.sum(Et_Phase, axis=0))) + Spec_Min , 0)

        return Et_Phase, plt_limits, T_Bins, E_Bins, t_line, Energy_line
        

    def Div_r_space(self, Coord: str, r_res, r_round, div_res, div_round, lineout_height: float = 0.2):
        """
        Calculate the transverse phase with lineouts.

        Parameters
        ----------
        Coord : str
            Transverse coordinate for calculating phase space
        r_res : float
            Restolution of histogram binning in transverse axis. Default unit: microns
        r_round : float
            Rounding value for transverse axis to determine max & min r histogram bin. Default unit: microns
        div_res : float
            Restolution of histogram binning in divergence axis. Default unit: mrad
        div_round : float
            Rounding value for transverse axis to determine max & min divergence histogram bin. Default unit: mrad
        lineout_height : float, optional
            Height of lineout relative to the plot axes limits, by default 0.2

        Returns
        -------
        Trans_phase: array_like
            The bi-dimensional histogram of samples r and divergence. Values in r are histogrammed along the first dimension and values in divergence are histogrammed along the second dimension.
        plt_limits: list
            List containing the maximum and minimum bin values in each axis. Useful for setting plot axes limits.
            [R_min, R_max, div_min, div_max]
        R_Bins
            Array containing histogram bins of r axis
        D_Bins: array_like
            Array containing histogram bins of divergence axis
        col_sum: array_like
            Lineout of phase space data in r axis
        Row_sum: array_like
            Lineout of phase space data in divergence axis
        """

        if Coord == "x":
            Div_r = self.div_x
            r     = magnitude_conversion(self.x, "", self.r_unit)
        elif Coord == "y":
            Div_r = self.div_y
            r     = magnitude_conversion(self.y, "", self.r_unit)
        else:
            print("Axis incorrectly defined. Choose x or y. x chosen as default.")
            Div_r = self.div_x
            r     = magnitude_conversion(self.x, "", self.r_unit)
           
        div_min, div_max, D_Bins = binning.get_bins_absolute(Div_r, div_res, div_round)
        R_min, R_max, R_Bins     = binning.get_bins_absolute(r, r_res, r_round)
        Trans_phase              = np.flipud(np.histogram2d(Div_r, r, bins=(D_Bins, R_Bins), weights=self.w )[0])

        plt_limits = [R_min, R_max, div_min, div_max]

        # The lineouts should the their own function which takes in the histograms and ouputs the lines.
        X_height   = (div_max - div_min) * lineout_height                       # lineout 20% of plot window in height by default
        Div_height = (R_max - R_min) * lineout_height                           # lineout 20% of plot window in height by default

        Row_sum = np.flip(np.append( (Div_height * rounding.normalise(np.sum(Trans_phase, axis=1))) - abs(R_min) , R_min))
        col_sum = np.append( (X_height * rounding.normalise(np.sum(Trans_phase, axis=0))) + (-1*div_max) , div_min)

        return Trans_phase, plt_limits, R_Bins, D_Bins, col_sum, Row_sum


class BeamProjection():

    def __init__(self, a, b, w, r_unit = "micro"):
        """
        Calculate 2D histograms of beam projections

        To Do:
            - Modify to allow for different order of magnitude in the units for a and b.
                This would be useful for xz plots where x is um and z is mm.

        Parameters
        ----------
        a : array_like
            Particle positions in axis a. Could be x, y, z, r (m)
            This will for the x-axis of the histogram
        b : array_like
            Particle positions in axis b. Could be x, y, z, r (m)
            This will for the y-axis of the histogram
        w : array_like
            Macroparticle weights
        r_unit : str, optional
            Order of magnitude of particle positions when returned, by default "micro"
        """

        self.a = a
        self.b = b
        self.w = w
        self.r_unit = r_unit
        self.a_converted  = magnitude_conversion(self.a, "", self.r_unit)
        self.b_converted  = magnitude_conversion(self.b, "", self.r_unit)

    def beam_projection(self, r_res, r_round, independant_bins: bool = False, equal_bins: bool = False):
        """
        Calculate 2D histograms of beam projection.

        To Do:
            - It would be nice to separate out equal_bins for the two axes.
            - equal_bins isn't a name where it is obvious what is does. It should be changed.

        Parameters
        ----------
        r_res : float
            Resolution of histogram bins. Default unit: microns
        r_round : float
            Rounding value for axes to determine max & min histogram bins in each axis. Default unit: microns
        independant_bins : bool, optional
            If true then the historam limits and bin values in the a and b axes are allowed to different.
            If false then the historam limits and bin values are the same for the a and b axes.
            Independant is good for x vs z (or y vs z) and global is good for x vs y.
            Default: False
        equal_bins : bool, optional
            If true then bin_min = -1*bin_max for both the a and b axes
            If false then bin_min and bin_max are automatically found using the data and rounding values
            Default: False

        Returns
        -------
        projection: array_like
            The bi-dimensional histogram of samples a and b. Values in a are histogrammed along the first dimension and values in b are histogrammed along the second dimension.
        plt_limits: array
            Max and min histogram bins for both axes, can be used for plot axes limits. [A_min, A_max, B_min, B_max].
        A_Bins: array_like
            Histogram bins in a-axis
        B_Bins: array_like
            Histogram bins in b-axis
        """

        if independant_bins is True:
            
            if equal_bins is True:
                A_min, A_max, A_Bins = binning.get_bins_absolute(self.a_converted, r_res, r_round)
                B_min, B_max, B_Bins = binning.get_bins_absolute(self.b_converted, r_res, r_round)
            else:
                A_min, A_max, A_Bins = binning.get_bins(self.a_converted, r_res, r_round)
                B_min, B_max, B_Bins = binning.get_bins(self.b_converted, r_res, r_round)

            projection = np.flipud(np.histogram2d(self.b_converted, self.a_converted, bins=(B_Bins, A_Bins), weights=self.w )[0])

            return projection, [A_min, A_max, B_min, B_max], A_Bins, B_Bins

        else:

            if equal_bins is True:
                r_min, r_max, R_Bins = binning.get_bins_absolute(np.concatenate((self.a_converted, self.b_converted)), r_res, r_round)
            else:
                r_min, r_max, R_Bins = binning.get_bins(np.concatenate((self.a_converted, self.b_converted)), r_res, r_round)

            projection = np.flipud(np.histogram2d(self.b_converted, self.a_converted, bins=(R_Bins, R_Bins), weights=self.w )[0])

            return projection, [r_min, r_max, r_min, r_max], R_Bins, R_Bins
        
    def beam_projection_fixed_window(self, window_extent, r_res):
        """
        Calculate beam projection with user defined histogram limits.
        Useful for calcuating beam projection with whole simulation box view.

        Parameters
        ----------
        window_extent : array
            Array containing limits of histogram bins. [A_min, A_max, B_min, B_max].
        r_res : float
            Resolution of histogram bins. Default unit: microns

        Returns
        -------
        projection: array_like
            The bi-dimensional histogram of samples a and b. Values in a are histogrammed along the first dimension and values in b are histogrammed along the second dimension.
        A_Bins: array_like
            Histogram bins in a-axis
        B_Bins: array_like
            Histogram bins in b-axis
        """

        A_Bins = np.arange(window_extent[0], window_extent[1]+r_res, r_res)
        B_Bins = np.arange(window_extent[2], window_extent[3]+r_res, r_res)

        projection = np.flipud(np.histogram2d(self.b_converted, self.a_converted, bins=(B_Bins, A_Bins), weights=self.w )[0])

        return projection, A_Bins, B_Bins
    

def get_normalised_momentum(energy, energy_unit: str = "Mega"):
    """
    Calculate the normalised momenta for a given electron energy

    Parameters
    ----------
    energy : float
        Electron energy. Default unit: MeV
    energy_unit : str, optional
        Order of magnitude of electron energy unit, by default "Mega"

    Returns
    -------
    energy_norm: float
        normalised electron energy
    """

    energy_eV      = magnitude_conversion(energy, energy_unit, "")
    energy_norm    = energy_eV/elec_rest_mass    # Threshold for including electrons in particle diagnostic (normalised units)

    return energy_norm