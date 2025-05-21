"""
This file contains functions relevant for analysing the particles in PIC simualtions.


Date created: 05/01/2024
Authors: Lewis R Reid

TO DO:
    - Take lineouts of histograms (summed and line) should be their own function/class
    - Shift x, z, time axis to centre of beam should be its own function
    - Add type hints to functions int, float, bool, str

"""

import numpy as np
from scipy.constants import c, m_e, e
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion
from PICAnalysisTools.utils import binning, statistics, rounding


class ParticleEnergy():
    # you should put types in here:  eg ux: float
    def __init__(self, ux, uy, uz, q, w, energy_unit: str = "mega", charge_unit: str = "pico"):
        """
        Parameters
        ----------
        ux: float
            Particle momenta in x.
        uy: float
            Particle momenta in y.
        uz: float
            Particle momenta in z.
        q: float
            particle charge (C).
        w: float
            particle weights.
        energy_unit: string
            Order of magnitude of particle energy units
        charge_unit: string
            Order of magnitude of particle charge units
        Ek: float
            Energy of each macroparticle (default: MeV)

        Returns
        -------
        None.
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
            Energy of each beam macroparticles (MeV).
        """
        
        gamma  = np.sqrt(1 + self.ux ** 2 + self.uy ** 2 + self.uz ** 2)        # Get weighted average of the normalised momenta of each macroparticle
        Ek = ((m_e*c**2)/e)*(gamma-1)                                           # Convert normalised momenta to energy in eV 

        return magnitude_conversion(Ek, "", self.energy_unit)
    

    def get_charge(self):
        """
        Returns
        -------
        charge: float
            Total beam charge (C).
        """

        charge = -1*np.sum(self.w * self.q)                               # Calculate charge (C) within selection

        return magnitude_conversion(charge, "", self.charge_unit)
    

    def get_energy_spectrum(self, Spec_Res, E_Round):
        """
        Parameters
        ----------
        Spec_Res: float
            Resolution of histogram.
        E_Round: float
            Rounding value for max/min energy.

        Returns
        -------
        e_spec: tuple
            Particle species energy spectra [0] - number of electrons. [1] - particle energy.
        Spec_Min: int
            Minimumum energy of histogram.
        Spec_Max: int
            Maximum energy of histogram.
        E_Bins: TYPE
            Energy bins of histogram.
        """

        ## TO THINK: 
        """
        If the energy spec is NOT a partical energy property - use a separate class eg:
        class energy_spec:
            def __init__(self, the thimgs you need):
                self.e_spec, self.Spec_Min, self.Spec_Max, self.E_Bins = self._some_function(Ek, Spec_Res, E_Round)

            def _some+function(self):
            ......
        
            in your main code somewhere
         you have this PE class

         spec = energy_spec(PE.Ek, spec_res, e_round, PE.w)

         now spec has spec.SpecMin 

        """

        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)
        e_spec                     = np.histogram(self.Ek, bins=E_Bins, weights=self.w)            # Spectrum with resolution Spec_Res

        return e_spec, Spec_Min, Spec_Max, E_Bins
    
    def beam_energy_properties(self):
        """
        Returns
        -------
        Ek_mean: float
            Average particle energy (MeV).
        de_rms: float
            RMS energy spread from statistical calculation (MeV).
        PC_spread: float
            Percentage RMS energy spread.
        """

        Ek_mean, de_rms = statistics.w_std(self.Ek, weights=self.w )        # Get mean energy and RMS energy spread
        PC_spread       = (de_rms/Ek_mean)*100                              # calculate percentage energy spread

        return Ek_mean, de_rms, PC_spread


class ParticleTransverseProperties:

    def __init__(self, r, ur, z, uz, w, r_unit: str = "micro", div_unit: str = "milli", emit_unit: str = "micro"):
        """
        Parameters
        ----------
        r: float
            Particle positions (m).
        ur: float
            Particle momenta.
        z: float
            Particle positions in z axis (m).
        uz: float
            Particle momenta in z axis.
        w: float
            particle weights.
        r_unit: string
            order of magnitude to return values.

        Returns
        -------
        None.
        """
        self.r  = r
        self.ur = ur
        self.z  = z
        self.uz = uz
        self.w  = w
        self.r_unit    = r_unit
        self.div_unit  = div_unit
        self.emit_unit = emit_unit
        

    def beam_size(self):
        r_mean        = np.average( self.r, weights=self.w)
        r2_mean       = np.average( ((self.r-r_mean)**2), weights=self.w)
        beam_r        = (np.sqrt( np.abs(r2_mean)))         # beam size (m)

        return magnitude_conversion(beam_r, "", self.r_unit)
    
    def beam_length(self):
        z_mean        = np.average( self.z, weights=self.w)
        z2_mean       = np.average( ((self.z-z_mean)**2), weights=self.w)
        beam_z        = (np.sqrt( np.abs(z2_mean)))         # beam size (m)

        return magnitude_conversion(beam_z, "", self.r_unit)

    def transverse_beam_properties(self):
        r_prime       = np.arctan2(self.ur, self.uz)                                    # Get divergence in r (rad)
        _, div_r      = statistics.w_std( np.arctan2(self.ur, self.uz), self.w )        # Get whole beam divergence in r (rad)

        r_mean        = np.average( self.r, weights=self.w)
        r2_mean       = np.average( ((self.r-r_mean)**2), weights=self.w)
        rPrime2_mean  = np.average( (r_prime**2), weights=self.w)
        rrPrime_mean  = np.average( (self.r*r_prime), weights=self.w)
        rrPrime_mean2 = (np.average( ((self.r-r_mean)*r_prime), weights=self.w))**2
        eta_tr_r      = np.sqrt( (r2_mean*rPrime2_mean)  -  rrPrime_mean2 )
        eta_tr_norm_r = (np.average( self.uz, weights=self.w))  * eta_tr_r
        beam_r        = (np.sqrt( np.abs(r2_mean)))         # beam size (m)
        eta_tr        = np.sqrt( (r2_mean*rPrime2_mean)  -  rrPrime_mean2 )

        Twiss_alpha = -1*(rrPrime_mean/eta_tr)
        Twiss_beta  = r2_mean/eta_tr
        Twiss_gamma = rPrime2_mean/eta_tr

        return magnitude_conversion(div_r, "", self.div_unit), magnitude_conversion(eta_tr_norm_r, "", self.emit_unit), magnitude_conversion(beam_r, "", self.r_unit), Twiss_alpha, Twiss_beta, Twiss_gamma
    

class PhaseSpace(ParticleEnergy):
    
    def __init__(self, x, y, z, ux, uy, uz, w, r_unit: str = "micro", energy_unit: str = "mega", div_unit: str = "milli", time_unit: str = "femto"):
        #super().__init__(self)
        """
        Parameters
        ----------
        x: float
            Particle positions in x (m).
        y: float
            Particle positions in y (m).
        z: float
            Particle positions in z (m).
        ux: float
            Particle momenta in x.
        uy: float
            Particle momenta in y.
        uz: float
            Particle momenta in z.
        w: float
            particle weights.
        Ek: float
            Energy of each macroparticle (default: MeV)

        Returns
        -------
        None.
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
            Energy of each beam macroparticles (default: MeV).
        """
        
        gamma  = np.sqrt(1 + self.ux ** 2 + self.uy ** 2 + self.uz ** 2)        # Get weighted average of the normalised momenta of each macroparticle
        Ek = ((m_e*c**2)/e)*(gamma-1)                                           # Convert normalised momenta to energy in eV 

        return magnitude_conversion(Ek, "", self.energy_unit)


    def get_divergence(self, Coord = 'x'):

        if Coord == 'x':
            Div_r = np.arctan2(self.ux, self.uz)                            # Get divergence in x (rad)
        elif Coord == 'y':
            Div_r = np.arctan2(self.uy, self.uz)                            # Get divergence in y (rad)
        else:
            print("Coordinate incorrectly defined. Choose x or y.\nWARNING: x coordinate chosen as default")
            Div_r = np.arctan2(self.ux, self.uz)                            # Get divergence in x (rad)

        return magnitude_conversion(Div_r, "", self.div_unit)

    def Energy_z_space(self, z_res, z_round, Spec_Res, E_Round, Centre_z: bool = True, Find_Lineouts: bool = True, lineout_height=0.2):

        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)

        z_converted = magnitude_conversion(self.z, "", self.r_unit)

        z_min, z_max, Z_Bins = binning.get_bins(z_converted, z_res, z_round)

        E_Phase, _, _ = np.histogram2d(self.Ek, z_converted, bins=(E_Bins, Z_Bins), weights=self.w )
        Order         = np.floor(np.log10(np.max(E_Phase)))

        if Centre_z is True:
            # Set the centre of the beam to r = 0
            # Particularly useful if plotting z axis.
            col_sum  = E_Phase.sum(axis=0)                          # Sum all columns of histogram
            PosZ     = np.argmax(col_sum)                           # Index of maximum value
            z_at_max = Z_Bins[PosZ] 

            z_max = z_max - z_at_max
            z_min = z_min - z_at_max

        else:
            z_at_max = 0

        if Find_Lineouts is True:
            # Lineout plots are: plt.plot(Row_sum, E_Bins[:-1]) and plt.plot(Z_line[:-1], Col_sum)
            # I don't like the return having two different lengths. This can't be "good" python. I do want the option for calculating the lineouts

            Ene_height = (Spec_Max-Spec_Min) * lineout_height       # lineout 20% of plot window in height
            R_height   = (z_max - z_min) * lineout_height

            Row_sum = (R_height * rounding.normalise(np.sum(E_Phase, axis=1))) + (z_min)
            Z_line  = Z_Bins - z_at_max                                                 # Set centre of beam to z = 0
            Col_sum = (Ene_height * rounding.normalise(np.sum(E_Phase, axis=0))) + Spec_Min

            return E_Phase, Order, Spec_Min, Spec_Max, E_Bins, z_min, z_max, Z_Bins, Row_sum, Z_line, Col_sum

        else:
            return E_Phase, Order, Spec_Min, Spec_Max, E_Bins, z_min, z_max, Z_Bins
        

    def Energy_time_space(self, Spec_Res, E_Round, t_res, t_round, Find_Lineouts: bool = True, lineout_height=0.2):

        # Need to add option to make the centre of the beam t = 0. Is it there already and always on as default?

        time_data                  = magnitude_conversion((self.z/c), "", self.time_unit)
        t_min, t_max, T_Bins       = binning.get_bins(time_data, t_res, t_round)
        Spec_Min, Spec_Max, E_Bins = binning.get_bins(self.Ek, Spec_Res, E_Round)

        long_Phase, _, _ = np.histogram2d(self.Ek, time_data, bins=(E_Bins, T_Bins), weights=self.w )
        Order            = np.floor(np.log10(np.max(long_Phase)))

        # Set the centre of the beam to r = 0
        # Particularly useful if plotting z axis.
        col_sum  = long_Phase.sum(axis=0)           # Sum all columns of histogram
        PosZ     = np.argmax(col_sum)               # Index of maximum value
        t_at_max = T_Bins[PosZ] 

        t_max = t_max - t_at_max
        t_min = t_min - t_at_max

        if Find_Lineouts is True:
            # Lineout plots are: plt.plot(Row_sum, E_Bins[:-1]) and plt.plot(Z_line[:-1], Col_sum)

            Ene_height = (Spec_Max-Spec_Min) * lineout_height       # lineout 20% of plot window in height
            R_height   = (t_max - t_min) * lineout_height

            Row_sum = (R_height * rounding.normalise(np.sum(long_Phase, axis=1))) + (t_min)
            T_line  = T_Bins - t_at_max                                                 # Set centre of beam to z = 0
            Col_sum = (Ene_height * rounding.normalise(np.sum(long_Phase, axis=0))) + Spec_Min

            return long_Phase, Order, Spec_Min, Spec_Max, E_Bins, t_min, t_max, T_Bins, Row_sum, T_line, Col_sum
        
        else:
            return long_Phase, Order, Spec_Min, Spec_Max, E_Bins, t_min, t_max, T_Bins


    def Div_r_space(self, Coord, r_res, r_round, div_round, div_res, Find_Lineouts: bool = True, lineout_height=0.2):

        if Coord == "x":
            Div_r = self.div_x
            r     = magnitude_conversion(self.x, "", self.r_unit)
        elif Coord == "y":
            Div_r = self.div_y
            r     = magnitude_conversion(self.y, "", self.r_unit)
        else:
            print("Axis incorrectly defined. Choose x ot y. x chosen as default.")
            Div_r = self.div_x
            r     = magnitude_conversion(self.x, "", self.r_unit)
           
        _, max_div, D_Bins   = binning.get_bins_absolute(Div_r, div_res, div_round)
        R_min, R_max, R_Bins = binning.get_bins_absolute(r, r_res, r_round)

        Trans_phase, _, _ = np.histogram2d(Div_r, r, bins=(D_Bins, R_Bins), weights=self.w )
        Order_Trans       = np.floor(np.log10(np.amax(Trans_phase)))

        if Find_Lineouts is True:
            # Lineout plots are: plt.plot(Row_sum_Trans, D_Bins[:-1]) and plt.plot(R_line, col_sum_Trans)
            # This only makes sense with x and y coords so no need to set centre of beam to zero, unlike Energy_r_space()
            X_height   = (2*max_div) * lineout_height                     # lineout 20% of plot window in height
            Div_height = (R_max - R_min) * lineout_height
    
            Row_sum_Trans = (Div_height * rounding.normalise(np.sum(Trans_phase, axis=1))) - abs(R_min)
            R_line        = R_Bins[:-1]
            col_sum_Trans = (X_height * rounding.normalise(np.sum(Trans_phase, axis=0))) + (-1*max_div)

            return Trans_phase, Order_Trans, max_div, D_Bins, R_min, R_max, R_Bins, Row_sum_Trans, R_line, col_sum_Trans
        
        else:
            return Trans_phase, Order_Trans, max_div, D_Bins, R_min, R_max, R_Bins


class BeamProjection():

    """
    TO DO:
        - create option (or second function) for the binning to cover the whole simulation window, not just around the beam limits.
        - What is the most efficient way to get simulation window limits?
    """

    def __init__(self, a, b, w, r_unit = "micro"):
        self.a = a
        self.b = b
        self.w = w
        self.r_unit = r_unit
        self.a_converted  = magnitude_conversion(self.a, "", self.r_unit)
        self.b_converted  = magnitude_conversion(self.b, "", self.r_unit)

    def beam_projection(self, r_res, r_round, independant_bins=False, equal_bins=False):
        # independant is good for x vs z (or y vs z)
        # global is good for x vs y

        # equal_bins=False best for for x vs z (or y vs z)
        # equal_bins=True for x vs y

        if independant_bins is True:
            
            if equal_bins is True:
                A_min, A_max, A_Bins = binning.get_bins_absolute(self.a_converted, r_res, r_round)
                B_min, B_max, B_Bins = binning.get_bins_absolute(self.b_converted, r_res, r_round)
            else:
                A_min, A_max, A_Bins = binning.get_bins(self.a_converted, r_res, r_round)
                B_min, B_max, B_Bins = binning.get_bins(self.b_converted, r_res, r_round)

            projection, _, _ = np.histogram2d(self.a_converted, self.b_converted, bins=(A_Bins, B_Bins), weights=self.w )

            return projection, [A_min, A_max, B_min, B_max], A_Bins, B_Bins

        else:

            if equal_bins is True:
                r_min, r_max, R_Bins = binning.get_bins_absolute(np.concatenate((self.a_converted, self.b_converted)), r_res, r_round)
            else:
                r_min, r_max, R_Bins = binning.get_bins(np.concatenate((self.a_converted, self.b_converted)), r_res, r_round)

            projection, _, _ = np.histogram2d(self.a_converted, self.b_converted, bins=(R_Bins, R_Bins), weights=self.w )

            return projection, [r_min, r_max, r_min, r_max], R_Bins, R_Bins
        
    def beam_projection_fixed_window(self, window_extent, r_res):
        # Calculate beam projection with user defined histogram limits.
        # Useful for calcuating beam projection with whole simulation box view.

        window_extent_converted  = magnitude_conversion(np.array(window_extent), "", self.r_unit)

        A_Bins = np.arange(window_extent_converted[0], window_extent_converted[1]+r_res, r_res)
        B_Bins = np.arange(window_extent_converted[2], window_extent_converted[3]+r_res, r_res)

        projection, _, _ = np.histogram2d(self.a_converted, self.b_converted, bins=(A_Bins, B_Bins), weights=self.w )

        return projection, window_extent_converted, A_Bins, B_Bins