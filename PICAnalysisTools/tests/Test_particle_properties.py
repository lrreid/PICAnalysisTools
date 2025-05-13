# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:25:55 2024

@author: ryi76833

Test object orientated version of PIC viewer

"""


##################################################################################################################
##################################################################################################################
####################### SCRIPT NOT COMPATIBLE WITH CHANGES MADE TO FUNCTIONS AND CLASSES #########################
##################################################################################################################
##################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from openpmd_viewer import OpenPMDTimeSeries        # This should be used for pmd viewer version 1.x.
#from Useful_Functions import UsefulFunctions
from Particle_Properties import ParticleEnergy, ParticleTransverseProperties
import Useful_Functions as uf

from matplotlib import use
use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster.

fsize = 12

#%% Read data from files

FilePath = r'D:\Lewis\fbpic\20231206_FEBE_Plasma_Channel_110TW\20231208_FEBE_Matched_100mm_0p53Lp_Lab\diags\hdf5'

ts = OpenPMDTimeSeries(FilePath)

K = 8 # snapshot to analyse

x, y, z, ux, uy, uz, q, w = ts.get_particle( ['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'w'], species='beam',
                                                iteration=ts.iterations[K], plot=False)

ctau = round(ts.t[K]*c*1e3,2)

PE = ParticleEnergy(ux, uy, uz, q, w) # Create instance of class

#%% Test beam energy properties

# this line is useless
Ek_mean, de_rms, PC_spread, charge = PE.beam_energy_properties()


print("Mean energy: %0.2f MeV" % round(Ek_mean,2) )
print("RMS energy spread: %0.2f MeV" % round(de_rms,2) )
print("Percentage energy spread: %0.2f %%" % round(PC_spread,2) )
print("Tracked charge: %0.2f pC" % round(charge*-1e12,2) )

#%% Plot beam energy spectrum

e_spec, Spec_Min, Spec_Max, E_Bins = PE.get_energy_spectrum(0.1, 10)

Order_elec = np.floor(np.log10( np.amax( e_spec[0]) ))    # Order of magnitude of peak of spectrum
lin_max    = uf.roundup(np.amax(e_spec[0]),10**(Order_elec-1))

fig, ax = plt.subplots()
plt.plot(E_Bins[:-1], e_spec[0], '-', color='firebrick', linewidth=1.0)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.axis([Spec_Min, Spec_Max, 0, lin_max ])
ax.set_title(('c$\\tau$ = %0.2f mm' % ctau ), fontsize=fsize)
ax.set_xlabel('Beam energy (MeV)', fontsize=fsize)
ax.set_ylabel('Number of electrons', fontsize=fsize)
plt.grid(True)
plt.show()


#%% Create class for transverse particle properties in x and y

TX = ParticleTransverseProperties(x, ux, z, uz, w)
TY = ParticleTransverseProperties(y, uy, z, uz, w)

div_x, eta_tr_norm_x, beam_x, Twiss_alpha_x, Twiss_beta_x, Twiss_gamma_x = TX.transverse_beam_properties()
div_y, eta_tr_norm_y, beam_y, Twiss_alpha_y, Twiss_beta_y, Twiss_gamma_y = TY.transverse_beam_properties()
beam_z = TX.beam_length()
