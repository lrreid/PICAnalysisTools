
import numpy as np
from scipy.constants import c, pi
import matplotlib.pyplot as plt
from PICAnalysisTools.utils.ionisation import read_ionization_energies, ionisation_intensity_theshold
from PICAnalysisTools.utils.unit_conversions import get_order_letter
from PICAnalysisTools.utils.plot_limits import plt_limits, plt_limits_log
from PICAnalysisTools.utils.elements import get_element_name
from PICAnalysisTools.utils.laser_calcs import Gaussian_laser_intensitiy
fsize = 12

#%% Check ionisation energies function
potentials, ion_q = read_ionization_energies( "He", "eV" )  # Extract potentials from text file

print("Ionisation potentials of He:\n%s eV\n" % str(potentials) )

#%% Get ionisation intensity thresholds for each electron in a given atom

Element         = "N"
wavelength      = 800
wavelength_unit = "nano"
int_unit        = "centi"
full_unit       = "W%sm-2" % get_order_letter(int_unit)
Element_Name, Intensities, a0s,  = ionisation_intensity_theshold(Element, lambda0=wavelength, wavelength_unit=wavelength_unit, int_unit=int_unit)           # Calculate ionisation intensities for given atom

print("Element %s has full name %s" % (Element, Element_Name) )
print("Ionisation intensities for %s:" % Element_Name)
for i in range(len(Intensities)):
    print("%d+\t\t%0.2e %s\t\t%0.3f" % ((i+1), Intensities[i], full_unit, np.round(a0s[i],3) ))
print("\n\n")

#%% Plot intesnsity thresholds for range of atoms compared to laser peak intensity

save_plots         = False
show_laser         = True
Exp_laser_name     = "Laser peak intensity" 
save_string        = "ionisation_thresholds_test" 
plot_label         = "element"                              # Choose whether element name or symbol is used in plot legend. element or symbol.
Ana_dir            = "path/to/analysis/directory"

Elements = ["H", "He", "N", "O", "Ne", "Ar"]
Elements.reverse()                                          # Not required but improves visibility of smaller atoms in plot

# Set of laser parameters to compare thresholds to. 
E           = 3.00                                          # Pulse energy (J).
energy_unit = ""
w0          = 30                                            # Laser waist to 1/e^2 radius (um)
spot_unit   = "micro"
tau_FWHM    = 30                                            # FWHM pulse duration (fs)
time_unit   = "femto"
laser_a0, laser_I = Gaussian_laser_intensitiy(Energy=E, tau_FWHM=tau_FWHM, w0=w0, lambda0=wavelength, energy_unit=energy_unit, time_unit = time_unit, spot_unit=spot_unit, wavelength_unit=wavelength_unit, int_unit=int_unit)

if show_laser is True:
    title_string = "Laser ionisation intensity theshold of atoms\n$\\lambda_{0}$ = %d %sm, E = %0.2f %sJ, $\\omega_{0}$ = %0.2f %sm, $\\tau_{FWHM}$ = %d %ss" % (int(wavelength), get_order_letter(wavelength_unit), np.round(E,2), get_order_letter(energy_unit), np.round(w0,2), get_order_letter(spot_unit, True), int(tau_FWHM), get_order_letter(time_unit) )
else:
    title_string = "Laser ionisation intensity theshold of atoms"

intensity_list = [0] * len(Elements)
a0_list        = [0] * len(Elements)
name_list      = [0] * len(Elements)
atomic_nos     = [0] * len(Elements)

for k in range(len(Elements)):
    name_list[k], intensity_list[k], a0_list[k] = ionisation_intensity_theshold(Elements[k], lambda0=wavelength, wavelength_unit=wavelength_unit, int_unit=int_unit)
    atomic_nos[k] = get_element_name(Elements[k])[0]

y_min, y_max = plt_limits_log(np.concatenate((intensity_list)), min_offset = 0, max_offset = 1)   
x_min, x_max = plt_limits(atomic_nos, 5)
x_points     = np.arange(0, x_max+2, 1)

if plot_label == "symbol":
    label_list = Elements
elif plot_label == "element":
    label_list = name_list
else:
    print("plot lable incorrectly defined. Element symbol used as default")
    label_list = Elements

colours = ['lime', 'darkorange', 'seagreen', 'magenta', 'gold', 'cyan', 'firebrick', 'dodgerblue', 'mediumblue']

fig, ax = plt.subplots(dpi=150)
if show_laser is True:
    plt.plot([0,x_max+2], [laser_I, laser_I], '-', color='red', linewidth=2, label="%s" % Exp_laser_name)

for i in range(len(Elements)):
    plt.plot(np.arange(0,len(intensity_list[i]),1)+1, intensity_list[i], '-o', color=colours[i], linewidth=2, label=("%s" % label_list[i] ) )

plt.yscale('log')
ax.set_yticks([1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22, 1e23, 1e24])
ax.set_xticks(np.arange(0,x_max+2,2))
plt.axis([x_min, x_max, y_min, y_max])         # set axes [xmin, xmax, ymin, ymax]
ax.set_title("%s" % title_string , fontsize=fsize)
ax.set_xlabel('Ion charge state', fontsize=fsize)
ax.set_ylabel('Intensity (W%sm$^{-2}$)' % (get_order_letter(int_unit)), fontsize=fsize)
plt.grid(True, which="both", ls="-")
plt.legend(loc='lower right', prop={'size': 8}, ncol = 2)

if save_plots is True:
    plt.savefig("%s/Laser_ionisation_thresholds_%s.png" % (Ana_dir, save_string), dpi=300, format="png")
    plt.close()
else:
    plt.show()