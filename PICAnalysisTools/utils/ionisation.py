'''
Get ionisation potentials and intensity thresholds of atoms.

read_ionization_energies() adapted from a function of the same name in fbpic (read_atomic_data.py) fbpic.particles.elementary_process.ionization read_atomic_data.py
'''

import numpy as np
import re
import os
from scipy.constants import c, e, pi, epsilon_0, m_e
from PICAnalysisTools.utils.elements import get_element_name
from PICAnalysisTools.utils.laser_calcs import a0_from_intensity
from PICAnalysisTools.utils.unit_conversions import magnitude_conversion_area

def read_ionization_energies( element: str, unit: str = "eV" ):
    """
    Read the ionization energies from a data file

    Parameters
    ----------
    element: string
        The atomic symbol of the considered ionizable species
        (e.g. 'He', 'N' ;  do not use 'Helium' or 'Nitrogen')
        
    unit: string
        The unit that the data is returned in. Choose "eV" or "Joule"
        Electron volts are default.

    Returns
    -------
    energies: array_like
        Ionisation potentials of each electron in chosen atom
    ion_charge: float
        Ion charge state corresponding to each ionisation energy
    """
    # Open and read the file atomic_data.txt
    filename = os.path.join( os.path.dirname(__file__), 'atomic_data.txt' )
    with open(filename) as f:
        text_data = f.read()
    # Parse the data using regular expressions (a.k.a. regex)
    # (see https://docs.python.org/2/library/re.html)
    # The regex command below parses lines of the type
    # '\n     10 | Ne IV         |         +3 |           [97.1900]'
    # and only considers those for which the element (Ne in the above example)
    # matches the element which is passed as argument of this function
    # For each line that satisfies this requirement, it extracts a tuple with
    # - the atomic number (represented as (\d+))
    # - the ionization level (represented as the second (\d+))
    # - the ionization energy (represented as (\d+\.*\d*))
    regex_command = \
        '\n\s+(\d+)\s+\|\s+%s\s+\w+\s+\|\s+\+*(\d+)\s+\|\s+\(*\[*(\d+\.*\d*)' \
        %element
    list_of_tuples = re.findall( regex_command, text_data )
    # Return None if the requested element was not found
    if list_of_tuples == []:
        return(None)
    # Go through the list of tuples and fill the array of ionization energies.
    atomic_number = int( list_of_tuples[0][0] )
    assert atomic_number > 0
    energies = np.zeros( atomic_number )
    ion_charge = np.zeros( atomic_number )
    for ion_level in range( atomic_number ):
        # Check that, when reading the file,
        # we obtained the correct ionization level
        assert ion_level == int( list_of_tuples[ion_level][1] )
        
        if unit == "Joule":
            # Get the ionization energy and convert in Joules using e
            energies[ ion_level ] = e * float( list_of_tuples[ion_level][2] )
        elif unit == "eV":
            energies[ ion_level ] = float( list_of_tuples[ion_level][2] )
        else:
            print("Warning: Unit incorrectly defined, electronvolts are used as default")
            energies[ ion_level ] = float( list_of_tuples[ion_level][2] )
            
        ion_charge[ ion_level ] = float( list_of_tuples[ion_level][1] ) # Save list of ion charges

    return( energies, ion_charge )


def ionisation_intensity_theshold(element: str, lambda0: float = 800, wavelength_unit: str = "nano", int_unit: str = "centi"):
    """
    Calculate the threshold laser intensity to ionise electrons off atoms

    Parameters
    ----------
    element: str
        The atomic name, symbol or atomic number of the considered ionizable species.
    lambda0: float
        central wavelength of laser radiation. Used to find a0 equivalent to intensity, by default "nano"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    int_unit : str, optional
        Order of magnitude of laser intensity unit, by default "centi"

    Returns
    -------
    element_name: string
        Full name of chosen element.
    Int : float
        Laser intensity threshold for over the barrier ionisation. Default unit: Wcm^-2
    a0  : float
        Normalised laser vector potential. Unit: dimentionless
    """

    element_name = get_element_name(str(element))       # Get element name from symbol
    Ip, ion_q    = read_ionization_energies( element_name[2], unit = "eV" )
    
    Z   = ion_q + 1                                                         # final charge state of the ion
    Int = ((pi**2 * c * epsilon_0**3 * (Ip*e)**4) / (2* Z**2 * e**6) )      # intensity Wm^-2
    a0  = a0_from_intensity(Int, lambda0=lambda0, int_unit="", wavelength_unit=wavelength_unit)
    
    return element_name[1], magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0