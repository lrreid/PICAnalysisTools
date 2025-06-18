"""
Functions for getting chemical element names, atomic numbers and symbols.

To Do:
    - Add error handeling. What happens if you put in an invalid argument?

"""


import os
from numpy import array, where
from csv import reader

def get_element_name(element: str):
    """
    Get the name, atomic number and symbol for a given chemical element.
    Data taken from chemical_element_list.txt

    Parameters
    ----------
    element : String
        Atomic number, name or symbol for chemical element.
        If using atomic number, it is automatically converted to a string.

    Returns
    -------
    atomic_number : float
        Atomic number corresponding to chosen element.
    name : String
        Name of chosen element.
    symbol : String
        Symbol corresponding to chosen element.

    """
    
    filename = os.path.join( os.path.dirname(__file__), 'chemical_element_list.txt' )
    
    with open(filename, 'r') as f:
        csv_data = array([l for l in reader(f, delimiter='\t', skipinitialspace=True)])
    
    index         = where(csv_data == str(element))[0][0] # Get the row number which matches the chosen element
    atomic_number = csv_data[index,0]
    name          = csv_data[index,1]
    symbol        = csv_data[index,2]
    
    return float(atomic_number), name, symbol