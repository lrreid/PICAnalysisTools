"""
Set paths to directories where data files are saved.
"""

from os.path import exists, join
from os import makedirs


def set_sim_path(FolderPath: str, Simulation: str, boosted_frame: bool = False):
    """
    Set the path to directory which contains .h5 files containing simulation data

    Parameters
    ----------
    FolderPath : str
        Path to parent diectory of simulations
    Simulation : str
        Name of simulation 
    boosted_frame : bool, optional
        If true, function will choose lab frame diagnostics directory (lab_diags), by default False

    Returns
    -------
    FilePath: str
        Full path to data files
    SimPath: str
        Full path to simulation directory
    """

    if boosted_frame is True:
        FilePath = join(FolderPath, Simulation, 'lab_diags', 'hdf5')
    else:
        FilePath = join(FolderPath, Simulation, 'diags', 'hdf5')

    SimPath = join(FolderPath, Simulation)
   
    FilePath = FilePath.replace('//', '/')      # Remove unceccessary Windows style double slashes
    SimPath  = SimPath.replace('//', '/')       # Remove unceccessary Windows style double slashes

    return FilePath, SimPath

def set_analysis_path(SimPath: str, Ana_name: str):
    """
    Set and create directory where analysed data is saved

    Parameters
    ----------
    SimPath : str
        Path to simulation directory
    Ana_name : str
        Name of directory where data will be saved

    Returns
    -------
    Analysis_path: str
        Full path to analysis directory
    """

    Analysis_path = join(SimPath, 'Analysed', '%s' % Ana_name)

    if exists(Analysis_path) is False:
        makedirs(Analysis_path)

    return Analysis_path
