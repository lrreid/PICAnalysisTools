"""
Set paths to directories where data files are saved.
"""

from os.path import exists, join
from os import makedirs, listdir
from sys import exit

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

    check_for_h5_files(FilePath)                 # Stop script if no files are in simulation directory

    # print("FilePath: %s" % FilePath)
    # print("SimPath: %s" % SimPath)

    return FilePath, SimPath


def set_sim_path_from_cwd(Path_to_Sim: str, boosted_frame: bool = False):
    """
    Set the path to directory which contains .h5 files containing simulation data

    Parameters
    ----------
    Path_to_Sim : str
        Path to diectory of simulation containing input script
    boosted_frame : bool, optional
        If true, function will choose lab frame diagnostics directory (lab_diags), by default False

    Returns
    -------
    FilePath: str
        Full path to data files
    """

    if boosted_frame is True:
        FilePath = join(Path_to_Sim, 'lab_diags', 'hdf5')
    else:
        FilePath = join(Path_to_Sim, 'diags', 'hdf5')
   
    FilePath = FilePath.replace('//', '/')      # Remove unceccessary Windows style double slashes

    check_for_h5_files(FilePath)                # Stop script if no files are in simulation directory

    return FilePath

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

def check_for_h5_files(h5path: str):
    """
    Check if Simulation path directory contains any files. If not, the analysis stops.

    Parameters
    ----------
    h5path : string
        Full path to directory containing simulation files to analyse

    """

    if exists(h5path) is False:
        print("Chosen directory for simulation files does not exist\nAnalysis stopped.")
        exit()


    if len([i for i in listdir(h5path) if i.endswith('.h5')]) == 0:
        print("No .h5 files found in SimPath directory!\nAnalysis stopped.")
        exit()

    return None