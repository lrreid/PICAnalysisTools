"""

"""

from os.path import join


def set_sim_path(FolderPath, Simulation, boosted_frame=False):

    if boosted_frame is True:
        FilePath = join(FolderPath, Simulation, 'lab_diags', 'hdf5')
    else:
        FilePath = join(FolderPath, Simulation, 'diags', 'hdf5')

    SimPath = join(FolderPath, Simulation)
   
    # FilePath = FilePath.replace('//', '/')      # Remove unceccessary Windows style double slashes
    # SimPath  = SimPath.replace('//', '/')       # Remove unceccessary Windows style double slashes

    return FilePath, SimPath