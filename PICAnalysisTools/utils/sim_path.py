"""

"""

from os.path import exists, join
from os import makedirs


def set_sim_path(FolderPath, Simulation, boosted_frame=False):

    if boosted_frame is True:
        FilePath = join(FolderPath, Simulation, 'lab_diags', 'hdf5')
    else:
        FilePath = join(FolderPath, Simulation, 'diags', 'hdf5')

    SimPath = join(FolderPath, Simulation)
   
    # FilePath = FilePath.replace('//', '/')      # Remove unceccessary Windows style double slashes
    # SimPath  = SimPath.replace('//', '/')       # Remove unceccessary Windows style double slashes

    return FilePath, SimPath

def set_Analysis_path(FolderPath, Simulation, Ana_name):

    Analysis_path = join(FolderPath, Simulation, 'Analysed', '%s' % Ana_name)

    if exists(Analysis_path) is False:
        makedirs(Analysis_path)

    return Analysis_path
