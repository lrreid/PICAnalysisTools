# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:44:01 2024

@author: ryi76833
"""

from openpmd_viewer.addons import LpaDiagnostics
from Laser_Properties import LaserProperties

fsize = 12

#%% Read data from files

FilePath = r'D:\Lewis\fbpic\20230830_FEBE_self-injection\20230830_2806mJ_23fs_w0_26um_ne_020e17\diags\hdf5'

ts = LpaDiagnostics(FilePath, check_all_files=False, backend='h5py')
LP = LaserProperties(ts)

out = LP.get_laser_properties_loop(coord='x', method='fit', profile_Method='projection')