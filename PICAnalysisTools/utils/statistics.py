'''
Statistical calculations
'''


import numpy as np

def w_std( a, weights ):                                # Define weighted standard deviation calculation
    average = np.average(a, weights=weights)
    variance = np.average((a - average) ** 2, weights=weights)
    return( average, np.sqrt(variance) )