
##################################################################

import numpy as np
import matplotlib.pyplot as plt

##################################################################

__all__ = ['energy_to_wavelength']

##################################################################

def energy_to_wavelength(energy):
    """
    Converts energy to wavelength for X-ray radiation. 
    :param energy: energy [eV]
    :return: wavelength [m]
    """
    return 1.23984193e-6/energy

##################################################################