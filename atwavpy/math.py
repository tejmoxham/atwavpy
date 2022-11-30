
##################################################################

import numpy as np
import matplotlib.pyplot as plt
from .utilities import energy_to_wavelength
from numpy.fft import ifftshift,fftshift,fftn,ifftn

##################################################################

__all__ = ['gaussian_function_1d']

##################################################################

def gaussian_function_1d(x, amplitude, mu, sigma):
    return amplitude*np.exp(-(x-mu)**2/(2.*sigma**2))

##################################################################

































