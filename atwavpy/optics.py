
##################################################################

import numpy as np
from .field import field1d
import matplotlib.pyplot as plt

##################################################################

__all__ = ['aperture_1d', 'ideal_lens_1d']

##################################################################

def aperture_1d(field, size, offset=None):
    nx = field.nx
    if offset is None:
        center = int(nx/2)
    size = int(size/(2*field.pixel))
    field.field[:center-size] = 0
    field.field[center+size:] = 0
    return field

def ideal_lens_1d(field, focal_length):
    nx = field.nx
    lam = field.wavelength
    factor = -np.pi*field.pixel**2/(focal_length*lam)
    ix = np.arange(-nx//2,nx//2)
    field.field *= np.exp(1j*((factor*ix**2)))
    return field

##################################################################