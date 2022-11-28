
##################################################################

import copy
import numpy as np
import matplotlib.pyplot as plt
from .utilities import energy_to_wavelength
from .prop import fraunhofer_propagation_1d, fresnel_propagation_1d, fresnel_kirchhoff_propagator_1d, angular_spectrum_rayleigh_sommerfeld_1d

##################################################################

__all__ = ['field1d']

##################################################################

class field1d:

    def __init__(self, nx, energy, pixel, zpos=0):
        self.nx = nx
        self.energy = energy
        self.pixel = pixel
        self.zpos = zpos
        self.wavelength = energy_to_wavelength(energy)
        self.xdim = nx*pixel
        self.field = np.ones(nx, dtype=np.complex128)

    def propagate(self, distance, method='fresnel', **args):
        if method=='fraunhofer':
            self.pixel, self.field = fraunhofer_propagation_1d(self.field, self.energy, self.pixel, distance, **args)
            self.zpos += distance
            self.nx = self.field.size
        elif method=='fresnel':
            self.pixel, self.field = fresnel_propagation_1d(self.field, self.energy, self.pixel, distance, **args)
            self.zpos += distance
            self.nx = self.field.size
        elif method=='fresnel-kirchoff':
            self.pixel, self.field = fresnel_kirchhoff_propagator_1d(self.field, self.energy, self.pixel, distance, **args)
            self.zpos += distance
            self.nx = self.field.size  
        elif method=='angular-spectrum':
            self.pixel, self.field = angular_spectrum_rayleigh_sommerfeld_1d(self.field, self.energy, self.pixel, distance, **args)
            self.zpos += distance
            self.nx = self.field.size  
        return self
    
    def beam_caustic(self, nz, prange, method='fresnel', start=0, verbose=0, **args):
        if method=='fraunhofer':
            propagator = lambda distance : fraunhofer_propagation_1d(self.field, self.energy, self.pixel, distance, **args)
        elif method=='fresnel':
            propagator = lambda distance : fresnel_propagation_1d(self.field, self.energy, self.pixel, distance, **args)
        elif method=='fresnel-kirchoff':
            propagator = lambda distance : fresnel_kirchhoff_propagator_1d(self.field, self.energy, self.pixel, distance, **args)
        elif method=='angular-spectrum':
            propagator = lambda distance : angular_spectrum_rayleigh_sommerfeld_1d(self.field, self.energy, self.pixel, distance, **args)
        prop_range = start+np.linspace(prange[0], prange[-1], nz)
        for i,distance in enumerate(prop_range):
            if verbose > 0:
                print('iteration:',str(i+1)+'/'+str(nz),'- propagating:',distance,'m')
            pixep, field = propagator(distance)
            if i==0:
                caustic =  np.zeros([nz, field.size], dtype=np.complex)
            caustic[i] = field
        self.caustic = caustic.T
        self.caustic_pixel = pixep
        self.caustic_range = np.array(prange)
        self.caustic_xdim = pixep*caustic.T.shape[0]
        return self

    def copy(self):
        return copy.deepcopy(self)

##################################################################