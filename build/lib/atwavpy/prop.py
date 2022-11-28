
##################################################################

import cv2
import numpy as np
import matplotlib.pyplot as plt
from .utilities import energy_to_wavelength
from numpy.fft import ifftshift,fftshift,fftn,ifftn

##################################################################

__all__ = ['fraunhofer_propagation_1d', 'fresnel_propagation_1d', 
'fresnel_kirchhoff_propagator_1d', 'angular_spectrum_rayleigh_sommerfeld_1d']

##################################################################

def fraunhofer_propagation_1d(efield, energy, pixel, distance, pfactor=1.0, refindx=1):
    """
    Propagates a 1d complex electric field using the Fraunhofer (far-field) approximation. 
    :param efield: field to propagate (1d complex array)
    :param energy: energy (float) [ev]
    :param pixel: pixel size (float) [m]
    :param distance: propagation distance (float) [m]
    :param pfactor: pixel size scaling factor
    :param refindx: refractive index of propagation medium (float)
    :return: pixel size in propagation plane (float) [m]
    :return: propagated field (1d complex array)
    """
    nx = efield.size
    if pfactor<1:
        efield = np.pad(efield, (int((nx*(1/pfactor)-nx)/2),int((nx*(1/pfactor)-nx)/2)), mode='constant')
    elif pfactor>1:
        efield = efield[int((1-(1/pfactor))*nx/2):-int((1-(1/pfactor))*nx/2)]
    nx = efield.size
    lam = energy_to_wavelength(energy)
    vx = np.arange(-nx//2,nx//2)
    pixelp = lam*np.abs(distance)/nx/pixel
    if distance>0:
        # Forward Propagation
        factor = (refindx*np.pi)/(lam*np.abs(distance))*pixel**2
        quad_factor = np.exp(1j*((factor*vx**2)))
        efieldp = ifftshift(fftn(fftshift(quad_factor*efield),norm='ortho'))
    elif distance<0:
        # Backward Propagation    
        factor = (refindx*np.pi)/(lam*np.abs(distance))*pixelp**2
        quad_factor = np.exp(-1j*((factor*vx**2)))
        efieldp = quad_factor*fftshift(ifftn(ifftshift(efield),norm='ortho'))
    return pixelp,efieldp

##################################################################

def fresnel_propagation_1d(efield, energy, pixel, distance, pfactor=1.0, dfactor=1.0, refindx=1.0):
    """
    Propagates a 1d complex electric field using the Frensel (near-field) approximation. 
    :param efield: field to propagate (1d complex array)
    :param energy: energy (float) [ev]
    :param pixel: pixel size (float) [m]
    :param distance: propagation distance (float) [m]
    :param pfactor: pixel size scaling factor (float)
    :param dfactor: dimension size scaling factor (float)
    :param refindx: refractive index of propagation medium (float)
    :return pixelp: pixel size in propagation plane (float) [m]
    :return efieldp: propagated field (1d complex array)
    """
    nx = efield.size
    if dfactor>1:
        efield = np.pad(efield, (int((nx*dfactor-nx)/2),int((nx*dfactor-nx)/2)), mode='constant')
    elif dfactor<1:
        efield = efield[int((1-dfactor)*nx/2):-int((1-dfactor)*nx/2)]
    nx = efield.size
    lam = energy_to_wavelength(energy)
    min_dist = max(np.abs((pfactor-1)/pfactor),np.abs(pfactor-1))*min(nx,nx)*min(pixel,pixel)**2/lam
    max_dist = np.abs(pfactor)*min(nx,nx)*min(pixel,pixel)/lam
    if (np.abs(distance)>max_dist) or (np.abs(distance)<min_dist):
        print('Propagation distance must be between |{:0.2e}|<z<|{:0.2e}|!'.format(min_dist,max_dist))
    pixelp = pfactor*pixel
    f1 = np.pi/(lam*distance)*(1-pfactor)*pixel**2
    f2 = -np.pi*lam*distance/pfactor/(nx*pixel)**2
    f3 = np.pi/(lam*distance)*pfactor*(pfactor-1)*pixel**2
    efield = ifftshift(efield)
    vx = ifftshift(np.arange(-nx//2,nx//2))
    quad = lambda x: np.exp(1j*refindx*(x*vx**2))
    efieldp = fftshift(quad(f3)*ifftn(quad(f2)*fftn(quad(f1)*efield, norm='ortho'), norm='ortho'))
    return pixelp,efieldp

##################################################################

def fresnel_kirchhoff_propagator_1d(field, energy, pixel, distance, pfactor=1.0, sfactor=1.0, refindx=1):
    """
    Propagates a 1d complex electric field using the Frensel-Kirchoff equation. 
    :param field: field to propagate (1d complex array)
    :param energy: energy (float) [ev]
    :param pixel: pixel size (float) [m]
    :param distance: propagation distance (float) [m]
    :param pfactor: pixel size scaling factor (float)
    :param sfactor: sampling size scaling factor (float)
    :param refindx: refractive index of propagation medium (float)
    :return: pixel size in propagation plane (float) [m]
    :return: propagated field (1d complex array)
    """
    lam = energy_to_wavelength(energy)
    k = 2.0*refindx*np.pi/lam if distance>0 else -2.0*refindx*np.pi/lam 
    nx = field.size
    x = np.arange(0, float(nx))-(nx/2)
    x *= pixel
    nxp = int(nx*sfactor)
    pixelp = pixel*pfactor
    fieldp = np.ones([nxp], dtype=np.complex128)
    xp = np.arange(0, float(nxp))-(nxp/2)
    xp *= pixelp
    """
    # Pure Numpy
    x1 = np.outer(x, np.ones(nxp))
    x2 = np.outer(np.ones(nx), xp)
    plen = np.sqrt(np.power(x1-x2,2)+np.power(distance,2))
    fieldp = np.dot(field, (np.exp(1j*k*plen[:,])))
    """
    fieldp = np.zeros(nxp, dtype=np.complex)
    for i in range (nxp):
        plen = np.sqrt(np.power(x-xp[i],2)+np.power(distance,2))
        fieldp[i] = np.sum(field*(np.exp(1j*k*plen)))
    fieldp *= pixel*np.sqrt(pfactor)/np.sqrt(lam*distance)
    return pixelp,fieldp

##################################################################

def angular_spectrum_rayleigh_sommerfeld_1d(efield, energy, pixel, distance, dfactor=1.0, refindx=1):
    """
    Propagates a 1d complex electric field using the Rayleigh-Sommerfeld angular specturm method. 
    :param efield: field to propagate (1d complex array)
    :param energy: energy (float) [ev]
    :param pixel: pixel size (float) [m]
    :param distance: propagation distance (float) [m]
    :param dfactor: dimension size scaling factor (float)
    :param refindx: refractive index of propagation medium (float)
    :return: pixel size in propagation plane (float) [m]
    :return efieldp: propagated field (1d complex array)
    """
    nx = efield.size
    if dfactor>1:
        efield = np.pad(efield, (int((nx*dfactor-nx)/2),int((nx*dfactor-nx)/2)), mode='constant')
    elif (dfactor is not None) and (dfactor<1):
        efield = efield[int((1-dfactor)*nx/2):-int((1-dfactor)*nx/2)]
    nx = efield.size
    lam = energy_to_wavelength(energy)
    k = 2.0*refindx*np.pi/lam if distance>0 else -2.0*refindx*np.pi/lam 
    ix = np.arange(-nx//2,nx//2)
    uu = ix*(lam/(pixel*nx)) 
    ang_mesh = np.sqrt((1-uu**2))
    ang_mesh = np.exp(1j*k*np.abs(distance)*ang_mesh)
    if distance>0:
        # Forward Propagation 
        efieldp = ifftn(ifftshift(ang_mesh)*fftn(efield))
    elif distance<0:
        # Backward Propagation 
        efieldp = fftn(ifftshift(ang_mesh)*ifftn(efield))
    elif distance==0:
        # No Propagation 
        efieldp = efield
    return pixel,efieldp

##################################################################

































