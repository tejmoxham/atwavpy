
##################################################################

import numpy as np
import matplotlib.pyplot as plt

##################################################################

__all__ = ['plot_intensity_phase', 'plot_beam_caustic', 'plot_peak_focal']

##################################################################

def plot_intensity_phase(field, show=False):
    xaxis = np.linspace(-field.xdim/2, field.xdim/2, field.nx)*1e6
    fig, [ax1,ax2] = plt.subplots(nrows=1, ncols=2, figsize=(8,5))
    fig.suptitle('Field Intensity & Phase', y=0.96, fontweight='bold', fontsize=16)
    ax1.plot(xaxis, np.abs(field.field)**2, 'b')
    ax1.set_xlabel('Position [$\mathrm{\mu m}$]', fontsize=16)
    ax1.set_ylabel('Intensity [$\mathrm{a.u}$]', fontsize=16)
    ax2.plot(xaxis, np.angle(field.field), 'b')
    ax2.set_xlabel('Position [$\mathrm{\mu m}$]', fontsize=16)
    ax2.set_ylabel('Phase [$\mathrm{rad}$]', fontsize=16)
    if show:
        plt.tight_layout()
        plt.show()

##################################################################

def plot_beam_caustic(field, show=False):
    ydim = field.caustic_xdim*1e6
    fig = plt.figure()
    plt.imshow(np.abs(field.caustic)**2, extent=[field.caustic_range[-1]*1e3, field.caustic_range[0]*1e3, -ydim/2, ydim/2], aspect='auto')
    fig.suptitle('Beam Caustic', y=0.96, fontweight='bold', fontsize=16)
    plt.xlabel('Longitudinal Position [$\mathrm{mm}$]', fontsize=16)
    plt.ylabel('Transverse Position [$\mathrm{\mu m}$]', fontsize=16)    
    if show:
        plt.tight_layout()
        plt.show()    

##################################################################

def plot_peak_focal(field, show=False):
    nz = field.caustic.shape[1]
    zrange = np.linspace(field.caustic_range[-1], field.caustic_range[0], nz)*1e3
    peak_intensity = np.max(np.abs(field.caustic)**2, axis=0)
    fig = plt.figure()
    plt.plot(zrange, peak_intensity, 'b')
    fig.suptitle('Peak Focal Position', y=0.96, fontweight='bold', fontsize=16)
    plt.xlabel('Longitudinal Position [$\mathrm{mm}$]', fontsize=16)
    plt.ylabel('Intensity [$\mathrm{a.u}$]', fontsize=16)   
    print("Peak Longitudinal Position:", zrange[peak_intensity.argmax()], 'mm')
    if show:
        plt.tight_layout()
        plt.show()

##################################################################