
##################################################################



##################################################################


























##################################################################

def fraunhofer_propagation_2d(field,energy,pixel_size,distance):
    field = np.fft.fftshift(field)
    nyx = np.array(field.shape)
    lam = energy_to_wavelength(energy)
    if (type(pixel_size) is float) or (type(pixel_size) is int):
        pixel_size = [pixel_size,pixel_size]
    pixel_size = np.array(pixel_size)
    ix = np.fft.fftshift(np.arange(-nyx[1]//2,nyx[1]//2))
    iy = np.fft.fftshift(np.arange(-nyx[0]//2,nyx[0]//2))
    vy,vx = np.meshgrid(iy,ix,indexing='ij')
    pixel_size_new = lam*np.abs(distance)/nyx/pixel_size
    if distance>0:
        factor = np.pi/(lam*distance)*pixel_size**2
        quad = np.exp(1j*((factor[1]*vx**2)+(factor[0]*vy**2)))
        field_prop = np.fft.fftshift(quad*np.fft.fft2(field))
    elif distance<0:
        factor = np.pi/(lam*distance)*pixel_size_new**2
        quad = np.exp(1j*((factor[1]*vx**2)+(factor[0]*vy**2)))
        field_prop = np.fft.fftshift(quad*np.fft.ifft2(field))
    return pixel_size_new,field_prop

##################################################################

def fresnel_propagation_2d(field,energy,pixel_size,distance,mag=1.0):
    field = np.fft.fftshift(field)
    nyx = np.array(field.shape)
    lam = energy_to_wavelength(energy)
    if (type(pixel_size) is float) or (type(pixel_size) is int):
        pixel_size = [pixel_size,pixel_size]
    pixel_size = np.array(pixel_size)
    min_dist = max(np.abs((mag-1)/mag),np.abs(mag-1))*min(nyx)*min(pixel_size)**2/lam
    max_dist = np.abs(mag)*min(nyx)*min(pixel_size)/lam
    if (np.abs(distance)>max_dist) or (np.abs(distance)<min_dist):
        print('Propagation distance must be between |{:0.2e}|<z<|{:0.2e}|!'.format(min_dist,max_dist))
    ix = np.fft.fftshift(np.arange(-nyx[1]//2,nyx[1]//2))
    iy = np.fft.fftshift(np.arange(-nyx[0]//2,nyx[0]//2))
    vy,vx = np.meshgrid(iy,ix,indexing='ij')
    quad = lambda x: np.exp(1j*((x[1]*vx**2)+(x[0]*vy**2)))
    f1 = np.pi/(lam*distance)*(1-mag)*pixel_size**2
    f2 = -np.pi*lam*distance/mag/(nyx*pixel_size)**2
    f3 = np.pi/(lam*distance)*mag*(mag-1)*pixel_size**2
    field_prop = quad(f3)*np.fft.ifft2(quad(f2)*np.fft.fft2(quad(f1)*field))
    pixel_size_new = mag*pixel_size
    field_prop = np.fft.fftshift(field_prop)
    pixel_size_new = mag*pixel_size
    return pixel_size_new,field_prop

##################################################################
# Propagator Functions
##################################################################




##################################################################

























def angular_spectrum_propagation_2d(array,energy,pixel_size,distance,dim_factor=1):
    if dim_factor!=1:
        ny,nx = array.shape
        padding = [int(ny*(dim_factor-1)),int(nx*(dim_factor-1))]
        array = np.pad(array,((padding[0],padding[1]),(padding[0],padding[1])),'constant')
    if (type(pixel_size) is float) or (type(pixel_size) is int):
        pixel_size = [pixel_size,pixel_size]
    pixel_size = np.array(pixel_size)
    wavelength = energy_to_wavelength(energy)
    du = 1/(pixel_size*np.size(array,0))
    ii, jj = np.float64(np.meshgrid(range(0,np.size(array,1)),range(0,np.size(array,0))))
    ii -= np.size(array,1)/2
    jj -= np.size(array,0)/2
    ii *= du[1]; jj *= du[0]
    ii *= wavelength; jj *= wavelength
    ang_mesh = np.square(ii)+np.square(jj)
    ang_mesh = np.sqrt(1-ang_mesh)
    ang_mesh = np.exp(1j*(2*np.pi*distance/wavelength)*ang_mesh)
    ang_mesh = np.fft.fftshift(ang_mesh,[-2,-1])
    prop_array = np.fft.ifftn(ang_mesh*np.fft.fftn(array,None,[-2,-1],norm='ortho'),norm='ortho')
    return prop_array

##################################################################

def modal_angular_spectrum_propagation(probe,energy,pixel_size,distance,dim_factor=1):
    for i,mode in enumerate(probe):
        prop_mode = angular_spectrum_propagation_2d(mode,energy,pixel_size,distance,dim_factor)
        if i==0:
            total_intensity = np.zeros_like(prop_mode,dtype=np.float)
        total_intensity += np.abs(prop_mode)**2
    return total_intensity

##################################################################
'''
def fresnel_propagation_2d(electric_field,energy,pixel_size,distance,dimension=1,magnification=1):
    if type(pixel_size) is float:
        pixel_size = [pixel_size,pixel_size]
    ny,nx = electric_field.shape
    padding = [int(ny*(dimension-1)),int(nx*(dimension-1))]
    electric_field = np.pad(electric_field,((padding[0],padding[1]),(padding[0],padding[1])),'constant')
    wavelength = energy_to_wavelength(energy)
    ny,nx = electric_field.shape
    min_dist = max(np.abs((magnification-1)/magnification),np.abs(magnification-1))*nx*pixel_size[0]**2/wavelength
    max_dist = np.abs(magnification)*nx*magnification**2/wavelength
    if np.abs(distance)<min_dist:
        print('Propagation distance must be greater!')
    elif np.abs(distance)>max_dist:
        print('Propagation distance must be smaller!')
    k = 2*np.pi/wavelength
    freq_nyquist = 0.5/pixel_size[1]
    freq_n = np.linspace(-1.0, 1.0, nx)
    freq_x = freq_n*freq_nyquist
    freq_nyquist = 0.5/pixel_size[0]
    freq_n = np.linspace(-1.0, 1.0, ny)
    freq_y = freq_n*freq_nyquist
    shift_half_pixel = True
    if shift_half_pixel:
        freq_x = freq_x-0.5*np.abs(freq_x[1]-freq_x[0])
        freq_y = freq_y-0.5*np.abs(freq_y[1]-freq_y[0])
    f_x, f_y = np.meshgrid(freq_x,freq_y, indexing='ij')
    fsq = np.fft.fftshift(f_x**2+f_y**2)
    xdim,ydim = nx*pixel_size[1],ny*pixel_size[0]
    xaxis = np.linspace(-xdim/2,xdim/2,nx)
    yaxis = np.linspace(-ydim/2,ydim/2,ny)
    x = np.meshgrid(xaxis,xaxis)[0].T
    y = np.meshgrid(yaxis,yaxis)[1].T
    x_rescaling = np.meshgrid(xaxis,xaxis)[0].T*magnification
    y_rescaling = np.meshgrid(yaxis,yaxis)[1].T*magnification
    r1sq = x**2+y**2
    r2sq = x_rescaling**2+y_rescaling**2
    Q1 = k/2*(1-magnification)/distance*r1sq
    Q2 = np.exp(-1.0j*np.pi*wavelength*distance/magnification*fsq)
    Q3 = np.exp(1.0j*k/2*(magnification-1)/(magnification*distance)*r2sq)
    electric_field *= np.exp(1.0j*Q1)
    fft = np.fft.fft2(electric_field)
    return np.fft.ifft2(fft*Q2)*Q3/magnification
'''
##################################################################

def modal_fresnel_propagation(probe,energy,pixel_size,distance,dimension=1,magnification=1):
    for i,mode in enumerate(probe):
        print(i)
        prop_mode = fresnel_propagation_2d(mode,probe,energy,pixel_size,distance,dimension,magnification)
        if i==0:
            total_intensity = np.zeros_like(prop_mode,dtype=np.float)
        total_intensity += np.abs(prop_mode)**2
    return total_intensity

##################################################################















##################################################################

def single_probe_propagation_2d(probe,energy,pixel_size,distance,dimension=1,magnification=1,method='fresnel'):
    probe = np.squeeze(probe)
    if method=='angular_spectrum':
        probe_propagated = angular_spectrum_propagation_2d(probe,energy,pixel_size,distance,dimension)
    elif method=='fresnel':
        probe_propagated = fresnel_propagation_2d(probe,energy,pixel_size,distance,dimension,magnification)
    return probe_propagated

##################################################################

def modal_probe_propagation_2d(probe,energy,pixel_size,distance,dimension=1,magnification=1,method='fresnel'):
    power,probe = ortho_multimodal_array(probe)
    if method=='angular_spectrum':
        for i,mode in enumerate(probe):
            print('Propagating Mode:',i)
            probe_propagated = angular_spectrum_propagation_2d(mode,energy,pixel_size,distance,dimension)
            if i==0:
                intensity_propagated = np.zeros_like(probe_propagated,dtype=np.float)
            intensity_propagated += power[i]*np.abs(probe_propagated)**2
    elif method=='fresnel':
        for i,mode in enumerate(probe):
            print('Propagating Mode:',i)
            probe_propagated = fresnel_propagation_2d(mode,energy,pixel_size,distance,dimension,magnification)
            if i==0:
                intensity_propagated = np.zeros_like(probe_propagated,dtype=np.float)
            intensity_propagated += power[i]*np.abs(probe_propagated)**2
        pixel_size *= magnification
    return intensity_propagated


##################################################################

def exact_prop(in_wave,out_wave,L_in,L_out,energy,z):
    wavel = energy_to_wavelength(energy)
    pi = np.pi
    N_in = np.shape(in_wave)[0]
    in_domain = np.linspace(-L_in/2,L_in/2,N_in)
    N_out = np.shape(out_wave)[0]
    out_domain = np.linspace(-L_out/2,L_out/2,N_out)
    step_in = L_in/N_in
    for i in range(N_out):
        for j in range(N_in):
            x = in_domain[j]
            f = in_wave[j]
            x1 = out_domain[i]
            out_wave[i] += f*np.exp((-1j*pi*x*x)/(wavel*z))*np.exp((-1j*2*pi*x*x1)/(wavel*z))
            #out_wave[i] += ne.evaluate('f*exp((-1j*pi*x*x)/(wavel*z))*exp((-1j*2*pi*x*x1)/(wavel*z))')
    out_wave *= (1/np.sqrt(1j*wavel*z))*step_in
    return out_wave


##################################################################

order = 2
if order == 0:
    def path_length (delx, d):
        return (0.0)
elif order == 1:
    def path_length (delx, d):
        return (delx**2/d/2.0)
elif order == 2:
    def path_length (delx, d):
        xsq = (delx / d)**2
        return (d * (xsq/2 + xsq**2/8.0))
elif order == 3:
    def path_length (delx, d):
        xsq = (delx / d)**2
        return (d * (xsq/2 + xsq**2/8.0 + xsq**3/16.0))
#
def propagation_phase_shift (lam, r):
    return (2.0*np.pi*r/lam)

def cis (z):
    return (np.exp (-1j * z))

def fresnel_kirchoff_propagator_1d(efield, energy, xaxis, xaxisp, d):
    nx = xaxis.size
    xmin = xaxis[0]
    xmax = xaxis[-1]
    dx = np.abs(xmax-xmin)/nx
    nxp = xaxisp.size
    xminp = xaxisp[0]
    xmaxp = xaxisp[-1]
    dxp = np.abs(xmaxp-xminp)/nxp
    lam = energy_to_wavelength(energy)
    prop = np.zeros(nxp,dtype=np.complex)

    for j in range(nxp):
        x1 = xaxisp[j]
        esum = 0.0
        for k in range(nx):
            x2 = xaxis[k]
            plen = path_length(x2-x1,d) #np.sqrt((x2-x1)**2+d**2) # p
            phase = cis(propagation_phase_shift(lam,plen)) #np.exp(1j*2*np.pi*plen/lam)#/(plen**2) #cis(propagation_phase_shift(lam,plen))
            esum += efield[k]*phase
        enew = esum * dx /np.sqrt(lam * d)
        prop[j] = enew
    return prop
    

##################################################################




def refine_focal_position(energy,efield,xdim,nx,xdimp,nxp,distance,range,nz,p0=[1,0,1e-6]):
    long_range = np.linspace(range[0],range[-1],nz)
    for i,position in enumerate(long_range):
        print(position)
        efieldp = fresnel_kirchoff_propagator_1d(energy,efield,xdim,nx,xdimp,nxp,distance+position)
        xaxisp = np.linspace(-xdimp/2,xdimp/2,nxp)
        intensity = np.abs(efieldp)**2
        norm_intensity = intensity/intensity.max()
        try:
            popt, pcov = curve_fit(gaussian_function_1d,xaxisp,norm_intensity,p0)
            print('Focal Size FWHM:',round(popt[2]*2.355*1e9,2),'nm')
            plt.plot(xaxisp*1e6,gaussian_function_1d(xaxisp,*popt),'r--')
        except:
            pass
        plt.plot(xaxisp*1e6,norm_intensity,'b')
        plt.title("Fresnel-Kirchhoff Diffraction")
        plt.xlabel("Position (um)")
        plt.ylabel("Normalised Intensity [a.u.]")
        plt.show()

##################################################################




















##################################################################

def angular_spectrum_propagation_1d(array,energy,pixel_size,distance):

    u = array
    z = -1*distance

    dx = pixel_size

    wavelength = energy_to_wavelength(energy)
    wavenumber = 2*np.pi/wavelength
    wavenum_sq = wavenumber*wavenumber

    kx = np.fft.fftfreq(u.shape[0],dx/(2*np.pi))
    kz_sq = kx*kx

    mask = wavenumber*wavenumber>kz_sq

    g = np.zeros([len(kx)], dtype=np.complex)
    g[mask] = np.exp(1j*np.sqrt(wavenum_sq-kz_sq[mask])*z)

    res = np.fft.ifft(g*np.fft.fft(u)) 

    return res
   

##################################################################
"""
def angular_spectrum_propagation_2d(array,energy,pixel_size,distance):
    wavelength = energy_to_wavelength(energy)
    du = 1/(pixel_size*np.size(array, 0))
    ii, jj = np.float64(np.meshgrid(range(0, np.size(array, 1)),range(0, np.size(array, 0))))
    ii -= np.size(array,1)/2
    jj -= np.size(array,0)/2
    ii *= du[1]; jj *= du[0]
    ii *= wavelength; jj *= wavelength
    ang_mesh = np.square(ii)+np.square(jj)
    ang_mesh = np.sqrt(1 - ang_mesh)
    ang_mesh = np.exp(1j*(2*np.pi*distance/wavelength)*ang_mesh)
    ang_mesh = np.fft.fftshift(ang_mesh, [-2, -1])
    prop_array = np.fft.ifftn(ang_mesh*np.fft.fftn(array,None,[-2, -1]))
    
    return prop_array

##################################################################

def exact_prop(in_wave,out_wave,L_in,L_out,energy,z):
    wavel = energy_to_wavelength(energy)
    pi = np.pi
    N_in = np.shape(in_wave)[0]
    in_domain = np.linspace(-L_in/2,L_in/2,N_in)
    N_out = np.shape(out_wave)[0]
    out_domain = np.linspace(-L_out/2,L_out/2,N_out)
    step_in = L_in/N_in
    for i in range(N_out):
        for j in range(N_in):
            x = in_domain[j]
            f = in_wave[j]
            x1 = out_domain[i]
            out_wave[i] += f*np.exp((-1j*pi*x*x)/(wavel*z))*np.exp((-1j*2*pi*x*x1)/(wavel*z))
            #out_wave[i] += ne.evaluate('f*exp((-1j*pi*x*x)/(wavel*z))*exp((-1j*2*pi*x*x1)/(wavel*z))')
    out_wave *= (1/np.sqrt(1j*wavel*z))*step_in
    return out_wave

##################################################################

def propTF(u,step,L1,energy,z):
    wavel = energy_to_wavelength(energy)

    N = np.shape(u)[0]
    pi = np.pi
    F = np.fft.fftfreq(N,step)
    u = np.fft.fft(u)
    u = np.exp(-1j*(2*pi*z/wavel)*np.sqrt(1-wavel**2*(F**2)))*u
    u = np.fft.ifft(u)
    return u,L1

##################################################################

"""