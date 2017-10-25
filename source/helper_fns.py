#helper_fns.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 28 September 2017.
#
#helper functions used by nbody.py

#angular frequency
def Omega(J2, Rp, a):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Omega2 = (GM/a2/a)*(1.0 + (1.5*J2)*Ra2)
    return np.sqrt(Omega2)

#epicyclic frequency
def Kappa(J2, Rp, a, kappa_squared=False):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Kappa2 = (GM/a2/a)*(1.0 - (1.5*J2)*Ra2)
    if (kappa_squared):
        return Kappa2
    else:
        return np.sqrt(Kappa2)

#adjust angles to live between -Pi and Pi
import numpy as np
def adjust_angle(angle):
    idx = angle > np.pi
    angle[idx] -= 2.0*np.pi
    idx = angle < -np.pi
    angle[idx] += 2.0*np.pi    
    return angle

#drift step advances M
def drift(a, M, J2, Rp, dt):
    return M + Kappa(J2, Rp, a)*dt

#radial acceleration due to ring self-gravity
def ring_gravity(lambda0, r):
    G = 1.0
    two_G_lambda = 2.0*G*lambda0
    Ar = np.zeros_like(r)
    Nr, Nt = r.shape
    for shft in range(1, Nr):
        dr = np.roll(r, -shft, axis=0) - r
        Ar += two_G_lambda/dr
    return Ar

#tangential acceleration due to ring viscosity
def ring_viscosity(shear_viscosity, r, vt):
    factor = 1.5*shear_viscosity*vt/r
    #acceleration that each streamline exerts on exterior neighbor
    dr = np.roll(r, -1, axis=0) - r
    At_ext = factor/dr
    At_ext[dr < 0.0] = 0.0
    At_ext = np.roll(At_ext, 1, axis=0)
    #acceleration that each streamline exerts on interior neighbor
    dr = r - np.roll(r, 1, axis=0)
    At_int = factor/dr
    At_int[dr < 0.0] = 0.0
    At_int = -np.roll(At_ext, -1, axis=0)
    At = At_ext + At_int
    return At

#compute semimajor drift due to orbit-averaged torque
def orbit_averaged_da(At, a, J2, Rp, dt):
    Omg = Omega(J2, Rp, a)
    Kap2 = Kappa(J2, Rp, a, kappa_squared=True)
    for At_streamline in At:
        At_streamline[:] = At_streamline.mean()
    da = (Omg/Kap2)*(2.0*dt)*At
    return da
    
#velocity kicks
def kick(lambda0, shear_viscosity, J2, Rp, r, vr, vt, a, dt):
    Ar = np.zeros_like(r)
    At = np.zeros_like(r)
    #acceleration due to streamline gravity
    Ar += ring_gravity(lambda0, r)
    #acceleration due to streamline viscosity
    At += ring_viscosity(shear_viscosity, r, vt)
    #kick velocity
    vr += Ar*dt
    vt += At*dt
    #orbit-averaged kick to a
    a += orbit_averaged_da(At, a, J2, Rp, dt)
    return vr, vt, a

#convert orbit elements to coordinates
def elem2coords(J2, Rp, a, e, wt, M, sort_particle_longitudes=True):
    e_sin_M = e*np.sin(M)
    e_cos_M = e*np.cos(M)
    r = a*(1.0 - e_cos_M)
    Omg = Omega(J2, Rp, a)
    Kap = Kappa(J2, Rp, a)
    t = adjust_angle(   (Omg/Kap)*(M + 2.0*e_sin_M) + wt   )
    vr = a*Kap*e_sin_M
    vt = a*Omg*(1.0 + e_cos_M)
    #sort each streamline's particles by longitude as needed
    if (sort_particle_longitudes):
        r, t, vr, vt = sort_particles(r, t, vr, vt)
    return r, t, vr, vt

#convert coordinates to orbit elements
def coords2elem(J2, Rp, r, t, vr, vt, a):
    Omg = Omega(J2, Rp, a)
    Kap = Kappa(J2, Rp, a)
    e_sin_M = vr/(a*Kap)
    e_cos_M = 1.0 - r/a
    e = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
    M = np.arctan2(e_sin_M, e_cos_M)
    wt = adjust_angle(   t - (Omg/Kap)*(M + 2.0*e_sin_M)   )
    return e, wt, M

#order particles in each streamline by their longitudes
def sort_particles(r, t, vr, vt):
    for streamline_idx in range(len(t)):
        longitude_idx = t[streamline_idx].argsort()
        r[streamline_idx] = r[streamline_idx][longitude_idx]
        t[streamline_idx] = t[streamline_idx][longitude_idx]
        vr[streamline_idx] = vr[streamline_idx][longitude_idx]
        vt[streamline_idx] = vt[streamline_idx][longitude_idx]
    return r, t, vr, vt

#append current a,e,wt,M to the end of lists az,ez etc
def save_arrays(az, ez, wtz, Mz, timestep, timestepz, a, e, wt, M):
    az.append(a)
    ez.append(e)
    wtz.append(wt)
    Mz.append(M)
    timestepz.append(timestep)
    return az, ez, wtz, Mz

#save orbit element arrays in files
def save_output(a, e, wt, M, times, output_folder):
    import os
    cmd = 'mkdir -p ' + output_folder
    r = os.system(cmd)
    np.save(output_folder + '/a.npy', a)
    np.save(output_folder + '/e.npy', e)
    np.save(output_folder + '/wt.npy', wt)
    np.save(output_folder + '/M.npy', M)
    np.save(output_folder + '/times.npy', times)

#restore orbit elements from files
def restore_output(output_folder):
    a = np.load(output_folder + '/a.npy')
    e = np.load(output_folder + '/e.npy')
    wt = np.load(output_folder + '/wt.npy')
    M = np.load(output_folder + '/M.npy')
    times = np.load(output_folder + '/times.npy')
    return a, e, wt, M, times

#initialize numpy arrays
def initialize_orbits(number_of_streamlines, particles_per_streamline, initial_orbits,
    initial_e, radial_width, total_ring_mass):
    
    #initialize particles in circular orbits
    a_streamlines = np.linspace(1.0, 1.0 + radial_width, num=number_of_streamlines)
    a_list = []
    for a_s in a_streamlines:
        a_list.append(np.zeros(particles_per_streamline) + a_s)
    a0 = np.array(a_list)
    e0 = np.zeros_like(a0)
    M0 = np.zeros_like(a0)
    wt_streamline = np.linspace(-np.pi, np.pi, num=particles_per_streamline, endpoint=False)
    if (particles_per_streamline > 1): 
        wt_streamline += (wt_streamline[1] - wt_streamline[0])/2.0
    else:
        wt_streamline = np.zeros(particles_per_streamline)
    wt_list = []
    for idx in range(number_of_streamlines):
        wt_list.append(wt_streamline)
    wt0 = np.array(wt_list)

    #alter initial orbits as needed
    if (initial_orbits == 'eccentric'):
        pass
    if (initial_orbits == 'breathing mode'):
        e0[:] = initial_e
        M0[:] = 0.0
    if (initial_orbits == 'random'):
        e0 = np.exp(   np.random.uniform(low=np.log(initial_e[0]), high=np.log(initial_e[1]), size=e0.shape)   )
        M0 = np.random.uniform(low=-np.pi, high=np.pi, size=M0.shape)
        #wt0 = np.random.uniform(low=-np.pi, high=np.pi, size=wt0.shape)
        wt0 = np.zeros_like(a0)
    
    #lambda0=streamline mass-per-lenth
    mass_per_streamline = total_ring_mass/number_of_streamlines
    twopi = 2.0*np.pi
    lambda0 = np.zeros_like(a0) + mass_per_streamline/(twopi*a0)
    if (total_ring_mass > 0):
        print 'this lambda-check should equal one = ', \
            (lambda0[:,0]*twopi*a_streamlines).sum()/total_ring_mass
    
    return a0, e0, M0, wt0, lambda0
