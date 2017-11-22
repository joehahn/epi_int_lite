#helper_fns.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 28 September 2017.
#
#helper functions used by nbody.py

#angular frequency
def Omega(J2, Rp, a, Ar=0.0):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Omega2 = (GM/a2/a)*(   1.0 + (1.5*J2)*Ra2 - Ar*(a2/GM)   )
    return np.sqrt(Omega2)

#epicyclic frequency
def Kappa(J2, Rp, a, Ar=0.0, kappa_squared=False):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Kappa2 = (GM/a2/a)*(   1.0 - (1.5*J2)*Ra2 - Ar*(a2*(3.0/GM))   )
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

#unwrap angle that lives between -Pi and Pi
def unwrap_angle(angle):
    angle_unw = angle.copy()
    twopi = 2.0*np.pi
    for idx in range(1, len(angle)):
        delta = angle_unw[idx] - angle_unw[idx - 1]
        if (delta < -np.pi):
            angle_unw[idx:] += twopi 
        if (delta > np.pi):
            angle_unw[idx:] -= twopi 
    return angle_unw

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

#velocity kicks due to ring gravity and viscosity
def kick(J2, Rp, lambda0, shear_viscosity, r, t, vr, vt, dt): 
    Ar = np.zeros_like(r)
    At = np.zeros_like(r)
    #acceleration due to streamline gravity
    Ar += ring_gravity(lambda0, r)
    #acceleration due to streamline viscosity
    At += ring_viscosity(shear_viscosity, r, vt)
    #kick velocity
    vr += Ar*dt
    vt += At*dt
    return vr, vt

#convert orbit elements to coordinates
def elem2coords(J2, Rp, a, e, wt, M, Ar=0.0, sort_particle_longitudes=True):
    e_sin_M = e*np.sin(M)
    e_cos_M = e*np.cos(M)
    r = a/np.sqrt(1.0 + 2.0*e_cos_M)
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    t = adjust_angle(   (Omg/Kap)*(M + 2.0*e_sin_M) + wt   )
    ra = r/a
    ra3 = ra*ra*ra
    vr = (a*Kap*ra3)*e_sin_M
    vt = a*a*Omg/r
    #sort each streamline's particles by longitude as needed
    if (sort_particle_longitudes):
        r, t, vr, vt = sort_particles(r, t, vr, vt)
    return r, t, vr, vt

#convert coordinates to orbit elements
def coords2elem(J2, Rp, r, t, vr, vt, Ar=0.0):
    GM = 1.0
    h = r*vt
    c = (h*h)/(2.0*GM*Rp)
    a = Rp*(   c + np.sqrt(c*c - 1.5*J2)   )
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    ar = a/r
    ar2 = ar*ar
    ar3 = ar2*ar
    e_cos_M = (ar2 - 1.0)/2.0
    aK = a*Kap
    e_sin_M = vr*(ar3/aK)
    e = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
    M = np.arctan2(e_sin_M, e_cos_M)
    wt = adjust_angle(   t - (Omg/Kap)*(M + 2.0*e_sin_M)   )
    return a, e, wt, M

#order particles in each streamline by their longitudes
def sort_particles(r, t, vr, vt):
    for streamline_idx in range(len(t)):
        longitude_idx = t[streamline_idx].argsort()
        r[streamline_idx] = r[streamline_idx][longitude_idx]
        t[streamline_idx] = t[streamline_idx][longitude_idx]
        vr[streamline_idx] = vr[streamline_idx][longitude_idx]
        vt[streamline_idx] = vt[streamline_idx][longitude_idx]
    return r, t, vr, vt

#append current r,t,vr,vt,a,timestep to lists rz,tz etc
def store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, timestep):
    rz.append(r)
    tz.append(t)
    vrz.append(vr)
    vtz.append(vt)
    timestepz.append(timestep)
    return rz, tz, vrz, vtz, timestepz

#save orbit element arrays in files
def save_output(r, t, vr, vt, times, output_folder):
    import os
    cmd = 'mkdir -p ' + output_folder
    q = os.system(cmd)
    np.save(output_folder + '/r.npy', r)
    np.save(output_folder + '/t.npy', t)
    np.save(output_folder + '/vr.npy', vr)
    np.save(output_folder + '/vt.npy', vt)
    np.save(output_folder + '/times.npy', times)

#restore orbit elements from files
def restore_output(output_folder):
    r = np.load(output_folder + '/r.npy')
    t = np.load(output_folder + '/t.npy')
    vr = np.load(output_folder + '/vr.npy')
    vt = np.load(output_folder + '/vt.npy')
    times = np.load(output_folder + '/times.npy')
    return r, t, vr, vt, times

#initialize numpy arrays
def initialize_orbits(number_of_streamlines, particles_per_streamline, initial_orbits,
    radial_width, total_ring_mass, J2, Rp, initial_e=None):
    
    #initialize particles in circular orbits
    a_streamlines = np.linspace(1.0, 1.0 + radial_width, num=number_of_streamlines)
    a_list = []
    for sma in a_streamlines:
        a_list.append(np.zeros(particles_per_streamline) + sma)
    a = np.array(a_list)
    e = np.zeros_like(a)
    M = np.zeros_like(a)
    
    #longitude of periapse wt
    wt = np.zeros_like(a)
    wt_streamline = np.linspace(-np.pi, np.pi, num=particles_per_streamline, endpoint=False)
    if (particles_per_streamline > 1): 
        wt_streamline += (wt_streamline[1] - wt_streamline[0])/2.0
    else:
        wt_streamline = np.zeros(particles_per_streamline)
    wt_list = []
    for idx in range(number_of_streamlines):
        wt_list.append(wt_streamline)
    wt = np.array(wt_list)
    
    #lambda0=streamline mass-per-lenth
    mass_per_streamline = total_ring_mass/number_of_streamlines
    twopi = 2.0*np.pi
    lambda0 = np.zeros_like(a) + mass_per_streamline/(twopi*a)
    if (total_ring_mass > 0):
        print 'this lambda-check should equal one = ', \
            (lambda0[:,0]*twopi*a_streamlines).sum()/total_ring_mass
    
    #alter initial orbits as needed
    if (initial_orbits == 'circular'):
        pass
    if (initial_orbits == 'eccentric'):
        pass
    if (initial_orbits == 'breathing mode'):
        #e0[:] = initial_e
        #M0[:] = 0.0
        pass
    if (initial_orbits == 'log-e'):
        #initial e is lograthmically distributed between initial_e[0] < e0 < initial_e[1]
        #while M0 and wt0 are randomized between -pi and pi
        e = np.exp(   np.random.uniform(low=np.log(initial_e[0]), high=np.log(initial_e[1]), size=a.shape)   )
        M = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
        wt = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
    
    #convert elements to coordinates
    Ar = ring_gravity(lambda0, a)
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M, Ar=Ar)
    return r, t, vr, vt, lambda0
