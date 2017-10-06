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
def Kappa(J2, Rp, a):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Kappa2 = (GM/a2/a)*(1.0 - (1.5*J2)*Ra2)
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
    return adjust_angle(M + Kappa(J2, Rp, a)*dt)

#velocity kicks
def kick(lambda0, r, vr, dt):
    #radial acceleration due to streamline gravity
    G = 1.0
    Nr, Nt = vr.shape
    Ar = np.zeros((Nr, Nt))
    for shft in range(1, Nr):
        dr = np.roll(r, -shft, axis=0) - r
        Ar += 2.0*G*lambda0/dr
    vr += Ar*dt
    return vr

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

#stash current a,e,wt,M at the end of list e_store etc
def save_arrays(az, ez, wtz, Mz, timestep, timestepz, a, e, wt, M):
    az.append(a)
    ez.append(e)
    wtz.append(wt)
    Mz.append(M)
    timestepz.append(timestep)
    return az, ez, wtz, Mz

#save orbit element arrays in files
def save_output(a, e, wt, M, times):
    np.save('output/a.npy', a)
    np.save('output/e.npy', e)
    np.save('output/wt.npy', wt)
    np.save('output/M.npy', M)
    np.save('output/times.npy', times)

#restore orbit elements from files
def restore_output():
    a = np.load('output/a.npy')
    e = np.load('output/e.npy')
    wt = np.load('output/wt.npy')
    M = np.load('output/M.npy')
    times = np.load('output/times.npy')
    return a, e, wt, M, times
