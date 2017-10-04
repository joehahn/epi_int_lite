#helper_fns.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 28 September 2017.
#
#helper functions used by nbody.py

#angular frequency
def Omega(a):
    GM = 1.0
    Omega_0 = np.sqrt(GM/a)/a
    return Omega_0

#epicyclic frequency
def Kappa(a):
    return Omega(a)

#adjust angles to live between -Pi and Pi
import numpy as np
def adjust_angle(angle):
    idx = angle > np.pi
    angle[idx] -= 2.0*np.pi
    idx = angle < -np.pi
    angle[idx] += 2.0*np.pi    
    return angle

#drift step advances M
def drift(a, M, dt):
    return adjust_angle(M + Kappa(a)*dt)

#convert orbit elements to coordinates
def elem2coords(a, e, wt, M):
    e_sin_M = e*np.sin(M)
    e_cos_M = e*np.cos(M)
    r = a*(1.0 - e_cos_M)
    Omg = Omega(a)
    Kap = Kappa(a)
    t = adjust_angle(   (Omg/Kap)*(M + 2.0*e_sin_M) + wt   )
    vr = a*Kap*e_sin_M
    vt = a*Omg*(1.0 + e_cos_M)
    return r, t, vr, vt

#convert coordinates to orbit elements
def coords2elem(r, t, vr, vt, a):
    Omg = Omega(a)
    Kap = Kappa(a)
    e_sin_M = vr/(a*Kap)
    e_cos_M = 1.0 - r/a
    e = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
    M = np.arctan2(e_sin_M, e_cos_M)
    wt = adjust_angle(   t - (Omg/Kap)*(M + 2.0*e_sin_M)   )
    return e, wt, M

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


