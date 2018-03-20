#level_curves.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 March 2018.
#
#these helper functions are used to analyze a librating ringlet

#imports
import numpy as np

#calculate da, de, dwt differences at inner & outer streamline's periapse
def calculate_deltas(r, a, e, wt):
    dwt_list = []
    de_list = []
    da_list = []
    a_avg_list = []
    e_avg_list = []
    N_times = r.shape[0]
    for tidx in range(N_times):
        r0 = r[tidx]
        r_inner = r0[0]
        r_outer = r0[-1]
        pidx_inner = r_inner.argmin()  
        pidx_outer = r_outer.argmin()    
        wt0 = wt[tidx]
        wt_inner = wt0[0]
        wt_outer = wt0[-1]
        dwt = wt_outer[pidx_outer] - wt_inner[pidx_inner]
        if (dwt > np.pi):
            dwt -= 2*np.pi
        if (dwt < -np.pi):
            dwt += 2*np.pi
        dwt_list += [dwt]
        e0 = e[tidx]
        e_inner = e0[0]
        e_outer = e0[-1]
        de = e_outer[pidx_outer] - e_inner[pidx_inner]
        de_list += [de]
        e_avg = (e_inner[pidx_inner] + e_outer[pidx_outer])/2
        e_avg_list += [e_avg]
        a0 = a[tidx]
        a_inner = a0[0]
        a_outer = a0[-1]
        da = a_outer[pidx_outer] - a_inner[pidx_inner]
        da_list += [da]
        a_avg = (a_inner[pidx_inner] + a_outer[pidx_outer])/2
        a_avg_list += [a_avg]
    dwt = np.array(dwt_list)
    de = np.array(de_list)
    da = np.array(da_list)
    a_avg = np.array(a_avg_list)
    e_avg = np.array(e_avg_list)
    return da, de, dwt, a_avg, e_avg

#compute q2 to lowest order and H(q2)
def H_q2(a_avg, e_avg, da, de, dwt):
    q2 = (a_avg*de/da)**2 + (a_avg*e_avg*dwt/da)**2
    q_factor = np.sqrt(1 - q2)
    H = (1 - q_factor)/q2/q_factor
    return H, q2

#unroll angle
def unroll_angle(angle):
    unrolled_angle = angle.copy()
    for idx in range(1, len(unrolled_angle)):
        delta = unrolled_angle[idx] - unrolled_angle[idx - 1]
        if (delta > np.pi):
            unrolled_angle[idx:] -= 2*np.pi
        if (delta < -np.pi):
            unrolled_angle[idx:] += 2*np.pi
    return unrolled_angle
