#libration.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 March 2018.
#
#these helper functions are used to analyze the motions (libration etc) of a self-confining ringlet

#imports
import numpy as np

#calculate da, de, dwt differences at inner & outer streamline's periapse
def orbit_deltas(times, r, t, vr, vt, a, e, wt, J2, Rp):
    a_inner_list = []
    a_outer_list = []
    e_inner_list = []
    e_outer_list = []
    wt_inner_list = []
    wt_outer_list = []
    for tidx, time_ in enumerate(times):
        r_now = r[tidx]
        r_inner = r_now[0]
        r_outer = r_now[-1]
        pidx_inner = r_inner.argmin()  
        pidx_outer = r_outer.argmin()    
        wt_now = wt[tidx]
        wt_inner_streamline = wt_now[0]
        wt_outer_streamline = wt_now[-1]
        wt_inner = wt_inner_streamline[pidx_inner]
        wt_outer = wt_outer_streamline[pidx_outer]
        e_now = e[tidx]
        e_inner_streamline = e_now[0]
        e_outer_streamline = e_now[-1]
        e_inner = e_inner_streamline.mean()
        e_outer = e_outer_streamline.mean()
        a_now = a[tidx]
        a_inner_streamline = a_now[0]
        a_outer_streamline = a_now[-1]
        a_inner = a_inner_streamline.mean()
        a_outer = a_outer_streamline.mean()
        a_inner_list += [a_inner]
        a_outer_list += [a_outer]
        e_inner_list += [e_inner]
        e_outer_list += [e_outer]
        wt_inner_list += [wt_inner]
        wt_outer_list += [wt_outer]
    a_inner = np.array(a_inner_list)
    a_outer = np.array(a_outer_list)
    e_inner = np.array(e_inner_list)
    e_outer = np.array(e_outer_list)
    wt_inner = np.array(wt_inner_list)
    wt_outer = np.array(wt_outer_list)
    a_mean = (a_outer + a_inner)/2
    da = a_outer - a_inner
    e_mean = (e_outer + e_inner)/2
    de = e_outer - e_inner
    dwt = wt_outer - wt_inner
    idx = (dwt > np.pi)
    dwt[idx] = dwt[idx] - 2*np.pi
    idx = (dwt < -np.pi)
    dwt[idx] = dwt[idx] + 2*np.pi
    return a_inner, a_outer, a_mean, da, e_inner, e_outer, e_mean, de, wt_inner, wt_outer, dwt

#compute e_prime, wt_prime, q, and H(q)
def calculate_Hq(a_mean, e_mean, da, de, dwt):
    e_prime = a_mean*de/da
    wt_prime = e_mean*a_mean*dwt/da
    q2 = e_prime**2 + wt_prime**2
    q = np.sqrt(q2)
    q_factor = np.sqrt(1 - q2)
    H = (1 - q_factor)/q2/q_factor
    return H, q, e_prime, wt_prime

