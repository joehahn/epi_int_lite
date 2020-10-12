#fluxes.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 6 June 2018
#
#compute viscous angular momentum flux and luminosity through ringlet's central streamline

import numpy as np
from helper_fns import accelerations

#compute angular momentum flux and luminosity
def calculate_angular_momentum_flux(lambda0, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, wt, 
        fast_gravity, confine_edges):
    #select two central streamlines
    total_number_of_outputs, number_of_streamlines, particles_per_streamline = r.shape
    middle_index = number_of_streamlines/2
    rc = r[:, middle_index-1:middle_index+1,:]
    tc = t[:, middle_index-1:middle_index+1,:]
    vrc = vr[:, middle_index-1:middle_index+1,:]
    vtc = vt[:, middle_index-1:middle_index+1,:]
    wtc = wt[:, middle_index-1:middle_index+1,:]
    lambdac = lambda0[middle_index-1:middle_index+1,:]
    #compute angular momentum flux and luminosity
    twopi = 2*np.pi
    angular_momentum_flux_list = []
    angular_momentum_luminosity_list = []
    for t_idx in range(total_number_of_outputs):
        #compute one-sided acceleration=acceleration of middle streamline due to interior neighbor
        r_now = rc[t_idx]
        t_now = tc[t_idx]
        vr_now = vrc[t_idx]
        vt_now = vtc[t_idx]
        Ar, At = accelerations(lambdac, G_ring, shear_viscosity, bulk_viscosity, c, r_now, t_now, 
           vr_now, vt_now, fast_gravity, confine_edges)
        angular_momentum_flux = lambdac*r_now*At
        angular_momentum_flux = angular_momentum_flux[-1]
        angular_momentum_flux_list += [angular_momentum_flux]
        #compute one-sided torque
        delta_ell = twopi*r_now[-1].mean()/particles_per_streamline
        torque = angular_momentum_flux.sum()*delta_ell
        angular_momentum_luminosity_list += [torque]
    #convert lists to arrays
    angular_momentum_flux = np.array(angular_momentum_flux_list)
    angular_momentum_luminosity = np.array(angular_momentum_luminosity_list)
    #select central streamline
    rc = rc[:, -1, :]
    tc = tc[:, -1, :]
    wtc = wtc[:, -1, :]
    return angular_momentum_flux, angular_momentum_luminosity, rc, tc, wtc

#compute energy momentum flux and luminosity
def calculate_energy_flux(lambda0, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, wt, 
        fast_gravity, confine_edges):
    #select two central streamlines
    total_number_of_outputs, number_of_streamlines, particles_per_streamline = r.shape
    middle_index = number_of_streamlines/2
    rc = r[:, middle_index-1:middle_index+1,:]
    tc = t[:, middle_index-1:middle_index+1,:]
    vrc = vr[:, middle_index-1:middle_index+1,:]
    vtc = vt[:, middle_index-1:middle_index+1,:]
    wtc = wt[:, middle_index-1:middle_index+1,:]
    lambdac = lambda0[middle_index-1:middle_index+1,:]
    #compute energy flux and luminosity
    twopi = 2*np.pi
    energy_flux_list = []
    energy_luminosity_list = []
    for t_idx in range(total_number_of_outputs):
        r_now = rc[t_idx]
        t_now = tc[t_idx]
        vr_now = vrc[t_idx]
        vt_now = vtc[t_idx]
        Ar, At = accelerations(lambdac, G_ring, shear_viscosity, bulk_viscosity, c, r_now, t_now, 
            vr_now, vt_now, fast_gravity, confine_edges)
        energy_flux = lambdac*(Ar*vr_now + At*vt_now)
        energy_flux = energy_flux[-1]
        energy_flux_list += [energy_flux]
        #compute one-sided work
        delta_ell = twopi*r_now[-1].mean()/particles_per_streamline
        work = energy_flux.sum()*delta_ell
        energy_luminosity_list += [work]
    #convert lists to arrays
    energy_flux = np.array(energy_flux_list)
    energy_luminosity = np.array(energy_luminosity_list)
    #select central streamline
    rc = rc[:, -1, :]
    tc = tc[:, -1, :]
    wtc = wtc[:, -1, :]
    return energy_flux, energy_luminosity, rc, tc, wtc


