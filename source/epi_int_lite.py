#!/usr/bin/env python

#epi_int_lite.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this simulates the dynamical evolution of self-gravitating planetary rings.
#A portion of this code was written while drinking and singing at
#the Water Tank karaoke bar in northwest Austin TX, so buyer beware.

#read input parameters
import numpy as np
execfile('inputs.py')
print 'number_of_streamlines =', number_of_streamlines
print 'particles_per_streamline =', particles_per_streamline
print 'dt =', dt
print 'timesteps_per_output = ', timesteps_per_output
print 'total_number_of_outputs =', total_number_of_outputs
print 'radial_width =', radial_width
print 'total_ring_mass =', total_ring_mass
print 'ring gravitation constant =', G_ring
print 'fast_gravity =', fast_gravity
print 'shear_viscosity =', shear_viscosity
print 'bulk_viscosity =', bulk_viscosity
print 'confine_edges =', confine_edges
print "Toomre's Q_ring =", Q_ring
print 'Rp =', Rp
print 'J2 =', J2
print 'initial_orbits =', initial_orbits
print 'output_folder =', output_folder

#initialize orbits
from helper_fns import *
r, t, vr, vt, lambda0, c = initialize_streamline(number_of_streamlines, particles_per_streamline,
    radial_width, total_ring_mass, G_ring, fast_gravity, shear_viscosity, bulk_viscosity, confine_edges,
    Q_ring, J2, Rp, initial_orbits)

#prep for main loop
timestep = 0
number_of_outputs = 0
null_coordinate = False
(rz, tz, vrz, vtz, timestepz) = ([r], [t], [vr], [vt], [timestep])
import time as tm
clock_start = tm.time()

#evolve system using Chamber's (1993) 2nd order kick-drift-kick scheme 
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    #kick velocities forwards by timestep +dt/2
    r, t, vr, vt = velocity_kick(J2, Rp, lambda0, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, r, t, vr, vt, dt/2.0, fast_gravity, confine_edges)
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #kick coordinates to account for central body's motion about center of mass
        r, t, vr, vt = coordinate_kick(dt/2, total_ring_mass, r, t, vr, vt)
        #convert coordinates to elements
        a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #convert orbit elements to coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick coordinates to account for central body's motion about center of mass
        r, t, vr, vt = coordinate_kick(dt/2, total_ring_mass, r, t, vr, vt)    
        #kick velocities
        r, t, vr, vt = velocity_kick(J2, Rp, lambda0, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, r, t, vr, vt, dt, fast_gravity, confine_edges)
        #updates
        timestep += 1
        timesteps_since_output += 1
        #alert if null is obtained
        if ((np.isnan(r).any() == True) and (null_coordinate == False)):
            print 'null coordinate at timestep = ', timestep
            null_coordinate = True
    #kick velocities backwards by timestep -dt/2
    r, t, vr, vt = velocity_kick(J2, Rp, lambda0, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, r, t, vr, vt, -dt/2.0, fast_gravity, confine_edges)
    #save output
    rz, tz, vrz, vtz, timestepz = store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, total_ring_mass, timestep)
    run_time_min = (tm.time() - clock_start)/60.0
    number_of_outputs += 1
    eta_min = int((total_number_of_outputs - number_of_outputs)*run_time_min/number_of_outputs)
    if (20*number_of_outputs%total_number_of_outputs == 0):
        print 'time = ' + str(timestep*dt) + \
            '    number of outputs = ' + str(number_of_outputs) + \
            '    number of orbits = ' + str(int(timestep*dt/2.0/np.pi)) + \
            '    eta (minutes) = ', eta_min

#save results
timez = np.array(timestepz)*dt
save_output(rz, tz, vrz, vtz, timez, lambda0, output_folder)
print 'execution time (minutes) = ', (tm.time() - clock_start)/60.0
