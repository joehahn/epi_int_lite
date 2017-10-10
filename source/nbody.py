#!/usr/bin/env python

#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#start time
import time as tm
time_start = tm.time()

#read input parameters
execfile('inputs.py')
print 'number_of_streamlines =', number_of_streamlines
print 'particles_per_streamline =', particles_per_streamline
print 'dt =', dt
print 'timesteps_per_output = ', timesteps_per_output
print 'total_number_of_outputs =', total_number_of_outputs
print 'radial_width =', radial_width
print 'total_ring_mass =', total_ring_mass
print 'shear_viscosity =', shear_viscosity
print 'Rp =', Rp
print 'J2 =', J2
print 'initial_orbits =', initial_orbits
print 'initial_e =', initial_e

#initialize orbits
from helper_fns import *
a0, e0, M0, wt0, lambda0 = initialize_orbits(number_of_streamlines, particles_per_streamline,
    initial_orbits, initial_e, radial_width, total_ring_mass)

#need to kick vt to balance Ar

#prep for main loop
timestep = 0
number_of_outputs = 0
(a, e, wt, M) = (a0, e0, wt0, M0)
(az, ez, wtz, Mz, timestepz) = ([a], [e], [wt], [M], [timestep])

#evolve system
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #update coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick velocities
        vr, vt = kick(lambda0, shear_viscosity, r, vr, vt, dt)
        #update a
        #convert coordinates to elements
        e, wt, M = coords2elem(J2, Rp, r, t, vr, vt, a)
        #updates
        timestep += 1
        timesteps_since_output += 1
    #save output
    number_of_outputs += 1
    az, ez, wtz, Mz = save_arrays(az, ez, wtz, Mz, timestep, timestepz, 
        a, e, wt, M)
    print 'number_of_outputs = ', number_of_outputs
    print 'number of timesteps = ', timestep
    print 'time = ', timestep*dt

#save results
times = np.array(timestepz)*dt
save_output(az, ez, wtz, Mz, times)
time_stop = tm.time()
print 'execution time (sec) = ', time_stop - time_start
