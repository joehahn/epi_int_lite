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
print 'output_folder =', output_folder

#initialize orbits
import numpy as np
from helper_fns import *
a0, e0, M0, wt0, lambda0 = initialize_orbits(number_of_streamlines, particles_per_streamline,
    initial_orbits, initial_e, radial_width, total_ring_mass, J2, Rp)

#boost vt and a? due to Ar

#prep for main loop
timestep = 0
number_of_outputs = 0
(a, e, wt, M) = (a0.copy(), e0.copy(), wt0.copy(), M0.copy())
(az, ez, wtz, Mz, timestepz) = ([a], [e], [wt], [M], [timestep])
        
#evolve system...this follows Chamber's (1993) 2nd order drift-kick scheme but assumes
#the central mass has no significant motion about system's center-of-mass ie the ring is nearly
#circular and there are no point-mass satellites such that Chamber's exp(tau*C/2)=1
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    #convert orbit elements to coordinates
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    #kick velocities and evolve a forwards by +dt/2
    vr, vt, a = kick(lambda0, shear_viscosity, J2, Rp, r, vr, vt, a, dt/2.0)
    #convert coordinates to elements
    e, wt, M = coords2elem(J2, Rp, r, t, vr, vt, a)
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #convert orbit elements to coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick velocities and evolve a
        vr, vt, a = kick(lambda0, shear_viscosity, J2, Rp, r, vr, vt, a, dt)
        #convert coordinates to elements
        e, wt, M = coords2elem(J2, Rp, r, t, vr, vt, a)
        #updates
        timestep += 1
        timesteps_since_output += 1
    #convert orbit elements to coordinates
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    #kick velocities and evolve a backwards by -dt/2
    vr, vt, a = kick(lambda0, shear_viscosity, J2, Rp, r, vr, vt, a, -dt/2.0)
    #convert coordinates to elements
    e, wt, M = coords2elem(J2, Rp, r, t, vr, vt, a)
    #save output
    number_of_outputs += 1
    az, ez, wtz, Mz = save_arrays(az, ez, wtz, Mz, timestep, timestepz, a, e, wt, M)
    if (20*number_of_outputs%total_number_of_outputs == 0):
        print 'time = ' + str(timestep*dt) + \
            '    number of outputs = ' + str(number_of_outputs) + \
            '    number of orbits = ' + str(int(timestep*dt/2.0/np.pi))

#save results
times = np.array(timestepz)*dt
save_output(az, ez, wtz, Mz, times, output_folder)
time_stop = tm.time()
print 'execution time (sec) = ', time_stop - time_start
