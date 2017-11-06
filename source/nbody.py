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
r0, t0, vr0, vt0, a0, lambda0 = initialize_orbits(number_of_streamlines, particles_per_streamline,
    initial_orbits, initial_e, radial_width, total_ring_mass, J2, Rp)

#prep for main loop
timestep = 0
number_of_outputs = 0
(r, t, vr, vt) = (r0.copy(), t0.copy(), vr0.copy(), vt0.copy())
(rz, tz, vrz, vtz, timestepz) = ([r], [t], [vr], [vt], [timestep])

#evolve system...this largely follows Chamber's (1993) 2nd order drift-kick scheme but assumes
#the central mass has negligable motion about system's center-of-mass ie the ring is nearly
#circular and there are no point-mass satellites such that Chamber's exp(tau*C/2)=1
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    #kick velocities forwards by timestep +dt/2
    vr, vt = kick(J2, Rp, lambda0, shear_viscosity, r, t, vr, vt, dt/2.0)
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #convert coordinates to elements
        a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #convert orbit elements to coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick velocities
        vr, vt = kick(J2, Rp, lambda0, shear_viscosity, r, t, vr, vt, dt)
        #updates
        timestep += 1
        timesteps_since_output += 1
    #kick velocities backwards by timestep -dt/2
    vr, vt = kick(J2, Rp, lambda0, shear_viscosity, r, t, vr, vt, -dt/2.0)
    #save output
    number_of_outputs += 1
    rz, tz, vrz, vtz, timestepz = store_system(rz, tz, vrz, vtz, timestepz, 
        r, t, vr, vt, timestep)
    if (20*number_of_outputs%total_number_of_outputs == 0):
        print 'time = ' + str(timestep*dt) + \
            '    number of outputs = ' + str(number_of_outputs) + \
            '    number of orbits = ' + str(int(timestep*dt/2.0/np.pi))

#save results
timez = np.array(timestepz)*dt
save_output(rz, tz, vrz, vtz, timez, output_folder)
time_stop = tm.time()
print 'execution time (sec) = ', time_stop - time_start
