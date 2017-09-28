#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#set number of streamlins and particles per streamline
number_of_streamlines = 3
particles_per_streamline = 2
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 1
total_number_of_outputs = 3

#radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.0e-10

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.0e-3

#choose initial orbits
initial_orbits = 'eccentric'

#choose initial orbits
initial_orbits = 'circular'

#start time
import time as tm
time_start = tm.time()

#initialize particles in circular orbits
import numpy as np
a_streamlines = np.linspace(1.0, 1.0 + radial_width, num=number_of_streamlines)
a_list = []
for a_s in a_streamlines:
    a_list.append(np.zeros(particles_per_streamline) + a_s)
a0 = np.array(a_list)
e0 = np.zeros_like(a0)
M0 = np.zeros_like(a0)
wt_streamline = np.linspace(-np.pi, np.pi, num=particles_per_streamline, endpoint=False)
wt_streamline += (wt_streamline[1] - wt_streamline[0])/2.0
wt_list = []
for idx in range(number_of_streamlines):
    wt_list.append(wt_streamline)
wt0 = np.array(wt_list)

#alter initial orbits as needed
if (initial_orbits == 'circular'):
    pass
if (initial_orbits == 'eccentric'):
    pass
if (initial_orbits == 'breathing mode'):
    e0[:] = initial_e
    M0[:] = 0.0

#lambda=streamline mass-per-lenth
mass_per_streamline = total_ring_mass/number_of_streamlines
lambda0 = np.zeros_like(a0) + mass_per_streamline/(2.0*np.pi*a0)
print 'this lambda-check should equal one = ', \
    (lambda0[:,0]*2.0*np.pi*a_streamlines).sum()/total_ring_mass

#prep for main loop
timestep = 0
the_time = 0.0
number_of_outputs = 0
(a, e, wt, M) = (a0, e0, wt0, M0)
(a_store, e_store, wt_store, M_store) = ([a], [e], [wt], [M])
from helper_fns import *

#evolve system
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #mean anomaly advances during drift step
        M = drift(a, M, dt)
        #update coordinates
        r, t, vr, vt = elem2coords(a, e, wt, M)
        #check longitude ordering
        #compute accelerations
        #update a
        #convert coordinates to elements
        e, wt, M = coords2elem(r, t, vr, vt, a)
        #updates
        timestep += 1
        timesteps_since_output += 1
    #save output
    a_store, e_store, wt_store, M_store = save_arrays(a_store, e_store, wt_store, \
        M_store, a, e, wt, M)
    number_of_outputs += 1
    print 'number_of_outputs = ', number_of_outputs
    print 'number of timesteps = ', timestep
    print 'time = ', timestep*dt
    print t

#save results
time_stop = tm.time()
print 'execution time (sec) = ', time_stop - time_start
