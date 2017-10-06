#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#set number of streamlins and particles per streamline
number_of_streamlines = 5
particles_per_streamline = 31
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 10
total_number_of_outputs = 650

#radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.0e-7
#total_ring_mass = 0.0

#oblateness parameters
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.5e-3

##choose initial orbits
#initial_orbits = 'eccentric'

##choose initial orbits
#initial_orbits = 'circular'

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

#lambda0=streamline mass-per-lenth
mass_per_streamline = total_ring_mass/number_of_streamlines
twopi = 2.0*np.pi
lambda0 = np.zeros_like(a0) + mass_per_streamline/(twopi*a0)
if (total_ring_mass > 0):
    print 'this lambda-check should equal one = ', \
        (lambda0[:,0]*twopi*a_streamlines).sum()/total_ring_mass

#prep for main loop
timestep = 0
number_of_outputs = 0
(a, e, wt, M) = (a0, e0, wt0, M0)
(az, ez, wtz, Mz, timestepz) = ([a], [e], [wt], [M], [timestep])
from helper_fns import *

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
        vr = kick(lambda0, r, vr, dt)
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

#restore saved data & compare
ar, er, wtr, Mr, timesr = restore_output()
rz, tz, vrz, vtz = elem2coords(J2, Rp, ar, er, wtr, Mr, sort_particle_longitudes=False)

#
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

#pad array
def pad_array(t, longitudes=False):
    Nr, Nt = t.shape
    tp = np.zeros((Nr, Nt+2))
    tp[:, 1:-1] = t
    if (longitudes == True):
        offset = 2.0*np.pi
    else:
        offset = 0.0
    tp[:, 0] = t[:, -1] - offset
    tp[:, -1] = t[:, 0] + offset
    return tp

#this function returns tuple of plot's xy=(x[i],y[i]) coordinates
def xyt(i):
    a = ar[i]
    e = er[i]
    wt = wtr[i]
    M = Mr[i]
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    tp = pad_array(t, longitudes=True)
    rp = pad_array(r, longitudes=False)
    x = tp/np.pi
    y = rp - 1.0
    y_mid = 0*y[len(y)/2].copy()
    for ys in y:
        ys -= y_mid
    tm = timesr[i]
    return (x, y, tm)

#this iterator provides the animation's xyt coordinates
def update():
    for idx in range(len(ar)):
        yield xyt(idx)

#draw frame
def draw(xyt):
    x, y, tm = xyt
    ax.set_title('t = ' + str(tm))
    for idx in range(len(x)):
        line = lines[idx]
        line.set_data(x[idx], y[idx])
    return lines

#show animation
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-0.002, 0.003), 
    xlabel='$\\theta/\pi$', ylabel='$(r - r_o)/r_o$', title='t = 0.0')
x, y, tm = xyt(0)
ax.set_title('t = ' + str(tm))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
lines = [ax.plot([],[], 'o-', markersize=3, color=colors[idx], linewidth=1)[0]
    for idx in range(number_of_streamlines)]
for line in lines:
    line.set_data([],[])
ani = animation.FuncAnimation(fig, draw, update, interval=1, blit=False, repeat=False)
plt.show()
