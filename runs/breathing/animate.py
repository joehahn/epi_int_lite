#!/usr/bin/env python

#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 7 October 2017.
#
#this animates the output of nbody.py

#restore saved data & compare
from helper_fns import *
execfile('inputs.py')
ar, er, wtr, Mr, timesr = restore_output(output_folder)
rz, tz, vrz, vtz = elem2coords(J2, Rp, ar, er, wtr, Mr, sort_particle_longitudes=False)

#get plotting packages
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
    xlabel='longitude   $\\theta/\pi$', ylabel='radius   $(r - r_o)/r_o$', title='t = 0.0')
x, y, tm = xyt(0)
ax.set_title('t = ' + str(tm))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
lines = [ax.plot([],[], 'o-', markersize=3, color=colors[idx], linewidth=1)[0]
    for idx in range(number_of_streamlines)]
for line in lines:
    line.set_data([],[])
ani = animation.FuncAnimation(fig, draw, update, interval=1, blit=False, repeat=False)
plt.show()
