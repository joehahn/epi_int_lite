#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 2
particles_per_streamline = 31
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 12000
total_number_of_outputs = 1000

#ring radial width assuming circular orbits
radial_width = 0.0005

#total ring mass
total_ring_mass = 1.5e-9

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity
G_ring = 1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = 1.0e-12

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits..adeda = eccentricity gradient = a*(de/da)
#    'adeda':0.00345       #I = 0.00410803185633 at fixed point
#    'adeda':0.050         #I = 0.0617254674082
#    'adeda':0.059         #I = 0.0711019702404
#    'adeda':0.062         #I = 0.0741398534508
initial_orbits = {
    'shape':'eccentric',
    'e':5.0e-3,
    'e_prime':0.05657,
    'w_prime':-0.005
}

#output folder
output_folder = 'output'
