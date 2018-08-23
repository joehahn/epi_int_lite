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
timesteps_per_output = 10000
total_number_of_outputs = 1000

#ring radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.5e-9

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity
G_ring = 1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
#shear_viscosity = 3.0e-13#...da=0.0018 and hair=-0.06 rad at t=10x10^4 orbits
shear_viscosity = 6.0e-13

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits..adeda = eccentricity gradient = a*(de/da)
initial_orbits = {
    'shape':'eccentric',
    'e':5.0e-3,
    'adeda':0.2263
}

#output folder
output_folder = 'output'

