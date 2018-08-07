#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 2
particles_per_streamline = 101
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 1000
total_number_of_outputs = 500

#ring radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.5e-10

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity
G_ring = 1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = -1.0e-11

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits..adeda = eccentricity gradient = a*(de/da)
#    'adeda':0.897         #I =  0.00110214925612 at fixed point
#    'adeda':0.900         #I =  0.00479292151221
#    'adeda':0.910         #I =  0.0229534093052
#    'adeda':0.923         #I =  0.0481953480501
#    'adeda':0.928         #I =  0.0584448849972
initial_orbits = {
    'shape':'eccentric',
    'e':5.0e-3,
    'adeda':0.928
}

#output folder
output_folder = 'output'