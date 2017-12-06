#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 15
particles_per_streamline = 1
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 17000
total_number_of_outputs = 200

#ring radial width assuming circular orbits
radial_width = 3.0e-3

#total ring mass
total_ring_mass = 1.0e-10

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity
G_ring = -1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = 1.0e-11

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'circular'
initial_e = 0.0

#output folder
output_folder = 'output'
