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
timesteps_per_output = 2500
total_number_of_outputs = 640

#ring radial width assuming circular orbits
radial_width = 0.1

#total ring mass
total_ring_mass = 0.0

#ring's gravitation constant is usually Q_ring=1 but set Q_ring < 0 to turn off ring gravity
G_ring = -1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = -1.0

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters, set J2=0 to turn off
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'log-e'
initial_e = (2.0e-8, 2.0e-2)

#output folder
output_folder = 'output'
