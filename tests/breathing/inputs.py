#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 5
particles_per_streamline = 31
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 4
total_number_of_outputs = 500

#ring radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 2.0e-7

#ring's gravitation constant is usually Q_ring=1 but set Q_ring < 0 to turn off ring gravity
G_ring = 1.0

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = -1.0

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.5

#oblateness parameters, set J2=0 to turn off
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.5e-3
initial_q = 0.0

#output folder
output_folder = 'output'
