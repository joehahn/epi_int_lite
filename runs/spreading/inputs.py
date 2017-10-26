#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 6
particles_per_streamline = 3
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 3500
total_number_of_outputs = 100

#ring radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.0e-10

#ring kinematic shear viscosity
shear_viscosity = 1.0e-11

#oblateness parameters
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'circular'
initial_e = 0.0

#output folder
output_folder = 'output'
