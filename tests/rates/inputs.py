#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 11
particles_per_streamline = 1
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2
timesteps_per_output = 1
total_number_of_outputs = 3200

#ring radial width assuming circular orbits
radial_width = 0.1

#total ring mass
total_ring_mass = 0.0

#ring kinematic shear viscosity
shear_viscosity = 0.0

#oblateness parameters
Rp = 0.5
J2 = 0.0

#choose initial orbits
initial_orbits = 'log-e'
initial_e = (1.0e-6, 1.0e-2)

#output folder
output_folder = 'output'
