#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 5
particles_per_streamline = 31
#number_of_streamlines = 3
#particles_per_streamline = 5
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 10
total_number_of_outputs = 650

#ring radial width assuming circular orbits
radial_width = 1.0e-3

#total ring mass
total_ring_mass = 1.0e-7
#total_ring_mass = 0.0

#ring kinematic shear viscosity
shear_viscosity = 0.0e-15

#oblateness parameters
Rp = 0.5
J2 = 0.02

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.5e-3

##choose initial orbits
#initial_orbits = 'eccentric'
#initial_e = 1.5e-3

##choose initial orbits
#initial_orbits = 'circular'
#initial_e = 0.0

#output folder ends in /
output_folder = 'output/'

