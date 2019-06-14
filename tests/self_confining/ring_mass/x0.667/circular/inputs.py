#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 2
particles_per_streamline = 61
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.2*5
timesteps_per_output = 36000/5
total_number_of_outputs = 1000

#ring radial width assuming circular orbits
radial_width = 0.0005

#total ring mass
total_ring_mass = 1.5e-09

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=True for approximate gravity that is 2x faster and almost as accurate.
G_ring = 1.0
fast_gravity = True

#ring kinematic shear viscosity, set shear_viscosity < 0 to turn off
shear_viscosity = 3.0e-13

#ring kinematic bulk viscosity, set bulk_viscosity < 0 to turn off
bulk_viscosity = shear_viscosity

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits
initial_orbits = {
    'shape':'circular'
}

#output folder
output_folder = 'output'
