#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline
number_of_streamlines = 130#130
particles_per_streamline = 50#50
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 560#2000
total_number_of_outputs = 100

#ring radial width assuming circular orbits
radial_width = 1.1e-3

#total ring mass
total_ring_mass = 3.454196830585073e-10

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=False since there is very little speed benefit when the fast_gravity approximation is used
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = 1.0e-13
bulk_viscosity = 1.0*shear_viscosity

#add fictitious torques at inner and/or outer streamlines, to oppose any radial spreading
confine_inner_edge = True
confine_outer_edge = True

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = 2.0

#oblateness parameters
Rp = 0.5120761875660467
J2 = 0.016932376630835

#perturbations from satellite's m^th Lindblad resonance, set mass_final < 0 to turn off
satellite = {
    'mass_final' : 6.5969e-8,
    'time_grow' : 5.0e3,
    'mass': 0.0,
    'r' : 1.578325914739137,
    't' : 0.0,
    'm' : 2,
    'lc' : 1.0,
    'dlc' : 1.5
}

#choose ringlet's initial orbits
initial_orbits = {
    'shape':'circular',
}

#output folder
output_folder = 'output'

#parse any optional commandline modifications to the above, note these modifications will be ignored by Jupyter
import argparse
import json
parser = argparse.ArgumentParser() 
parser.add_argument('-m', '--modified_params', type=str, dest='modified_params', required=False)
args, unknown_args = parser.parse_known_args()
modified_params = None
if (args.modified_params):
    modified_params_str = args.modified_params
    modified_params = json.loads(modified_params_str)
    print 'modified_params =', modified_params
    for key, val in modified_params.iteritems():
        exec(key + '=val')
