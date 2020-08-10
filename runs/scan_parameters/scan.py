#!/usr/bin/env python

#scan.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 31 July 2020.
#
#this launches multiple concurrent epi_int simulations having varied input parameters

#set number of concurrent processes
N_processes = 7

#all possible parameter variations
params = {
    'total_ring_mass':[1.0e-12, 3.3e-12, 1.0e-11, 3.3e-11, 1.0e-10, 3.3e-10, 1.0e-9, 3.3e-9, 1.0e-8, 3.3e-8, 1.0e-7],
    'radial_width':[0.0001, 0.0002, 0.0003, 0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.003],
    'shear_viscosity':[1.0e-11, 1.0e-10, 1.0e-9]
}

#all possible parameter variations
params = {
    'total_ring_mass':[8.0e-11, 4.0e-10, 2.0e-9, 1.0e-8, 5.0e-8],
    'radial_width':[0.0002, 0.0004, 0.0006, 0.0008, 0.0012],
    'shear_viscosity':[1.0e-12, 5.0e-12, 3.0e-11, 2.0e-10, 1.0e-9]
}

#execution start time
import time as tm
clock_start = tm.time()

#calculate nominal nominal_viscous_timescale
execfile('inputs.py')
import numpy as np
nominal_timesteps_per_output = timesteps_per_output
nominal_total_ring_mass = total_ring_mass
nominal_radial_width = radial_width
nominal_shear_viscosity = shear_viscosity

#generate all permutations of the above
keys, values = zip(*params.items())
import itertools
permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]

#adjust timesteps_per_output to scale as 1/shear_viscosity, 1/total_ring_mass
for idx, p in enumerate(permutations):
    total_ring_mass = p['total_ring_mass']
    radial_width = p['radial_width']
    shear_viscosity = p['shear_viscosity']
    factor = total_ring_mass/nominal_total_ring_mass
    factor *= (nominal_shear_viscosity/shear_viscosity)**0.5
    factor *= (nominal_radial_width/radial_width)**0.5
    timesteps_per_output = int(nominal_timesteps_per_output*factor)
    print total_ring_mass, nominal_total_ring_mass, factor
    if (timesteps_per_output < 1):
        timesteps_per_output = 1
    p['timesteps_per_output'] = timesteps_per_output
    p['sim_id'] = idx

#add output_folders to each permutation
for p in permutations:
    output_folder = 'permutations/'
    for k,v in p.iteritems():
        output_folder += k + '=' + str(v) + '_'
    p['output_folder'] = output_folder

#generate list of commands to be executed in parallel
import json
commands = []
for p in permutations:
    p_json = json.dumps(p)
    cmd = "./epi_int_lite.py -m '" + p_json + "'"
    commands += [cmd]

#remove previous scans
import os
cmd = 'rm -rf permutations/*'
r = os.system(cmd)

#execute simulations in parallel
print '********'
print 'N_processes = ', N_processes
print 'total number of simulations = ', len(commands)
print '********'
from multiprocessing import Pool
pool = Pool(processes=N_processes)
results = pool.map(os.system, commands)
print 'number of simulations executed = ', len(results)
print 'execution time (minutes) = ', (tm.time() - clock_start)/60.0
