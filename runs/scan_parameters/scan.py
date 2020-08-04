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
    'total_ring_mass':[5.0e-10, 2.0e-09, 8.0e-9],
    'radial_width':[0.0003, 0.0005, 0.0007],
    'shear_viscosity':[1.0e-12, 3.3e-12, 1.0e-11, 3.3e-11, 1.0e-10]
}

#execution start time
import time as tm
clock_start = tm.time()

#calculate nominal nominal_timesteps_per_output and viscous_timescale
execfile('inputs.py')
nominal_timesteps_per_output = timesteps_per_output
import numpy as np
nominal_viscous_timescale = (radial_width**2)/(12*np.abs(shear_viscosity))

#generate all permutations of the above
keys, values = zip(*params.items())
import itertools
permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]

#adjust timesteps_per_output to scale with viscous_timescale
for idx, p in enumerate(permutations):
    radial_width = p['radial_width']
    shear_viscosity = p['shear_viscosity']
    viscous_timescale = (radial_width**2)/(12*np.abs(shear_viscosity))
    timesteps_per_output = int(nominal_timesteps_per_output*viscous_timescale/nominal_viscous_timescale)
    if (timesteps_per_output < 1):
        timesteps_per_output = 1
    p['timesteps_per_output'] = timesteps_per_output
    p['sim_id'] = idx

#manually tweak selected sims' timesteps_per_output
ids = [40, 21, 44, 36, 18, 25]
for idx, p in enumerate(permutations):
    if (p['sim_id'] in ids):
        p['timesteps_per_output'] *= 2
ids = [43, 24, 39, 27, 42, 12]
for idx, p in enumerate(permutations):
    if (p['sim_id'] in ids):
        p['timesteps_per_output'] *= 4
ids = [4]
for idx, p in enumerate(permutations):
    if (p['sim_id'] in ids):
        p['timesteps_per_output'] //= 2

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
