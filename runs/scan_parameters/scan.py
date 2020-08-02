#!/usr/bin/env python

#scan.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 31 July 2020.
#
#this launches multiple concurrent epi_int simulations having varied input parameters

#set number of concurrent processes
N_processes = 7

#set nominal_viscosity and timesteps_per_output
nominal_viscosity = 1.0e-11
nominal_timesteps_per_output = 420

#all possible parameter variations
params = {
    'total_ring_mass':[2.0e-10, 2.0e-09, 2.0e-8],
    'radial_width':[0.0001, 0.0002, 0.0005, 0.001],
    'viscosity_multiplier':[0.2, 1.0, 5.0]
}

#execution start time
import time as tm
clock_start = tm.time()

#generate all permutations of the above
keys, values = zip(*params.items())
import itertools
permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]

#update viscosities as needed, with timesteps_per_output scaling as 1/shear_viscosity
if ('viscosity_multiplier' in params.keys()): 
    for p in permutations:
        viscosity_multiplier = p.pop('viscosity_multiplier')
        p['shear_viscosity'] = viscosity_multiplier*nominal_viscosity
        p['bulk_viscosity'] = p['shear_viscosity']
        timesteps_per_output = int(nominal_timesteps_per_output*nominal_viscosity/p['shear_viscosity'])
        if (timesteps_per_output < 1):
            timesteps_per_output = 1
        p['timesteps_per_output'] = timesteps_per_output

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
print ('execution time (minutes) = ', (tm.time() - clock_start)/60.0
