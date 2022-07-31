#!/usr/bin/env python

#scan.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 31 July 2020.
#
#this launches multiple concurrent epi_int simulations having varied input parameters, executes in 36 hours 

#set number of concurrent processes
N_processes = 7

#generate range of logarithmically-spaced values for total_ring_mass
import numpy as np
mass_min = 3.33e-13
mass_max = 1.0e-10
N_masses = 6#11
total_ring_mass = np.exp(np.linspace(np.log(mass_min), np.log(mass_max), num=N_masses))
print 'total_ring_mass = ', total_ring_mass.tolist()

#values for radial_width
radial_width = np.array([0.000025, 0.00005, 0.0001, 0.0002])
print 'radial_width = ', radial_width.tolist()

#generate range of logarithmically-spaced values for shear_viscosity
viscosity_min = 3.3e-15
viscosity_max = 3.3e-12
N_viscosities = 7#19
shear_viscosity = np.exp(np.linspace(np.log(viscosity_min), np.log(viscosity_max), num=N_viscosities))
print 'shear_viscosity = ', shear_viscosity.tolist()

#gather the parameter ranges in the params dict
params = {
    'total_ring_mass':total_ring_mass,
    'radial_width':radial_width,
    'shear_viscosity':shear_viscosity
}

#set power laws employed below
mass_power_law = 1.0
viscosity_power_law = -1.0
width_power_law = 0.0

#execution start time
import time as tm
clock_start = tm.time()

#get nominal input parameters
execfile('inputs.py')
nominal_timesteps_per_output = timesteps_per_output
nominal_total_ring_mass = 3.3e-11
nominal_radial_width = 0.00005
nominal_shear_viscosity = 1.0e-13

#generate all possible permutations of the values in params
keys, values = zip(*params.items())
import itertools
permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]

#adjust timesteps_per_output to scale as total_ring_mass^0.5, shear_viscosity^(-0.5), and radial_width^(0.0),
#such that execution_time > 10*viscous_timescale. Also set bulk_viscosity=shear_viscosity
for sim_id, p in enumerate(permutations):
    total_ring_mass = p['total_ring_mass']
    radial_width = p['radial_width']
    shear_viscosity = p['shear_viscosity']
    factor = (total_ring_mass/nominal_total_ring_mass)**mass_power_law
    factor *= (shear_viscosity/nominal_shear_viscosity)**viscosity_power_law
    factor *= (radial_width/nominal_radial_width)**width_power_law
    timesteps_per_output = int(nominal_timesteps_per_output*factor)
    if (timesteps_per_output < 1):
        timesteps_per_output = 1
    execution_time = dt*timesteps_per_output*total_number_of_outputs
    ten_viscous_timescales = 10*(radial_width**2)/(12*shear_viscosity)
    if (execution_time < ten_viscous_timescales):
        print 'resetting timesteps_per_output from ', timesteps_per_output, 
        timesteps_per_output = int(timesteps_per_output*ten_viscous_timescales/execution_time)
        print 'to ', timesteps_per_output
    p['timesteps_per_output'] = timesteps_per_output
    p['sim_id'] = sim_id
    p['bulk_viscosity'] = shear_viscosity
    print sim_id, total_ring_mass, radial_width, shear_viscosity, timesteps_per_output

#add output_folders to each permutation
for p in permutations:
    output_folder = 'permutations/'
    for k,v in p.iteritems():
        output_folder += k + '=' + str(v) + '!'
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
cmd = 'touch permutations/nothing'
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
