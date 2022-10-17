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

#execution start time
import time as tm
clock_start = tm.time()

#get nominal simulation's input parameters
execfile('inputs.py')

#generate all possible permutations of the values, initially assume all sims' dynamical_timescale=max(times)/20
keys, values = zip(*params.items())
import itertools
permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]
import pandas as pd
df = pd.DataFrame(data=permutations)
df['sim_id'] = df.index
times_max = dt*timesteps_per_output*total_number_of_outputs
dynamical_timescale = times_max/20.0
df['dynamical_timescale'] = dynamical_timescale
df_permutations = df

#get dynamical_timescale from file df_results.parquet if it exists
#unconfined sims with dynamical_timescale<viscous_timescale get dynamical_timescale=viscous_timescale
import os
file = 'df_results.parquet'
if (os.path.exists(file)):
    df = pd.read_parquet(file)
    df = df[['sim_id', 'viscous_timescale', 'dynamical_timescale', 'outcome']]
    idx = (df.outcome == 'unconfined?') & (df.dynamical_timescale < df.viscous_timescale)
    df.loc[idx, 'dynamical_timescale'] = df.loc[idx, 'viscous_timescale']
    df = df[['sim_id', 'dynamical_timescale']]
    idx = (df.dynamical_timescale > 0)
    df = df[idx]
    df = df.rename({'dynamical_timescale':'dynamical_timescale_obs'}, axis=1)
    df_dynamical = df
else:
    df_dynamical = None

#update dynamical_timescale as needed
df = df_permutations
if (df_dynamical is not None): 
    df = df.merge(df_dynamical, on='sim_id', how='left')
    idx = (df.dynamical_timescale_obs > 0)
    df.loc[idx, 'dynamical_timescale'] = df.loc[idx, 'dynamical_timescale_obs']
    df = df.drop('dynamical_timescale_obs', axis=1)
df_update = df

#set timesteps_per_output
df = df_update
#execution_time = 20*df.dynamical_timescale
execution_time = 10*df.dynamical_timescale
df['timesteps_per_output'] = (execution_time/(dt*total_number_of_outputs)).astype(int)
print 'min timesteps_per_output = ', df.timesteps_per_output.min()
df_timesteps = df

#set output_folder, note that bulk_viscosity=shear_viscosity
df = df_timesteps
df['output_folder'] = 'permutations/'
df.output_folder += 'sim_id=' + df.sim_id.astype(str) + '!'
df.output_folder += 'total_ring_mass=' + df.total_ring_mass.astype(str) + '!'
df.output_folder += 'shear_viscosity=' + df.shear_viscosity.astype(str) + '!'
df.output_folder += 'bulk_viscosity=' + df.shear_viscosity.astype(str) + '!'
df.output_folder += 'radial_width=' + df.radial_width.astype(str) + '!'
df.output_folder += 'timesteps_per_output=' + df.timesteps_per_output.astype(str) + '!'
df_output = df

#create the commands that will execute each epi_int job, with bulk_viscosity=shear_viscosity
df = df_output
df['command'] = ''
import json
for idx, row in df.iterrows():
    d = {'sim_id':row.sim_id, 'total_ring_mass':row.total_ring_mass, 'shear_viscosity':row.shear_viscosity,
        'bulk_viscosity':row.shear_viscosity, 'radial_width':row.radial_width, 'timesteps_per_output':row.timesteps_per_output,
        'output_folder':row.output_folder}
    d_json = json.dumps(d)
    command = "./epi_int_lite.py -m '" + d_json + "'"
    df.loc[idx, 'command'] = command
df_command = df

#remove previous scans
cmd = 'rm -rf permutations/*'
r = os.system(cmd)
cmd = 'touch permutations/nothing'
r = os.system(cmd)

#make list of commands that will be executed in parallel
df = df_command
commands = df.command.values.tolist()

#execute simulations in parallel
print '********'
print 'N_processes = ', N_processes
print 'total number of simulations = ', len(commands)
print '********'
from multiprocessing import Pool
pool = Pool(processes=N_processes)
results = pool.map(os.system, commands)
print 'number of simulations executed = ', len(results)
print 'execution time (hours) = ', (tm.time() - clock_start)/3600.0
