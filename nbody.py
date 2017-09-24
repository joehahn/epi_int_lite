#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#to execute in pyspark's ipython shell:
#    PYSPARK_DRIVER_PYTHON=/emr/miniconda2/bin/ipython pyspark --master yarn --num-executors 6 --executor-cores 1 --executor-memory 4G
#to run locally on master node....can have memory issues...:
#    PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[*] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#to submit spark job to yarn:
#74 sec#cache every output# PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master yarn --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#53 sec#cache every output# PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[*] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#48 sec#cache every output# PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[1] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#42 sec#cache every dt    # PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[1] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#42 sec#no cache          # PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[1] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py

#set number of streamlins and particles per streamline
number_of_streamlines = 2
particles_per_streamline = 5
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 3
total_number_of_outputs = 10

#radial width assuming circular orbits
radial_width = 1.0e-3

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.0e-3

#start time
import time
time_start = time.time()

#angular frequency
def Omega(a):
    GM = 1.0
    Omega2 = GM/a/a/a
    Omega = Omega2**0.5
    return Omega
from pyspark.sql.functions import udf
from pyspark.sql.types import *
Omega_udf = udf(Omega, DoubleType())

#epicyclic frequency
def Kappa(a):
    return Omega(a)
Kappa_udf = udf(Kappa, DoubleType())

#adjust angles to live between -Pi and Pi
import math
def adjust_angle(angle):
    twopi = 2.0*math.pi
    if (angle > math.pi):
        angle -= twopi
    if (angle < -math.pi):
        angle += twopi     
    return angle
adjust_angle_udf = udf(adjust_angle, DoubleType())

#drift step advances M and timestep
def drift(df, dt):
    df = df.withColumn('Kappa', Kappa_udf(df['a']))
    M = adjust_angle_udf(df['M'] + df['Kappa']*dt)
    df = df.withColumn('M', M)
    df = df.withColumn('timestep', df['timestep'] + 1)
    return df

#convert orbit elements to coordinates
from pyspark.sql.functions import sin, cos
def elem2coord(df):
    df = df.withColumn('e_sin_M', df['e']*sin(df['M']))
    df = df.withColumn('e_cos_M', df['e']*cos(df['M']))
    df = df.withColumn('r', df['a']*(1 - df['e_cos_M']))
    df = df.withColumn('Omega', Omega_udf(df['a']))
    df = df.withColumn('Kappa', Kappa_udf(df['a']))
    df = df.withColumn('t', (df['Omega']/df['Kappa'])*(df['M'] + 2.0*df['e_sin_M']) + df['wt'])
    df = df.withColumn('vr', df['a']*df['Kappa']*df['e_sin_M'])
    df = df.withColumn('vt', df['a']*df['Omega']*(1.0 + df['e_cos_M']))
    cols = ['id', 'timestep', 'streamline', 'a', 'r', 't', 'vr', 'vt']
    return df.select(cols)

#convert coordinates to elements assuming a remains constant
from pyspark.sql.functions import sqrt, atan2
def coords2elem(df):
    df = df.withColumn('Omega', Omega_udf(df['a']))
    df = df.withColumn('Kappa', Kappa_udf(df['a']))
    df = df.withColumn('e_sin_M', df['vr']/(df['a']*df['Kappa']))
    df = df.withColumn('e_cos_M', 1.0 - df['r']/df['a'])
    e2 = df['e_sin_M']**2 + df['e_cos_M']**2
    df = df.withColumn('e', sqrt(e2))
    df = df.withColumn('M', atan2(df['e_sin_M'], df['e_cos_M']))
    wt = df['t'] - (df['Omega']/df['Kappa'])*(df['M'] + 2.0*df['e_sin_M'])
    df = df.withColumn('wt', adjust_angle_udf(wt))
    cols = ['id', 'timestep', 'streamline', 'a', 'e', 'M', 'wt']
    return df.select(cols)

#create SparkSession
from pyspark.sql import SparkSession
spark = SparkSession.builder.appName('nbody').getOrCreate()

#initialize particles in circular orbits
from pyspark.sql.types import Row
ids = [Row(id=idx) for idx in range(number_of_particles)]
df = spark.createDataFrame(ids)
from pyspark.sql.functions import lit
df = df.withColumn('timestep', lit(0))
df = df.withColumn('streamline', (df['id']/particles_per_streamline).cast(IntegerType()) )
df = df.withColumn('longitude_index', df['id']%particles_per_streamline)
df = df.withColumn('a', lit(1.0))
if (number_of_streamlines > 1):
    df = df.withColumn('a', 1.0 + df['streamline']*radial_width/(number_of_streamlines - 1))
df = df.withColumn('e', lit(0.0))
df = df.withColumn('M', lit(0.0))
df = df.withColumn('wt', math.pi*(2.0*df['longitude_index']/particles_per_streamline - 1.0))
cols = ['id', 'timestep', 'streamline', 'a', 'e', 'M', 'wt']
particles_init = df.select(cols)

#alter initial orbits as needed
if (initial_orbits == 'breathing mode'):
    particles_init = particles_init.withColumn('e', lit(initial_e))
    particles_init = particles_init.withColumn('M', lit(0.0))

#evolve system
print 'evolving system...'
timestep = 0
number_of_outputs = 0
particles = particles_init
while (number_of_outputs < total_number_of_outputs):
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #drift orbit elements and advance timestep
        particles_drifted = drift(particles, dt)
        #update coordinates
        particles_coords = elem2coord(particles_drifted)
        #compute accelerations
        #convert coordinates to elements
        particles_elem = coords2elem(particles_coords)
        #updates
        particles = particles_elem
        timesteps_since_output += 1
    #particles.cache()
    number_of_outputs += 1
    print 'number_of_outputs = ', number_of_outputs
    #print particles.show()

#display results
cols = ['id', 'timestep', 'a', 'e', 'wt', 'M']
particles.show()
time_stop = time.time()
print 'execution time (sec) = ', time_stop - time_start
