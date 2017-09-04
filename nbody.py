#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#to execute in pyspark's ipython shell:
#    PYSPARK_DRIVER_PYTHON=/emr/miniconda2/bin/ipython pyspark
#to run locally on master node....can have memory issues...:
#    PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[*] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#to submit spark job to yarn:
#    PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master yarn --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py

#set number of streamlins and particles per streamline
number_of_streamlines = 1
particles_per_streamline = 5
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 10
total_number_of_outputs = 11

#radial width assuming circular orbits
radial_width = 1.0e-3

#choose initial orbits
initial_orbits = 'breathing mode'
initial_e = 1.0e-3

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

#drift step
def drift(particles, dt):
    particles = particles.withColumn('M', particles['M'] + particles['Kappa']*dt)
    particles = particles.withColumn('M', adjust_angle_udf(particles['M']))
    return particles

#convert orbit elements to coordinates
from pyspark.sql.functions import sin, cos
def elem2coord(particles):
    particles = particles.withColumn('e_sin_M', particles['e']*sin(particles['M']))
    particles = particles.withColumn('e_cos_M', particles['e']*cos(particles['M']))
    particles = particles.withColumn('r', particles['a']*(1 - particles['e_cos_M']))
    particles = particles.withColumn('t', (particles['Omega']/particles['Kappa'])*\
        (particles['M'] + 2.0*particles['e_sin_M']) + particles['wt'])
    particles = particles.withColumn('vr', particles['a']*particles['Kappa']*particles['e_sin_M'])
    particles = particles.withColumn('vt', particles['a']*particles['Omega']*(1.0 + particles['e_cos_M']))
    return particles

#create SparkSession
from pyspark.sql import SparkSession
spark = SparkSession.builder.appName('nbody').getOrCreate()

#initialize particles in circular orbits
from pyspark.sql.types import Row
particle_ids = [Row(particle_id=idx) for idx in range(number_of_particles)]
particles = spark.createDataFrame(particle_ids)
particles = particles.withColumn('streamline', (particles['particle_id']/particles_per_streamline).cast(IntegerType()) )
particles = particles.withColumn('longitude_index', particles['particle_id']%particles_per_streamline)
from pyspark.sql.functions import lit
particles = particles.withColumn('a', lit(1.0))
if (number_of_streamlines > 1):
    particles = particles.withColumn('a', 1.0 + particles['streamline']*radial_width/(number_of_streamlines - 1))
particles = particles.withColumn('e', lit(0.0))
particles = particles.withColumn('M', lit(0.0))
particles = particles.withColumn('wt', math.pi*(2.0*particles['longitude_index']/particles_per_streamline - 1.0))
particles = particles.withColumn('Omega', Omega_udf(particles['a']))
particles = particles.withColumn('Kappa', Omega_udf(particles['a']))

#alter initial orbits as needed
if (initial_orbits == 'breathing mode'):
    particles = particles.withColumn('e', lit(initial_e))
    particles = particles.withColumn('M', lit(initial_e))

#evolve system
print 'evolving system...'
t = 0.0
number_of_outputs = 0
while (number_of_outputs < total_number_of_outputs):
    number_of_timesteps = 0
    while (number_of_timesteps < timesteps_per_output):
        particles = drift(particles, dt)
        particles = elem2coord(particles)
        t += dt
        number_of_timesteps += 1
    number_of_outputs += 1
    particles.cache()
    print 't, number_of_outputs = ', t, number_of_outputs
particles.show()
