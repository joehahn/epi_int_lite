#nbody.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this...

#to execute in pyspark's ipython shell:
#    PYSPARK_DRIVER_PYTHON=/emr/miniconda2/bin/ipython pyspark
#to run locally on master node:
#    PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master local[*] --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py
#to submit spark job to yarn:
#    PYSPARK_PYTHON=/emr/miniconda2/bin/python spark-submit --master yarn --conf "spark.driver.extraJavaOptions=-Dlog4j.configuration=file:./log4j.properties" nbody.py

#set number of streamlins and particles per streamline
number_of_streamlines = 3
particles_per_streamline = 5
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.1
timesteps_per_output = 10
total_number_of_outputs = 3

#create SparkSession
from pyspark.sql import SparkSession
spark = SparkSession.builder.appName('nbody').getOrCreate()

#initialize
from pyspark.sql.types import Row
ids = [Row(id=idx) for idx in range(number_of_particles)]
coords = spark.createDataFrame(ids)
coords.show()

#
t = 0.0
number_of_outputs = 0
while (number_of_outputs < total_number_of_outputs):
    number_of_timesteps = 0
    while (number_of_timesteps < timesteps_per_output):
        t += dt
        number_of_timesteps += 1
        print number_of_outputs, number_of_timesteps, t
    number_of_outputs += 1
