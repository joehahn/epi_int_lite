#!/bin/bash
for entry in $(ls -d *_*/)
do
	echo $entry
	/bin/cp animate.py $entry
	/bin/cp epi_int_lite.py $entry
	/bin/cp submit.csh $entry
	/bin/cp helper_fns.py $entry
	/bin/cp libration.py $entry
	cd $entry
	sbatch -A astronomy-hi submit.csh
	cd -
done

