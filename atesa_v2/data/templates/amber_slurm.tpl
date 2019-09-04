#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ taskspernode }}
#SBATCH --export=ALL
#SBATCH -t {{ walltime }}

set -x
module load amber
module load python
export PYTHONPATH=/opt/amber/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/opt/amber/lib:$LD_LIBRARY_PATH
export AMBERHOME=/share/apps/compute/amber16.dat

ibrun /opt/amber/16/bin/{{ solver }}.MPI -ng 1 -groupfile none -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }}
