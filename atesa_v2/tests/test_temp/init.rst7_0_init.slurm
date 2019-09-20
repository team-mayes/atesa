#!/bin/bash
#SBATCH --job-name="init.rst7_0_init"
#SBATCH --output="init.rst7_0_init.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 00:30:00

set -x
module load amber
module load python
export PYTHONPATH=/opt/amber/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/opt/amber/lib:$LD_LIBRARY_PATH
export AMBERHOME=/share/apps/compute/amber16.dat

ibrun /opt/amber/16/bin/sander.MPI -ng 1 -groupfile none -O -i /Users/tburgin/Documents/PycharmProjects/atesa_v2/atesa_v2/data/input_files/aimless_shooting_init_amber.in -o init.rst7_0_init.out -p topology.prmtop -c init.rst7 -r init.rst7_0_init.rst7 -x init.rst7_0_init.nc