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

ibrun {{ solver }}.MPI -ng 1 -groupfile none -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }}
