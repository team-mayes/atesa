#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ taskspernode }}
#SBATCH --export=ALL
#SBATCH -t {{ walltime }}

set -x

srun $CP2KBIN -i {{ inp }} -o {{ out }}

mv *{{ nc }}* {{ nc }}
mv *{{ rst }}* {{ rst }}