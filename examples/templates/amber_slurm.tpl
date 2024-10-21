#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --output="{{ name }}.%j.%N.out"
#SBATCH -p shared
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ taskspernode }}
#SBATCH --time {{ walltime }}
#SBATCH --mem-per-cpu={{ mem }}

set -x

{{ solver }} -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }}
