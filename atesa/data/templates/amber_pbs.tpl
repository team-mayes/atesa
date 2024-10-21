### PBS preamble

#PBS -N {{ name }}
#PBS -M tburgin@umich.edu
#PBS -m ae

#PBS -A hbmayes_fluxod
#PBS -l qos=flux
#PBS -q fluxod

#PBS -l nodes={{ nodes }}:ppn={{ taskspernode }},mem={{ mem }}
#PBS -l walltime={{ walltime }}
#PBS -j oe
#PBS -V

### End PBS preamble

if [ -s "$PBS_NODEFILE" ] ; then
    echo "Running on"
    cat $PBS_NODEFILE
fi

if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
    echo "Running from $PBS_O_WORKDIR"
fi

### Start job commands below this line

module load hdf5/1.8.16/intel/17.0.1
module load intel/17.0.1 openmpi/1.10.2/intel/17.0.1 cuda/8.0.44 netcdf/4.4.1/intel/17.0.1 PnetCDF/1.9.0/intel/17.0.1
module load amber/18

cd {{ working_directory }}

mpirun {{ solver }}.MPI -ng 1 -groupfile none -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }}
