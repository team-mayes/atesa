mpirun {{ solver }}.MPI -ng 1 -groupfile none -O -i {{ inp }} -o {{ out }} -p {{ prmtop }} -c {{ inpcrd }} -r {{ rst }} -x {{ nc }}
