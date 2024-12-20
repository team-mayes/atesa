Getting Started with ATESA
==========================

.. _Installation:

Installation
------------

ATESA should be installed directly on the high-performance computing (HPC) resource that you intend to use. The simplest way to do so is with pip::

	pip install atesa
	
You may need to append the ``--user`` option depending on how your HPC resource is configured. ATESA requires Python 3, and the above pip command may fail if your active python environment defaults to Python 2.

Alternatively, you can download or clone ATESA directly from `GitHub <https://github.com/team-mayes/atesa>`_ and install using the ``setup.py`` script in its root directory directly::

	python setup.py install --user
	
If you're using a custom python environment, remember to activate it before installing!

Once you've installed ATESA, you may want to head over to the :ref:`ExampleStudy` page to take a look at what ATESA can do and get a feel for its basic operation.

If you want to perform umbrella sampling simulations with ATESA, you need to make sure that PLUMED is available. Versions of Amber released since 2015 have native PLUMED support, so once it has been installed (see `this link <https://www.plumed.org/doc-v2.5/user-doc/html/_installation.html>`_ for instructions), all you need to do is set the PLUMED_KERNEL environment variable to the appropriate path for ``libplumedKernel.so``, like so::

	export PLUMED_KERNEL=/path/to/libplumedKernel.so

Usage
-----

ATESA is designed to dynamically handle jobs on a PBS/Torque or Slurm batch system on a high-performance computing cluster. It is invoked as:

.. code-block:: bash

   atesa <config_file> [<working_directory>]
   
The config_file parameter is required and supports absolute or relative paths. The working_directory parameter is optional, but if it is provided it overrides the value in the configuration file. This option is made available for use on systems that only allocate working space for jobs after they have been initialized. For detailed documentation on the contents of the configuration file, see :ref:`TheConfigFile` (but it's easier to get started by modifying the examples in the :ref:`ExampleStudy`).

Although ATESA can be run directly from the command line, because its process continues for the entirety of the job it is usually best to submit it as its own batch job, or to run it on an interactive resource allocation. A single core and a modest allocation of memory should be sufficient to run ATESA on most platforms (it is typically neither memory- nor processor-intensive, although many features involve a significant amount of I/O). ATESA supports multiprocessing during aimless shooting on UNIX-based systems; this may improve performance in some cases, depending on ATESA settings, simulation speed, and the batch queue. Simply allocate the desired number of cores and ATESA will use them as efficiently as it can.

Note that almost all ATESA jobs will produce a large amount of data, at least transiently. For this reason, the working directory should probably be set to a path inside the "scratch" filesystem or equivalent on your HPC cluster. If you aren't sure what this means, consult the documentation or support staff for your specific resource.

.. _SettingUpSimulationFiles:

Setting Up Simulation Files
---------------------------

Beyond the installation of ATESA itself and the creation of the configuration file, the only other step required before performing simulations is to provide the program with the necessary input files and batch script templates for those simulations. ATESA looks for two directories:

* 'input_files', which contains input files for various types of simulations; and
* 'templates', which contains batch script files formatted for use as templates

ATESA comes packaged with examples of both of these directories, located inside the 'atesa/data' directory in the ATESA installation folder. The contents of these directories will be used by ATESA by default, although the user can specify other directories instead (see :ref:`FilePathSettings`). Note that the default scripts are certainly **not** appropriate for any given user's needs; they are provided as examples only. It is suggested that the user take a look at these files and the following documentation on this page before beginning to construct their own versions, in order to better understand what each file is for.

Input Files
~~~~~~~~~~~

The 'input_files' directory contains input files for molecular simulations. Each file should be named according to the scheme::

	<job_type>_<step_type>_<md_engine>.in
	
Each of the bracketed values should be replaced by the appropriate name. The "job_type" is identical to the chosen value of the option of the same name in the configuration file (see :ref:`CoreSettings`). The step type is either "init" or "prod": "prod" is used in every job type and is the primary simulation step, while "init" is used only in aimless shooting and equilibrium path sampling. Finally, "md_engine" should be the name of the simulations engine to be used (same as the option of the same name in the configuration file). For example, the file for the "prod" step of an "aimless_shooting" job performed with Amber would be::

	aimless_shooting_prod_amber.in

Every input file must be appropriate for your specific molecular model, and also for the specific job and step types. Following are details for each job type describing what each input file is used for in that job, and any necessary characteristics of that input file. This may seem daunting at first, but once you are comfortable with the basic usage of the MD engine you wish to use, it should be quite straightforward to construct input files appropriate for your model!

``aimless_shooting``

Aimless shooting input files for the following step types are required for jobs with job_type "aimless_shooting" or "find_ts":

* **aimless_shooting_init_amber.in**: Aimless shooting "init" steps are extremely short simulations whose only purpose is to obtain a fresh set of initial velocities for the following "prod" step. To this end, the "init" input file should be configured to generate new initial velocities from the Boltzmann distribution (as opposed to using velocities from the input coordinate file), and to only perform a single simulation step with an extremely small time step. In Amber, the following settings should be specified in the &cntrl namelist, in addition to any other model-specific settings ("!" denotes a comment in Amber input files)::
		
	ntx=1,		! read coordinates but not velocities from input coordinate file
  	ntxo=1,		! ASCII-formatted restart file (required for ATESA)
  	nstlim=1,	! one simulation step total
	dt=0.00001,	! extremely short time step (too small can cause velocity overflow errors)
  	tempi=300.0,	! or whatever temperature (same as temp0)
  	temp0=300.0,	! or whatever temperature (same as tempi)
	ig=-1		! random initial velocities; should be default but may not be in some versions/patches of Amber

* **aimless_shooting_init_cp2k.in**: Unlike in Amber, CP2K input files specify the names of input and output files directly (this is done in the batch script/on the command line with Amber). As such template slots for these are required in the CP2K input files. The box dimensions are also templated. Unless otherwise described, the following lines are required (note that depending on use-case, you may need additional settings within these sections)::

    &GLOBAL
      PROJECT AS_INIT
      PRINT_LEVEL LOW
      RUN_TYPE MD
      SEED {{ seed }}
    &END GLOBAL
    [...]
    &CELL   ! within the &FORCE_EVAL > &SUBSYS section
      ABC [angstrom] {{ box_xyz }}
      ALPHA_BETA_GAMMA {{ box_abc }}
    &END CELL
    &TOPOLOGY   ! within the &FORCE_EVAL > &SUBSYS section
      COORD_FILE_FORMAT pdb     ! ATESA will convert input files to .pdb format
      COORD_FILE_NAME {{ inpcrd }}
      CONN_FILE_FORMAT amber    ! ATESA may support other topology formats for use with CP2K, but only Amber has been tested
      CONN_FILE_NAME {{ prmtop }}
    &END TOPOLOGY
    [...]
    &MD     ! within the &MOTION > &MD section
      STEPS 1   ! one simulation step for initialization jobs
      TIMESTEP 0.001    ! extremely short timestep
    &END MD
    [...]
    &TRAJECTORY   ! within the &MOTION > &PRINT section
      FILENAME {{ nc }}
      FORMAT DCD    ! ATESA expects output files from CP2K to be in DCD format
    &END TRAJECTORY
    &RESTART   ! within the &MOTION > &PRINT section
      FILENAME {{ rst }}
      &EACH
        MD 1    ! so that a restart file is produced after the first (and only) simulation step
      &END
    &END RESTART

* **aimless_shooting_prod_amber.in**: Aimless shooting "prod" steps are the primary simulation steps for each shooting move. They take the initial coordinates and velocities from an "init" step (with velocities reversed in the case of backward trajectories) and run until the simulation commits to one of the stable states defined in the configuration file. Therefore, the time step and number of simulation steps should be much larger than in an "init" simulation. In Amber, the following settings should be specified in the &cntrl namelist, instead of the above "init" settings and in addition to any other model-specific settings::

	ntx=5,		! read coordinates AND velocities from input coordinate file
  	ntxo=1,		! ASCII-formatted restart file (required for ATESA)
  	nstlim=5000,	! a large maximum number of steps; will probably be terminated early
  	dt=0.001,	! or whatever desired simulation time step
  	irest=1,	! restart simulation from preceding "init" step
  	temp0=300.0,	! or whatever temperature
  	ntwx=1,		! or whatever trajectory write frequency, but not only at the end of the simulation
  	ntwv=-1,	! include velocities in trajectory files (required if the option "include_qdot" is True (which is default)

* **aimless_shooting_prod_cp2k.in**: In addition to template slots for box dimensions and input and output files, aimless shooting "prod" steps require initial velocities, which in CP2K are given explicitly by-atom in the input file. ATESA handles this with a velocities template slot::

    &GLOBAL
      PROJECT AS_PROD
      PRINT_LEVEL LOW
      RUN_TYPE MD
      SEED {{ seed }}
    &END GLOBAL
    [...]
    &CELL   ! within the &FORCE_EVAL > &SUBSYS section
      ABC [angstrom] {{ box_xyz }}
      ALPHA_BETA_GAMMA {{ box_abc }}
    &END CELL
    &TOPOLOGY   ! within the &FORCE_EVAL > &SUBSYS section
      COORD_FILE_FORMAT pdb     ! ATESA will convert input files to .pdb format
      COORD_FILE_NAME {{ inpcrd }}
      CONN_FILE_FORMAT amber    ! ATESA may support other topology formats for use with CP2K, but only Amber has been tested
      CONN_FILE_NAME {{ prmtop }}
    &END TOPOLOGY
    &VELOCITY   ! within the &FORCE_EVAL > &SUBSYS section
      {{ velocities }}
    &END VELOCITY
    [...]
    &MD     ! within the &MOTION > &MD section
      STEPS 1000   ! one simulation step for initialization jobs
      TIMESTEP 1    ! or whatever desired simulation timestep
    &END MD
    [...]
    &TRAJECTORY   ! within the &MOTION > &PRINT section
      FILENAME {{ nc }}
      FORMAT DCD    ! ATESA expects output files from CP2K to be in DCD format
    &END TRAJECTORY

``committor_analysis``

Only a "prod" committor analysis input file is required for jobs with job_type "committor_analysis":

* **committor_analysis_prod_amber.in**: These jobs can use exactly the same settings as aimless shooting "prod" steps, except that each simulation should obtain new initial velocities, as in an aimless shooting "init" steps. In Amber, that means that these three options should be set as follows::

	ntx=1,		! read coordinates but not velocities from input coordinate file
	tempi=300.0,	! or whatever temperature (same as temp0)
	irest=0,	! do not restart, use new velocities (this is the default)

* **committor_analysis_prod_cp2k.in**: Identical to an aimless shooting CP2K input file, but without the *velocities* template (new velocities should be drawn each time).

``umbrella_sampling``

Only a "prod" umbrella sampling input file is required for jobs with job_type "umbrella_sampling":

* 	**umbrella_sampling_prod_amber.in**: By default, umbrella sampling restraints are applied using PLUMED. The umbrella sampling input file can be almost identical to a committor analysis "prod" file, with the following changes::

		nstlim=10000,	! a "large" maximum number of steps (enough for convergent sampling, may take some trial and error)
		plumed=1,		! enable plumed backend
		plumedfile={{ plumedfile }},		! template slot for declaring plumed file

*   **umbrella_sampling_prod_cp2k.in**: CP2K umbrella sampling supports PLUMED only. The file can be identical to a committor analysis "prod" input file, with the following changes::

        &GLOBAL
          PROJECT US_PROD
          PRINT_LEVEL LOW
          RUN_TYPE MD
          SEED {{ seed }}
        &END GLOBAL
        [...]
        &METADYN    ! within the &MOTION > &FREE_ENERGY section
          USE_PLUMED .TRUE.
          PLUMED_INPUT_FILE {{ plumedfile }}
        &END METADYN

In either case, ATESA will write the appropriate plumed file automatically and insert a reference to it into the input file as needed.
	
``equilibrium_path_sampling``

Equilibrium path sampling input files for the following step types are required for jobs with job_type "equilibrium_path_sampling":

* **equilibrium_path_sampling_init_amber.in**: Equilibrium path sampling "init" steps are functionally identical to aimless shooting "init" steps and can use an identical input file.

* **equilibrium_path_sampling_init_cp2k.in**: Equilibrium path sampling "init" steps are functionally identical to aimless shooting "init" steps and can use an identical input file.

* **equilibrium_path_sampling_prod_amber.in**: This input file can be functionally identical to an aimless shooting "prod" input file, with two key exceptions: the number of simulation steps must be replaced with the exact string ``{{ nstlim }}`` and the frequency of writes to the output trajectory must be replaced with the exact string ``{{ ntwx }}``. In Amber::
	
    nstlim={{ nstlim }},
    ntwx={{ ntwx }},

* **equilibrium_path_sampling_prod_cp2k.in**: This input file can be functionally identical to an aimless shooting "prod" input file, with two key exceptions: the number of simulation steps must be replaced with the exact string ``{{ nstlim }}`` and the frequency of writes to the output trajectory must be replaced with the exact string ``{{ ntwx }}``. In CP2K::

    STEPS {{ nstlim }}      ! within the &MOTION > &MD section
    [...]
    &EACH       ! within the &MOTION > &PRINT section
      MD {{ ntwx }}
    &END

		
``find_ts``

In addition to making sure that appropriate input files for aimless shooting are available, jobs with job_type "find_ts" require their own "prod" input file:

* **find_ts_prod_amber.in** This file can be mostly identical to the "aimless_shooting" prod input file, with two key additions: there must be a restraint specified using the file "find_ts_restraints.disang", and the weight of the restraint must be set to steadily increase over time (beginning from zero). An example of a working implementation of this in Amber is as follows. Options in the &cntrl namelist that can be the same as in aimless shooting are here replaced by an elipse ([...]) for brevity, but they must still be explicitly specified in the input file. Other than that, it should be quite safe to copy the rest of this exactly into your Amber "find_ts" input file, or customize it as you see fit::

	 &cntrl
	  [...]
	  nmropt=1,		! turn on restraints
	 &end
	 &wt
  	  type="REST",
  	  istep1=0,
  	  istep2=1000,
  	  value1=0,
  	  value2=1,
 	 &end
 	 &wt
  	  type="REST",
  	  istep1=1001,
  	  istep2=2000,
  	  value1=1,
  	  value2=1,
 	 &end
 	 &wt
  	  type="END",
 	 &end
	DISANG=find_ts_restraints.disang

* **find_ts_prod_cp2k.in** This file can be mostly identical to the "aimless_shooting" prod input file, with three changes: the number of simulation steps should be templated as "nstlim", and template slots are required to define collective variables ("colvars") and constraints ("collective") as follows::

    &GLOBAL
      PROJECT FIND_TS_PROD
      PRINT_LEVEL LOW
      RUN_TYPE MD
      SEED {{ seed }}
    &END GLOBAL
    [...]
    {{ colvars|safe }}   ! within the &FORCE_EVAL > &SUBSYS section
    [...]
    &CONSTRAINT     ! within the &MOTION section
      CONSTRAINT_INIT TRUE
      {{ collective|safe }}
    &END CONSTRAINT
    [...]
    STEPS {{ nstlim }}  ! within the &MOTION > &MD section

Templates
~~~~~~~~~

The 'templates' directory contains files that ATESA will automatically customize for each individual simulation. It is used for templated batch scripts that will be filled using the :ref:`BatchTemplateSettings` in the configuration file, in addition to several other keywords, described below.

Batch script templates should be named according to the scheme::

	<md_engine>_<batch_system>.tpl
	
Each of the bracketed values should be replaced by the appropriate name. The "md_engine" and "batch_system" should each be the same as the setting picked for the option of the same name in the configuration file. For example, the Slurm batch template for a simulation with Amber would be named::

	amber_slurm.tpl
	
In general, the overwrite flag ("-O") should always be present in the command invoking Amber with ATESA, or certain features may not work.

Template slots are delimited by double curly braces, as in "{{ example }}". Anything not delimited in this way will be taken as literal. The user should provide batch files that will work for their particular software environment, making use of the templates wherever the batch file might differ between simulations -- please refer to the 'atesa/data/templates' directory for examples. In addition to the relevant configuration file settings (again, see :ref:`BatchTemplateSettings`), the following keywords should be included in batch script templates for ATESA to fill out automatically. It is safe to omit any of these keywords if you are sure that a fixed value (or no value at all) is appropriate instead. These values are also used to fill the template slots in input files for CP2K.

``{{ name }}``

The name of the batch job. This will be unique to each step of each thread.

``{{ inp }}``

The input file for this simulation (*e.g.*, one of the files from the 'input_files' directory).

``{{ out }}``

The output/log file for this simulation.

``{{ prmtop }}``

The parameter/topology file for this simulation (the file indicated for the "topology" option in the configuration file).

``{{ inpcrd }}``

The initial coordinate file for this simulation.

``{{ rst }}``

The output coordinate or restart file from this simulation.

``{{ nc }}``

The output trajectory file from this simulation.
