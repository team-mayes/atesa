Getting Started with ATESA
==========================

Installation
------------

[Insert finalized installation instructions here]

Usage
-----

ATESA is designed to dynamically handle jobs on a PBS/Torque or Slurm batch system. It is invoked as:

.. code-block:: bash

   atesa <config_file> [<working_directory>]
   
The config_file parameter is required and supports absolute or relative paths. The working_directory parameter is optional, but if it is provided it overrides the value in the configuration file. This option is made available for use on systems that only allocate working space for jobs after they have been initialized. For details on the contents of the configuration file, see :ref:`TheConfigFile`.

Although ATESA can be run directly from the command line, because its process continues for the entirety of the job it is usually recommended to submit it as its own batch job. A single core and a modest allocation of memory should be sufficient to run ATESA on most platforms (it is neither memory- nor processor-intensive, although certain features involve a significant amount of I/O, and more memory may be necessary when using the information error termination criterion).

.. _SettingUpSimulationFiles:

Setting Up Simulation Files
---------------------------

Beyond the installation of ATESA itself and the creation of the configuration file, the only other step required before performing simulations is to provide the program with the necessary input files and batch scripts for those simulations. ATESA looks for two directories:

* 'input_files', which contains inputs for the molecular simulations software of choice (*e.g.*, Amber); and
* 'templates', which contains files formatted for use as templates (see below for details)

ATESA comes packaged with examples of both of these directories, located inside the 'data' directory in the ATESA installation folder. The contents of these directories will be used by ATESA by default, although the user can specify other directories instead (see :ref:`FilePathSettings`). Note that the default scripts are certainly **not** appropriate for any given user's needs; they are provided as examples only. It is suggested that the user take a look at these files before beginning to construct their own versions, in order to better understand what each file is for.

Input Files
~~~~~~~~~~~

The 'input_files' directory contains input files for molecular simulations. Each file should be named according to the scheme::

	<job_type>_<step_type>_<md_engine>.in
	
Each of the bracketed values should be replaced by the appropriate name. The "job_type" is identical to the chosen value of the option of the same name in the configuration file (see :ref:`CoreSettings`). The step type is either "init" or "prod": "prod" is used in every job type and is the primary simulation step, while "init" is used only in aimless shooting and equilibrium path sampling. Finally, "md_engine" should be the name of the simulations engine to be used (at this time, only "amber" is supported). For example, the file for the "prod" step of an "aimless_shooting" job performed with Amber would be::

	aimless_shooting_prod_amber.in

Every input file must be appropriate for your specific molecular model, and also for the specific job and step types. As Amber is the only simulations software currently supported, users who are not familiar with Amber are directed to the `relevant tutorials <https://ambermd.org/tutorials/>`_. Following are details for each job type describing what each input file is used for in that job, and any necessary characteristics of that input file:

``aimless_shooting``

Aimless shooting input files for the following step types are required for jobs with job_type "aimless_shooting" or "find_ts".

* **init**: Aimless shooting "init" steps are extremely short simulations whose only purpose is to obtain a fresh set of initial velocities for the following "prod" step. To this end, the "init" input file should be configured to generate new initial velocities from the Boltzmann distribution (as opposed to using velocities from the input coordinate file), and to only perform a single simulation step with an extremely small time step. In Amber, the following settings should be specified in the &cntrl namelist, in addition to any other model-specific settings ("!" denotes a comment in Amber)::
		
	ntx=1,		! read coordinates but not velocities from input coordinate file
  	ntxo=1,		! ASCII-formatted restart file
  	nstlim=1,	! one simulation step total
	dt=0.0000001,	! extremely short time step
  	tempi=300.0,	! or whatever temperature (same as temp0)
  	temp0=300.0,	! or whatever temperature (same as tempi)
  		
* **prod**: Aimless shooting "prod" steps are the primary simulation steps for each shooting move. They take the initial coordinates and velocities from an "init" step (with velocities reversed in the case of backward trajectories) and run until the model commits to one of the stable states defined in the configuration file. Therefore, the time step and number of simulation steps should be much larger than in an "init" simulation. In Amber, the following settings should be specified in the &cntrl namelist, instead of the above "init" settings and in addition to any other model-specific settings::

	ntx=5,		! read coordinates AND velocities from input coordinate file
  	ntxo=1,		! ASCII-formatted restart file
  	nstlim=5000,	! a large maximum number of steps; will be terminated early
  	dt=0.001,	! or whatever desired simulation time step
  	irest=1,	! restart simulation from preceding "init" step
  	temp0=300.0,	! or whatever temperature
  	
``committor_analysis``

A committor analysis input file for the following step type is required for jobs with job_type "committor_analysis" only.

* **prod**: Committor analysis only consists of "prod" steps. These jobs can use exactly the same settings as aimless shooting "prod" steps, except that each simulation should obtain new velocities, as in an aimless shooting "init" steps. In Amber, that means that these three options should be set as follows::

	ntx=1,		! read coordinates but not velocities from input coordinate file
	tempi=300.0,	! or whatever temperature (same as temp0)
	irest=0,	! do not restart, use new velocities (this is the default)
	
``equilibrium_path_sampling``

Equilibrium path sampling input files for the following step types are required for jobs with job_type "equilibrium_path_sampling" only.

* **init**: Equilibrium path sampling "init" steps are functionally identical to aimless shooting "init" steps and can use an identical input file.

*
	**prod**: Equilibrium path sampling "prod" steps are the only type of job currently supported by ATESA that does *not* take its input file from the "input_files" directory. Instead, the input file is constructed from the file in the "templates" directory named as:
	
		::
			
			<md_engine>_eps_in.tpl
	
	This input file can be functionally identical to an aimless shooting "prod" input file, with two key exceptions: the number of simulation steps must be replaced with the exact string ``{{ nstlim }}`` and the frequency of writes to the output trajectory must be replaced with the exact string ``{{ ntwx }}``. In Amber::
	
		nstlim={{ nstlim }},
		ntwx={{ ntwx }},
		
``find_ts``

A find_ts input file for the following step type is required for jobs with job_type "find_ts" only.

* **prod** "find_ts" jobs consist only of "prod" steps. This file can be mostly identical to the "aimless_shooting" prod input file, with two key additions: there must be a restraint specified using the file "find_ts_restraints.disang", and the weight of the restraint must be set to steadily increase over time (beginning from zero). An example of a working implementation of this in Amber is as follows (options in the &cntrl namelist that can be the same as in aimless shooting are here replaced by an elipse (...) for brevity, but they must still be explicitly specified in the input file)::

	 &cntrl
	  ...
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
	
Templates
~~~~~~~~~

The 'templates' directory contains files that ATESA will automatically customize for each individual simulation. It is primarily used for templated batch scripts that will be filled using the :ref:`BatchTemplateSettings` in the configuration file, in addition to several internal keywords.

Batch script templates should be named according to the scheme::

	<md_engine>_<batch_system>.tpl
	
Each of the bracketed values should be replaced by the appropriate name. The "md_engine" should be the name of the simulations engine to be used (at this time, only "amber" is supported). The "batch_system" should be the same as the setting picked for the option of the same name in the configuration file (either "slurm" or "pbs" are currently supported). For example, the Slurm batch template for a simulation with Amber would be::

	amber_slurm.tpl
	
Template slots are delimited by double curly braces, as in "{{ example }}". Anything not delimited in this way will be taken as literal. The user should provide batch files that will work for their particular software environment, making use of the templates wherever the call to the molecular simulation software would differ between simulations. In addition to the configuration file settings (again, see :ref:`BatchTemplateSettings`), the following keywords should be included in batch script templates for ATESA to fill out automatically. It is safe to omit any of these keywords if you are sure that a fixed value is appropriate instead.

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

The output (final) coordinate file from this simulation.

``{{ nc }}``

The output trajectory file from this simulation.

As indicated in the preceding section, the 'templates' directory should also include a template file for the equilibrium path sampling "prod" step input file, if equilibrium path sampling is to performed.