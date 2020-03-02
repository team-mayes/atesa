.. _TheConfigFile:

The Configuration File
======================

The configuration file is the primary means of controlling the behavior of ATESA. In order to support the wide array of functionality that any given user may need, the configuration file supports many options and can be quite long; however, in most cases a relatively short configuration file will be sufficient. This page provides some recommendations for building the configuration file for a handful of common use cases, and provides detailed documentation for each setting.

The contents of the configuration file are read into ATESA as literal python code, which enables invocation of python built-in functions as well as methods of pytraj and numpy. **Warning**: This input is not sanitized in any way. For this reason among others, "shutil.rmtree('/')" makes for a poor reaction coordinate!

.. _CoreSettings:

Core Settings
-------------

Certain settings should be given regardless of the job type. The following settings do not have valid default values and should be set in every configuration file:

.. code-block:: python

	batch_system		# Valid options: 'slurm', 'pbs'
	restart			# Valid options: True, False
	overwrite		# Valid options: True, False
	topology		# Valid options: Absolute or relative path as a string
	working_directory	# Valid options: Absolute or relative path as a string
	
**batch_system**
	Indicates the type of batch system on which ATESA is running. Supported options are 'slurm' and 'pbs' (the latter is also known as TORQUE).
	
**restart**
	Indicates whether this is a new job (False) or a continuation of an old one *in the same working directory* (True).
	
**overwrite**
	Indicates whether to delete the existing working directory (if one exists) and create a new one *if and only if* restart = False.
	
**topology**
	An absolute or relative path given as a string and pointing to the simulation topology file.
	
**working_directory**
	An absolute or relative path given as a string and pointing to the desired working directory (this can be omitted if the working directory is set in the command line). This is the directory in which all of the simulations will be performed. It will be created if it does not exist.

Common Use Cases
----------------

Aimless Shooting
~~~~~~~~~~~~~~~~

The most central function of ATESA is aimless shooting. If you have one or more initial transition structure guesses and want to begin transition path sampling, you should call ATESA with the following settings (in addition to the :ref:`CoreSettings` above):

.. code-block:: python

	job_type = 'aimless_shooting'
	initial_coordinates = [<coord_file_1>, <coord_file_2>, ...]
	cvs = [<cv1>, <cv2>, ...]
	commit_fwd = [...]
	commit_bwd = [...]
	
See :ref:`CommitmentBasinDefinitions` for details on the `commit_fwd` and `commit_bwd` options and :ref:`CVDefinitions` for details on the `cvs` option.

These settings will automatically use :ref:`InformationError` as a termination criterion with the default settings. These can be modified as described in :ref:`AimlessShootingSettings`

Find Transition State
~~~~~~~~~~~~~~~~~~~~~

If you want to perform aimless shooting but don't have a candidate transition state structure from which to begin (or just want to obtain more), ATESA can automatically build
them from a given product or reactant state structure. In addition to the :ref:`CoreSettings`, such a job should be set up as:

.. code-block:: python

	job_type = 'find_ts'
	initial_coordinates = [<coord_file_1>]
	commit_fwd = [...]
	commit_bwd = [...]
	max_moves = 5
	
See :ref:`CommitmentBasinDefinitions` for details on the `commit_fwd` and `commit_bwd` options.

Committor Analysis
~~~~~~~~~~~~~~~~~~

After completing aimless shooting, the next step is to obtain a reaction coordinate and verify it with committor analysis. Before running committor analysis, the user should call  likelihood maximization with their preferred settings (see :ref:`LikelihoodMaximization`) and then reaction coordinate evaluation using the aimless shooting working directory and the resulting reaction coordinate (see :ref:`RCEval`). 

Then, committor analysis is performed through the main ATESA executable with the following settings:

.. code-block:: python

	job_type = 'committor_analysis'
	path_to_rc_out = <rc_out_file>
	rc_definition = <rc_definition>
	cvs = [<cv1>, <cv2>, ...]
	commit_fwd = [...]
	commit_bwd = [...]
	
See :ref:`CommitmentBasinDefinitions` for details on the `commit_fwd` and `commit_bwd` options, :ref:`CVDefinitions` for details on the `cvs` option, and :ref:`ReactionCoordinateDefinition` for details on the `rc_definition` option.

The working directory here should NOT be the same as the aimless shooting directory containing the data to perform committor analysis with. The aimless shooting directory will be identified using the `path_to_rc_out` setting (which should be inside the aimless shooting working directory). The working directory for a committor analysis job should be a new directory, though it can be a subdirectory of the aimless shooting directory it operates on.

Equilibrium Path Sampling
~~~~~~~~~~~~~~~~~~~~~~~~~

The final analysis step after a satisfactory committor analysis run is to obtain the free energy profile along the reaction coordinate. ATESA supports equilibrium path sampling (EPS) to obtain this profile through the main executable, using the following settings:

.. code-block:: python

	job_type = 'equilibrium_path_sampling'
	rc_definition = <rc_definition>
	cvs = [<cv1>, <cv2>, ...]
	
See  :ref:`CVDefinitions` for details on the `cvs` option and :ref:`ReactionCoordinateDefinition` for details on the `rc_definition` option.

EPS is a highly generalized free energy method that does not rely on restraints or biases of any kind. The cost of this benefit is that it is also among the least efficient free energy methods available, requiring a relatively large amount of simulation to acquire comparable sampling coverage to, for example, umbrella sampling. For this reason, EPS is recommended for use only in cases where other methods are unsuitable, such as for example in cases of highly complex reaction coordinates that do not lend themselves to restraints.

CAUTION: Because equilibrium path sampling measures the full energy profile instead of merely assessing the endpoints of simulations (as in aimless shooting and committor analysis), it is very sensitive to errors in the evaluation of the energy of any given state. For this reason, it is completely possible to have obtained reasonable aimless shooting and committor analysis results with a system or simulation parameters that are not suitable for equilibrium path sampling, for example owing to poor SCF convergence in QM calculations along portions of the RC. ATESA can NOT identify such errors on its own, and may produce EPS results that are not correct (but may appear reasonable at first glance)! It is the responsibility of the user to ensure that the EPS simulations are well-behaved and do not suffer from severe energetic errors.

The raw output data from an EPS run is stored in the working directory as 'eps.out'. This data can be converted into an energy profile using boltzmann_weight.py (see :ref:`BoltzmannWeight`), which calculates the relative probabilities of states within each bin and converts these into relative free energies.
	

Full Configuration Options
--------------------------

Here, the full list of valid configuration file options are documented (excluding the :ref:`CoreSettings`, documented above.)

Batch Template Settings
~~~~~~~~~~~~~~~~~~~~~~~

These settings control how batch file templates are filled out. In general, 'init' simulations are very short (one step) and only used to prepare something for a longer simulation, whereas 'prod' simulations are where the primary data collection of a job takes place. When and whether 'init' or 'prod' variables are used is controlled by the job type.

These settings are used to fill in the template slots in the user-provided template files. If you do not wish for ATESA to use an option, you can simply omit its template slot from the appropriate file and leave it unset in the configuration file.

**init_nodes**
	The number of compute nodes to request for 'init' simulations, given as an integer. Default = 1

**init_ppn**
	The number of cores or processes to request per node (ppn: "processes per node") for 'init' simulations, given as an integer. Default = 1
	
**init_mem**
	The amount of RAM to request for 'init' simulations, given as a string of appropriate format for the batch system. Depending on the batch system, this may be interpreted as total memory, or as memory per core. Default = '4000mb'
	
**init_walltime**
	The amount of walltime (real time limit for the batch job) to request for 'init' simulations, given as a string of appropriate format for the batch system. 'init' simulations are very quick, but if one does not produce the necessary results in this time and is cancelled early, it will be resubmitted; err on the side of more time. Default = '00:30:00'
	
**init_solver**
	The name of the executable to use to perform 'init' simulations, given as a string. Default = 'sander' (which is specific to Amber)
	
**init_extra**
	An additional template slot for 'init' simulations to be used however the user sees fit. This option is provided in case a user has an unforseen need to template something other than the above options. Default = '' (an empty string)
	
**prod_nodes**
	The number of compute nodes to request for 'prod' simulations, given as an integer. Default = 1

**prod_ppn**
	The number of cores or processes to request per node (ppn: "processes per node") for 'prod' simulations, given as an integer. Default = 8
	
**prod_mem**
	The amount of RAM to request for 'prod' simulations, given as a string of appropriate format for the batch system. Depending on the batch system, this may be interpreted as total memory, or as memory per core. Default = '4000mb'
	
**prod_walltime**
	The amount of walltime (real time limit for the batch job) to request for 'prod' simulations, given as a string of appropriate format for the batch system. Err on the side of more time. Default = '02:00:00'
	
**prod_solver**
	The name of the executable to use to perform 'prod' simulations, given as a string. Default = 'sander' (which is specific to Amber)
	
**prod_extra**
	An additional template slot for 'prod' simulations to be used however the user sees fit. This option is provided in case a user has an unforseen need to template something other than the above options. Default = '' (an empty string)

Path Settings
~~~~~~~~~~~~~

These settings define the paths where ATESA will search for user-defined input files and template files.

**path_to_input_files**
	Path (as a string) to the directory containing the input files. Default = sys.path[0] + '/data/input_files'

**path_to_templates**
	Path (as a string) to the directory containing the template files. Default = sys.path[0] + '/data/templates'

.. _CVDefinitions:

CV Settings
~~~~~~~~~~~

These settings define the combined variables (CVs) for the job. In aimless shooting, these are the values that are written to the output file for interpretation by likelihood maximization in building the reaction coordinate (RC). In equilibrium path sampling and committor analysis, they are used to evaluate the RC (see :ref:`ReactionCoordinateDefinition`).

**cvs**
	A list of CV definitions, given as strings, as in: [<cv1>, <cv2>, ... <cvN>] (where the contents of each pair of angled braces is a string). Each item is interpreted as raw python code (caution: unsanitized) that returns the desired CV value (in a format that can cast to a float). In addition to built-in python functions, calls to pytraj, mdtraj, and numpy (as both 'numpy' and 'np') are supported. In support of pytraj and mdtraj calls, the following variables are available for use:
	
		*traj*: the coordinate file being evaluated, as a pytraj.iterload object
	
		*traj_name*: the name of the coordinate file as a string
	
		*settings.topology*: the name of the topology file, as a string
	
	These CV evaluations are only ever performed on coordinate files with a single frame (not multi-frame trajectories). It is up to the user to ensure that each item in *cvs* returns exactly the desired CV value (and not, for example, a one-length list containing that value). For example, the following value of *cvs* would interpret the interatomic distance between atoms 1 and 2 as CV1, the difference between the distances 3-to-4 and 5-to-6 as CV2, and the angle formed by atoms 7-8-9 as CV3 (all on one line):
	
        ::

	    cvs = ['pytraj.distance(traj, \'@1 @2\')[0]',
	    'pytraj.distance(traj, \'@3 @4\')[0] - pytraj.distance(traj, \'@5 @6\')[0]',
	    'pytraj.angle(traj, \'@7 @8 @9\')[0]']

	Notice in particular the usage of "traj" as the first argument in these pytraj function calls, the escaped single-quote characters within each function call, and the call to the zero'th index of each function call (as these pytraj functions return one-length lists). Default = ['']

**include_qdot**
	A boolean indicating whether the instantaneous rate of change of each CV defined in *cvs* should also be counted as a CV. These values can be used in inertial likelihood maximization to favor reaction coordinates with high transmission coefficients (see `Peters 2012 <https://doi.org/10.1016/j.cplett.2012.10.051>`_). Rate of change CVs come after the rest of the CVs in the aimless shooting output file, as in:
	
	::
	
	<basin> <- <CV1> <CV2> ... <CVn> <qdot1> <qdot2> ... <qdotn>
	
	Default = True


Initial Coordinates
~~~~~~~~~~~~~~~~~~~
        initial_coordinates: typing.List[str] = ['']

.. _CommitmentBasinDefinitions:

Commitment Basin Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        commit_fwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])
        commit_bwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])
        
    	commit_fwd = ([<mask_f_1_1>, <mask_f_1_2>, ...], [<mask_f_2_1>, <mask_f_2_2>, ...], [<dist_f_1>, <dist_f_2>, ...], [<'gt'/'lt'>, <'gt'/'lt'>, ...])
		commit_bwd = ([<mask_b_1_1>, <mask_b_1_2>, ...], [<mask_b_2_1>, <mask_b_2_2>, ...], [<dist_b_1>, <dist_b_2>, ...], [<'gt'/'lt'>, <'gt'/'lt'>, ...])

.. _ReactionCoordinateDefinition:

Reaction Coordinate Definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rc_definition: str = ''
        as_out_file: str = 'as_raw.out'
        rc_reduced_cvs: bool = True

.. _AimlessShootingSettings:

Aimless Shooting Settings
~~~~~~~~~~~~~~~~~~~~~~~~~
        min_dt: int = 1
        max_dt: int = 10
        always_new: bool = True
        resample: bool = False
        degeneracy: int = 1
        cleanup: bool = True
        information_error_checking: bool = True
        information_error_freq: int = 250
        information_error_override: bool = False
        information_error_max_dims: int = 6
        max_moves: int = -1     # also used by find_ts
        max_consecutive_fails: int = 4

Committor Analysis Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~
        committor_analysis_n: int = 10
        committor_analysis_use_rc_out: bool = True
        path_to_rc_out: str = sys.path[0] + '/atesa_v2/tests/test_data/rc.out'
        rc_threshold: float = 0.05

Equilibrium Path Sampling Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eps_rc_min: float = -12
        eps_rc_max: float = 12
        eps_rc_step: float = 1
        eps_rc_overlap: float = 0.1
        eps_n_steps: int = 6
        eps_out_freq: int = 1
        eps_dynamic_seed: typing.Union[int, list] = 20  # int or list (int -> [int for window in eps_windows]; 0 or empty list turns off)
        samples_per_window: int = -1

Other Settings
~~~~~~~~~~~~~~
        restart_terminated_threads: bool = False
