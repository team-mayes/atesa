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
	Indicates the type of batch system on which ATESA is running. Supported options are 'slurm' and 'pbs' (also known as TORQUE).
	
**restart**
	Indicates whether this is a new job (False) or a continuation of an old one *in the same working directory* (True).
	
**overwrite**
	Indicates whether to delete the existing working directory (if one exists) and create a new one *if and only if* restart = False.
	
**topology**
	Absolute or relative path given as a string and pointing to the simulation topology file.
	
**working_directory**
	Absolute or relative path given as a string and pointing to the desired working directory (this can be omitted if the working directory is set in the command line). This is the directory in which all of the simulations will be performed. It will be created if it does not exist.

Common Use Cases
----------------

Aimless Shooting
~~~~~~~~~~~~~~~~

The most central function of ATESA is aimless shooting. If you have one or more initial transition structure guesses and want to begin transition path sampling, you should call ATESA with the following settings (in addition to the :ref:`CoreSettings` above):

.. code-block:: python

	job_type = aimless_shooting
	initial_coordinates = [<coord_file_1>, <coord_file_2>, ...]
	
These settings will automatically use :ref:`InformationError` as a termination criterion with the default settings. These can be modified as described in :ref:`AimlessShootingSettings`

Find Transition State
~~~~~~~~~~~~~~~~~~~~~

If you want to perform aimless shooting but don't have a candidate transition state structure from which to begin (or just want to obtain more), ATESA can automatically build
them from a given product or reactant state structure. In addition to the :ref:`CoreSettings`, such a job should be set up as:

.. code-block:: python

	job_type = find_ts
	initial_coordinates = [<coord_file_1>]
	commit_fwd = [...]
	commit_bwd = [...]
	max_moves = 5
	
See :ref:`CommitmentBasinDefinitions` for details on the `commit_fwd` and `commit_bwd` options.

Committor Analysis
~~~~~~~~~~~~~~~~~~

After completing aimless shooting, the next step is to obtain a reaction coordinate and verify it with committor analysis. Before running committor analysis, the user should call  likelihood maximization (lmax.py, see :ref:`LikelihoodMaximization`) and reaction coordinate evaluation (rc_eval.py, see :ref:`RCEval`). 

Then, 

Full Configuration Options
--------------------------

Here, the full list of valid configuration file options are documented (excluding the :ref`CoreSettings`, documented above.)

Batch Template Settings
~~~~~~~~~~~~~~~~~~~~~~~
        init_nodes: int = 1
        init_ppn: int = 1
        init_mem: str = '4000mb'
        init_walltime: str = '00:30:00'
        init_solver: str = 'sander'
        prod_nodes: int = 1
        prod_ppn: int = 8
        prod_mem: str = '4000mb'
        prod_walltime: str = '02:00:00'
        prod_solver: str = 'sander'

Path Settings
~~~~~~~~~~~~~
        path_to_input_files: str = sys.path[0] + '/data/input_files'    # todo: fix this for final publication
        path_to_templates: str = sys.path[0] + '/data/templates'

CV Settings
~~~~~~~~~~~
        cvs: typing.List[str] = ['']
        include_qdot: bool = True

Initial Coordinates
~~~~~~~~~~~~~~~~~~~
        initial_coordinates: typing.List[str] = ['']

.. _CommitmentBasinDefinitions

Commitment Basin Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        commit_fwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])
        commit_bwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])

.. _ReactionCoordinateDefinition

Reaction Coordinate Definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rc_definition: str = ''
        as_out_file: str = 'as_raw.out'
        rc_reduced_cvs: bool = True

.. _AimlessShootingSettings

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