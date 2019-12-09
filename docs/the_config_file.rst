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