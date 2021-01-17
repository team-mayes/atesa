.. atesa documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ATESA: Aimless Transition Ensemble Sampling and Analysis
=========================================================

A Python program for automating transition path sampling with aimless shooting, suitable for experts and novices alike.

ATESA automates a particular Transition Path Sampling (TPS) workflow that uses the flexible-length aimless shooting algorithm of `Mullen et al. 2015 <http://doi.org/10.1021/acs.jctc.5b00032>`_. ATESA interacts directly with a batch system or job manager to dynamically submit, track, and interpret various simulation and analysis jobs based on one or more initial structures provided to it. The flexible-length implementation periodically checks simulations for commitment to user-defined reactant and product states in order to maximize the acceptance ratio and minimize wasted computational resources.

ATESA implements automation for obtaining a suitable initial transition state, flexible-length aimless shooting, inertial likelihood maximization, committor analysis, umbrella sampling (and analysis with the Multistate Bennett Acceptance Ratio), and equilibrium path sampling. These components constitute a near-complete automation of the workflow between identifying the reaction of interest, and obtaining, validating, and analyzing the energy profile along an unbiased and *bona fide* reaction coordinate that describes it.

The batch systems and molecular simulations packages currently supported by ATESA (please raise an issue with the "enhancement" label on `our GitHub page <https://github.com/team-mayes/atesa>`_ if you'd like to see something added to this list!):

**Batch Systems**
	* Slurm
	* PBS/TORQUE

**Simulations Packages**
	* Amber

.. toctree::
   :maxdepth: 2
   :caption: For Users:

   theory_of_operation
   getting_started
   the_config_file
   auxiliary_scripts
   example_study
   on_termination_criteria
   troubleshooting
   
.. toctree::
   :maxdepth: 1
   :caption: API:
   
   apidoc/atesa.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
