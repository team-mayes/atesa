.. atesa documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ATESA's documentation!
=========================================================

Python program for automating transition path sampling with aimless shooting using Amber.

ATESA automates a particular Transition Path Sampling (TPS) workflow that uses the flexible length aimless shooting algorithm of `Mullen et al. 2015 <https://pubs.acs.org/doi/10.1021/acs.jctc.5b00032>`_ as its main workhorse. TPS is an approach to obtaining and analyzing a reaction coordinate that describes a given chemical transformation on a computationally tractable timescale. An excellent primer on the method can be found in `Beckham and Peters, 2010 <https://pubs.acs.org/doi/abs/10.1021/bk-2010-1052.ch013>`_. ATESA interacts directly with a batch system (either TORQUE (PBS) or Slurm are supported) to dynamically submit, track, and interpret various simulation and analysis jobs based on one or more initial structures provided to it. It employs a variable length approach that periodically checks simulations for commitment to user-defined reactant and product states in order to maximize the acceptance ratio and minimize wasted computational resources.

ATESA contains methods to automate initial transition state guessing, flexible length aimless shooting, intertial likelihood maximization, committor analysis, and equilibrium path sampling. In combination, these components constitute a near-complete automation of the workflow between identifying the reaction of interest, and obtaining, validating, and analyzing the energy profile along a *bona fide* reaction coordinate that describes it.

.. toctree::
   :maxdepth: 2
   :caption: For Users:

   getting_started
   
.. toctree::
   :maxdepth: 2
   :caption: API:
   
   apidoc/atesa_v2.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
