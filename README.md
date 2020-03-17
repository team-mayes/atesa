atesa
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/team-mayes/atesa_v2.png)](https://travis-ci.org/team-mayes/atesa)
[![codecov](https://codecov.io/gh/team-mayes/atesa_v2/branch/master/graph/badge.svg)](https://codecov.io/gh/team-mayes/atesa/branch/master)

A flexible and extensible program for automating transition path sampling with aimless shooting.

ATESA automates a particular Transition Path Sampling (TPS) workflow that uses the flexible length aimless shooting algorithm of [Mullen *et al.* 2015](http://doi.org/10.1021/acs.jctc.5b00032) as its main workhorse. ATESA interacts directly with a batch system or job manager to dynamically submit, track, and interpret various simulation and analysis jobs based on one or more initial structures provided to it. It employs a variable length approach that periodically checks simulations for commitment to user-defined reactant and product states in order to maximize the acceptance ratio and minimize wasted computational resources.

ATESA contains methods to automate initial transition state guessing, flexible length aimless shooting, intertial likelihood maximization, committor analysis, and equilibrium path sampling. In combination, these components constitute a near-complete automation of the workflow between identifying the reaction of interest, and obtaining, validating, and analyzing the energy profile along a *bona fide* reaction coordinate that describes it.

### Copyright

Copyright (c) 2019, Tucker Burgin


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
