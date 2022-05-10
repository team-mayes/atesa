<img src="docs/_images/atesa_logo.png" alt="ATESA_logo" width="400"/>

[//]: # (Badges)
[![Tests](https://github.com/team-mayes/atesa/actions/workflows/CI.yml/badge.svg)](https://github.com/team-mayes/atesa/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/team-mayes/atesa/branch/master/graph/badge.svg)](https://codecov.io/gh/team-mayes/atesa/branch/master)
[![Documentation Status](https://readthedocs.org/projects/atesa/badge/?version=latest)](https://atesa.readthedocs.io/en/latest/?badge=latest)


A Python program for automating transition path sampling with aimless shooting, suitable for experts and novices alike.

Full documentation available [here](https://atesa.readthedocs.io/en/latest/).

ATESA automates a particular Transition Path Sampling (TPS) workflow that uses the flexible-length aimless shooting algorithm of [Mullen *et al.* 2015](http://doi.org/10.1021/acs.jctc.5b00032). ATESA interacts directly with a batch system or job manager to dynamically submit, track, and interpret various simulation and analysis jobs based on one or more initial structures provided to it. The flexible-length implementation periodically checks simulations for commitment to user-defined reactant and product states in order to maximize the acceptance ratio and minimize wasted computational resources.

ATESA implements automation for obtaining a suitable initial transition state, flexible-length aimless shooting, inertial likelihood maximization, committor analysis, umbrella sampling (and analysis with the Multistate Bennett Acceptance Ratio), and equilibrium path sampling. These components constitute a near-complete automation of the workflow between identifying the reaction of interest, and obtaining, validating, and analyzing the energy profile along an unbiased and *bona fide* reaction coordinate that describes it.

At present, ATESA only supports simulations with Amber, and TORQUE/PBS or Slurm batch schedulers. If you are interested in using ATESA with another simulation engine or batch scheduler, please raise an "enhancement" issue describing your needs.

### Copyright

Copyright (c) 2019, Tucker Burgin


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.

Special thanks to Samuel Ellis and the Molecular Sciences Software Institute (MolSSI).
