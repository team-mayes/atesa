A Quick Theory Primer
=====================

This page is intended to briefly introduce readers to the theory of transition path sampling (TPS) with aimless shooting, as well as to how ATESA works. A thorough explanation of TPS as a technique and its relationship to other methods is beyond the scope of this document, and the reader is instead directed to `Beckham and Peters, 2010 <https://pubs.acs.org/doi/abs/10.1021/bk-2010-1052.ch013>`_.

Transition Path Sampling
------------------------

Transition path sampling refers to any method of rare event sampling that aims to characterize the ensemble of "transition paths" (that is, paths through phase space that connect one state to another) as a proxy for characterizing the event as a generalized whole.

Aimless Shooting
----------------

Aimless shooting is a transition path sampling method for performing efficient, unbiased sampling of the region of phase space corresponding to a putative transition state. Because by definition the transition state is a local maximum in energy along at least one dimension, this region is difficult to sample using conventional simulations (that is, a transition is a rare event). The aimless shooting approach is to leverage one or more putative or “guess” transition state structures (which are obtained by other methods as a prerequisite to beginning aimless shooting (though ATESA is equipped with a tool to help do so)) as a seed, which will be “aimlessly” “shot” through phase space using unbiased initial velocities chosen from the appropriate Boltzmann distribution. The resulting trajectory is at the same time also simulated in reverse (using initial velocities of opposite direction and equal magnitude), and if after simulations the two trajectories converge to different energetic basins (one products, one reactants), then the reactive trajectory connecting them is considered a success. New starting points are daisy-chained from older successful ones by taking an early frame from the reactive trajectory as the initial coordinates, and in this way it is ensured that sampling remains nearby the transition state separatrix (that is, the surface in phase space that divides the products from the reactants with equal committor probability in either direction).

ATESA
-----

ATESA automates the aimless shooting process with a system of independent “threads” representing one particular path in the search through phase space. A thread has a given set of initial coordinates, which it repeatedly “shoots” until it finds a successful reactive trajectory, at which point it picks a new shooting point on the reactive trajectory and continues. Because threads run entirely in parallel, aimless shooting with ATESA scales perfectly so long as sufficient computational resources are available.

ATESA also features a suite of analysis and utility tools that run in much the same fashion. Most importantly, once aimless shooting has been completed (see On Termination Criteria), ATESA automates likelihood maximization to mine the data for a reaction coordinate that describes the transition path, committor analysis to verify that reaction coordinate, and equilibrium path sampling to obtain the free energy profile along it.

Likelihood Maximization
-----------------------

The product of aimless shooting is a (large) set of combined variable (CV) values paired with corresponding commitment basins (products or reactants). In order to convert this information into a usable form, the method of likelihood maximization can be used to select a model that describes the reaction progress in terms of relatively few parameters. ATESA uses the intertial likelihood maximization procedure first published in `Peters 2012 <https://doi.org/10.1016/j.cplett.2012.10.051>`_.

Committor Analysis
------------------

Once a reaction coordinate has been obtained, it should be verified using new, unbiased simulations that were not included in the model training dataset. The method of committor analysis is to simply select a large number (hundreds) of initial coordinates with reaction coordinate values very close to zero (the predicted transition state) and run several unbiased simulations starting from each of them to verify that they are as likely on average to proceed towards the reactants as towards the products.

Equilibrium Path Sampling
-------------------------

The free energy profile along the reaction coordinate can be determined by several methods, but the most robust and general of them is equilibrium path sampling, wherein the reaction coordinate is divided into bins and the unbiased distribution of simulations within those bins is converted directly into free energy. ATESA automates collection of equilibrium path sampling data from an arbitrary array of initial coordinates, filling in gaps automatically using the tails of simulations from adjacent windows.