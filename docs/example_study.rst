.. _ExampleStudy:

Example Study
============

This page will detail a complete workflow in ATESA for an example chemical reaction. Not every setting or option applied in this workflow will be appropriate for every application of ATESA; rather than being prescriptive, the purpose of this page is to illustrate a particular application of the software to help readers better grasp what running ATESA is actually like. If you're instead looking for an introduction to the theory underlying any of these steps, check out the :ref:`TheoryAndDefinitions` page.

Where appropriate, we have included complete examples of input, template, and configuration files as well as key results for each step. They can be found in the 'examples' folder downloaded with ATESA. **Not all of the contents of these files will be appropriate for your application**. You will almost certainly need to make changes to one or more of these files to suit your particular model and computing cluster environment. Also note that because of the non-deterministic nature of these methods, you should not expect to get exactly identical results if you follow along with these steps yourself. Remember that configuration files are used by submitting an ATESA job as follows::

	atesa example.config

Finally, we have included the approximate computational resources consumed during each step. As with all molecular simulations, your results may vary widely based on the details of your simulation and the computational environment, so be cautious in applying these estimates. In our case, all performance statistics are based on simulations on the San Diego Supercomputer Center's Comet platform, using single cores (except where otherwise specified) on compute nodes with the following characteristics::

	Intel Haswell Standard Compute Nodes
	Clock speed	2.5 GHz
	Cores/node	24
	DRAM/node	128 GB
	SSD memory/node	320 GB

Initial Setup and the Model
---------------------------

We assume here that ATESA has already been installed on your system. If that is not the case, you should instead start at the :ref:`Installation` page.

The model we will be working with here is the gas-phase decomposition of ethyl chlorosulfite into chloroethane and sulfur dioxide. The mechanism via S\ :sub:`N`\ i reconfiguration has been studied by `Schreiner et al. 1995 <https://pubs.acs.org/doi/pdf/10.1021/jo00086a041>`_. We chose this reaction because the small number of atoms involved facilitate quick simulations and easy visualizations, but ATESA has also been successfully applied for much larger systems, such as enzyme reaction modeling.

	.. figure:: _images/reaction_pathway.png

	The reaction pathway via S\ :sub:`N`\ i reconfiguration for gas phase ethyl chlorosulfite, per Schreiner *et al.* 1995, who demonstrated that this “frontside” attack (where the chlorine bonds to the same side of the carbon as the oxygen departs from) is energetically favorable compared to attacking from the opposite side. Teal: carbon; white: hydrogen; red: oxygen; yellow: sulfur; green: chlorine.

Setup of ATESA for a new system begins with obtaining initial coordinates through whatever means, such as download from a repository like the Protein Databank, or created bespoke using appropriate software. In either case, as with most molecular simulations, the system should first be minimized, heated to the desired temperature, and equilibrated in the desired relaxed state. Excellent tutorials for these steps (among others) using Amber are available `here <https://ambermd.org/tutorials/basic/tutorial0/index.htm>`_. Fluency in basic molecular simulations is assumed of users of ATESA, so be sure you're comfortable before moving forward!

In this case, we produced initial coordinates using `OpenBabel <http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html>`_ with the SMILES string "CCOS(=O)Cl", and then equilibrated the system at 300 K over 10 ps in Amber (a very short time, appropriate for such a small system). Because we are modeling a chemical reaction, this equilibration is done with a hybrid QM/MM model. The topology and equilibration input files and the resulting coordinates can be found in 'examples/model_setup'.

Finding a Transition State
--------------------------

ATESA automates the discovery of suitable initial transition state models using gentle restraints to force reactions to take place, and then identifying the transition state(s) along that forced pathway. The restraints are based only on user-defined definitions of the two stable states that the transition state connects. The complete configuration file used for this job was as follows (remember that many of these options are user- and model-specific!)::

	job_type = 'find_ts'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = False
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_find_ts'
	overwrite = False

	initial_coordinates = ['ethyl_chlorosulfite.inpcrd']

	commit_fwd = ([3,1,1],[5,5,2],[2.2,2.0,3.0],['gt','lt','gt'])
	commit_bwd = ([3,1,1],[5,5,2],[1.5,3.5,2.4],['lt','gt','lt'])

	max_moves = 5

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '04:00:00'
	prod_ppn = 1

Based on the five aimless shooting moves with which each candidate transition state frame from the forced trajectory was tested, two frames were selected as suitable aimless shooting initial coordinate files, as indicated in ATESA's output following this job. The coordinates for those frames can be found in 'examples/find_ts'. As shown below, these transition states (they are nearly identical) are very close to the one proposed by Schreiner *et al.* Under our test conditions, this job took 18 minutes to complete.

	.. figure:: _images/find_ts.png

	Definitions of stable states and initial and final structures from the example transition state search. The stable state definitions are read by inner index; for example, the first element of the definition of the “bwd” state is read as “the distance between atom 3 and atom 5 is less than (‘lt’) 1.5 Å”. Based on these definitions, the initial coordinates (at left) occupy the “bwd” state, and restraints are automatically constructed to build a putative transition state (at right) that has significantly non-zero probabilities of relaxing to either state. The narrow, transparent bonds in the transition state structure show the original topology of the model, for comparison.

Aimless Shooting
----------------

Once a model has been set up near the transition state, aimless shooting can proceed. In this case we used the two transition state models identified in the previous step as the initial coordinates, with 12 copies ('degeneracy = 12') each to speed up the sampling::

	job_type = 'aimless_shooting'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = False
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_as'
	overwrite = False

	initial_coordinates = ['/scratch/tburgin/ethyl_chlorosulfite_find_ts/as_test/ethyl_chlorosulfite.inpcrd_0_ts_guess_134.rst7', '/scratch/tburgin/ethyl_chlorosulfite_find_ts/as_test/ethyl_chlorosulfite.inpcrd_0_ts_guess_136.rst7']
	degeneracy = 12

	commit_fwd = ([3,1,1],[5,5,2],[2.2,2.0,3.0],['gt','lt','gt'])
	commit_bwd = ([3,1,1],[5,5,2],[1.5,3.5,2.4],['lt','gt','lt'])
	
	information_error_freq = 2500

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '00:30:00'
	prod_ppn = 1
	
Note that the absence of any specific CVs in this file results in the default behavior, which is building CVs automatically based on the atoms indicated in the commitment definitions. In this case, ATESA derived 156 CVs to sample at each shooting point. We also set the number of steps between assessments of the information error termination criterion 10 times higher than the default since our system is very small and we'll accumulate hundreds of simulations very rapidly. Similarly, we set a short walltime and allocate only a single core to the production simulations to reflect their low computational requirements.

Note that because the simulations for this job are so short, it is best to take advantage of ATESA's built-in multiprocessing support for this task. The optimal number of cores to allocate will depend greatly on your platform, but using roughly as many cores as you have aimless shooting threads is a reasonable starting point. In this case, we selected 24 aimless shooting threads (12 for each initial coordinate file) to make optimal use of one 24-core node.

During our testing, this job collected data at the rate of approximately 410 shooting moves per hour; the bottleneck in this case was waiting for the batch system to allocate resources for individual simulations. For longer individual simulations, this overhead will be greatly diminished.

This job collected 15,142 shooting moves before terminating automatically using based on the :ref:`InformationError` termination criterion with the default settings. An average acceptance ratio of 31.76% (per ``status.txt`` in the working directory) reflects a very healthy level of efficiency (10-30% is average). ATESA also automatically generates a version of the aimless shooting file that has been decorrelated from the initial state(s) when assessing the information error termination criterion. Both the raw and decorrelated output files can be found in 'examples/aimless_shooting'.

..
	A visualization of the sampled data projected onto the three dimensions making up the reaction coordinate produced in the next step is shown here to help readers who may not be familiar with aimless shooting visualize the data, but this plot is not produced automatically by ATESA. The raw aimless shooting output file can be found in 'examples/aimless_shooting'.

..
	.. figure:: _images/as_data.gif

..	
	A visualization of the aimless shooting data. Each point represents a single shooting move, colored according to its fate: blue for reactant, orange for products.

Likelihood Maximization and rc_eval.py
--------------------------------------

After aimless shooting terminates, we pass the results to the auxiliary script :ref:`LikelihoodMaximization` in order to obtain a model reaction coordinate that describes the probable fate of a simulation beginning from a given set of initial conditions.

In order to minimize the influence of the initial coordinates chosen to begin aimless shooting with on this result, we use the largest decorrelated aimless shooting output file available, which in this case is named ``as_decorr_30000.out``. Decorrelated output files include only the shooting points after the point where all of the CVs have no correlation with their initial values for that thread with at least 95% confidence, or in other words when the autocorrelation of each CV is less than or equal to 1.96 / sqrt(n) for n shooting moves in the thread. These files are built automatically by ATESA when evaluating the information error termination criterion, but otherwise can be produced manually by running a repeat of the aimless shooting job with ``resample = True``.

	lmax.py -i /scratch/tburgin/ethyl_chlorosulfite_as/as_decorr_30000.out --two_line_test --plots

The ``--plots`` option produces the sigmoid committor plot ((b), at right) and, when automatic is used as is the case here, the two-line test plot (see :ref:`LikelihoodMaximization`). The good relationship between the modeled and ideal committor sigmoids is a necessary, but not a sufficient, condition for a good reaction coordinate:

	.. figure:: _images/lmax.png
	
ATESA does not automatically generate the figure in panel (a), but it's a useful way to visualize the reaction coordinate. Although the chosen dimensions are non-obvious, they nonetheless have a strong rationale: the first term describes the relative closeness of the reactive carbon to either of the other atoms it may bond to, while the two angle terms describe the direction of the partially charged face of the carbon atom, using the sulfur atom as a reference point.

After selecting a reaction coordinate, we need to call the auxiliary script :ref:`RCEval`. This will build a file named ``rc.out`` in the working directory, which we need in the next step::

	rc_eval.py /scratch/tburgin/ethyl_chlorosulfite_as/ -5.68+10.57*CV156-1.86*CV30+1.68*CV29

Committor Analysis
------------------

Having obtained what appears to be a suitable reaction coordinate, the final step in validating it before using it to analyze the energy profile is to perform committor analysis. By performing new simulations (*i.e.*, simulations whose results were not included in the likelihood maximization training data) with initial reaction coordinate values of approximately zero, we can confirm that the reaction coordinate is an accurate descriptor of the transition state (at least within the context of our particular simulation conditions).

Committor analysis is again called through the main ATESA script. Our complete configuration file for this job is as follows::

	job_type = 'committor_analysis'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = False
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_as/committor_analysis'
	overwrite = True

	as_settings_file = '/scratch/tburgin/ethyl_chlorosulfite_as/settings.pkl'

	committor_analysis_use_rc_out = True
	path_to_rc_out = '/scratch/tburgin/ethyl_chlorosulfite_as/rc.out'
	rc_threshold = 0.0025
	committor_analysis_n = 20

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '01:00:00'
	prod_ppn = 1
	
The use of ``as_settings_file`` to point to the ``settings.pkl`` file produced during aimless shooting ensures that the same commitment basin and CV definitions are used. The next block of options specifies how committor analysis will be carried out: all of the shooting points identified in ``/scratch/tburgin/ethyl_chlorosulfite_as/rc.out`` (the file produced just before by ``rc_eval.py``) as having a reaction coordinate absolute value of less than or equal to the threshold value of 0.005 will be used to seed 20 individual committor analysis simulations, the results of which when taken together make up the output of committor analysis.

Plotting the contents of the output file produced by this job (``/scratch/tburgin/ethyl_chlorosulfite_as/committor_analysis/committor_analysis.out``) as a histogram, we see that it is roughly even and centered towards the middle, which affirms that our reaction coordinate is a reasonably good one.

	.. figure:: _images/ethyl_chlorosulfite_comana.png

Umbrella Sampling
------------------------------------

Finally, we're ready to evaluate the energy profile along our reaction coordinate. ATESA features two separate job types for this purpose: equilibrium path sampling, and umbrella sampling. If you have access to it, the latter is usually strongly preferable, so we'll focus on that here.

First, we need to identify the appropriate minimum and maximum RC values to sample over. ATESA has a built-in script to facilitate this: ``rc_eval.py`` (see :ref:`RCEval` for more details)::

	rc_eval.py /scratch/tburgin/ethyl_chlorosulfite_as/ -1.600+2.053*CV156+0.576*CV1+0.655*CV22 /scratch/tburgin/ethyl_chlorosulfite_as_/as_decorr_15000.out True
	
This script completes in a matter of seconds, and simply returns the ending RC values observed at the ends of both trajectories in an accepted shooting move. In other words, it approximates the RC values of the stable states. In this case, it returned::

	 Shooting move name: ethyl_chlorosulfite.inpcrd_0_ts_guess_134.rst7_0_663_init.rst7
	 extrema: [8.196560999999999, -3.9342559999999995]

The first line indicates the shooting move that was selected, and the second indicates the RC extrema. To be sure we include the full stable state energy basins, we'll extend our umbrella sampling past these values by about 10%.

Finally, we need to select the appropriate spacing and restraint weights for our umbrella sampling windows. Since the applied restraints are harmonic, the expected width of the sampled distribution is be approximately proportional to the inverse square root of the restraint weight. We'll select a restraint weight of 10 times kT (the Boltzmann constant times the temperature, equal to about 0.6 kcal/mol at 300 K), which should produce a distribution of approximate width 1.5 along the dimensionless reaction coordinate, and space the windows 0.5 apart to ensure substantial overlap between adjacent windows::

	job_type = 'umbrella_sampling'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = False
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_as/umbrella_sampling'
	overwrite = True

	initial_coordinates = ['fwd.nc','bwd.nc']

	rc_definition = '-5.68 + 10.57*CV156 - 1.86*CV30 + 1.68*CV29'

	as_out_file = '/scratch/tburgin/ethyl_chlorosulfite_as/as_decorr_30000.out'
	as_settings_file = '/scratch/tburgin/ethyl_chlorosulfite_as/settings.pkl'

	us_rc_step = 0.5
	us_restraint = 6
	us_rc_min = -4.2
	us_rc_max = 9

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '04:00:00'
	prod_ppn = 1
	
The choice of a stepsize of 0.25 and a restraint weight of 5 is not arbitrary; these values must be carefully matched so as to ensure overlap among all adjacent windows. It is usually wise to run a pilot study with only a single window to verify the approximate width of the sampling histogram for your particular settings (it is safe to assume that each window will be approximately even in width, though they may be shifted from their centers somewhat).

This job produces a large number of output files named with the suffix "_us.out" in the working directory. When it's finished, we can use the auxiliary script ``mbar.py`` to analyze it::

	mbar.py --decorr -k 5
	
Here we use the `--decorr` flag to specify that we have not checked the data for decorrelation or equilibration, so pyMBAR will do that work for us. We also set `-k 5` to indicate that the umbrella sampling restraint is 5 kcal/mol. After a mean-value plot (see :ref:`UmbrellaSamplingTroubleshooting` for more information), ``mbar.py`` produces the data histograms and free energy profile shown below (though I have of course added the reference activation energy after the fact):

	.. figure:: _images/umbrella_sampling.png
	
That the histograms overlap everywhere leaving no gaps is a necessary-but-not-sufficient condition for an acceptable dataset. Comparison of the activation energy to the theoretical value goes a long way towards validating this result.
