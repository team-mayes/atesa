.. _ExampleStudy:

Example Study
============

This page will detail a complete workflow in ATESA for an example chemical reaction. Not every setting or option applied in this workflow will be appropriate for every application of ATESA; rather than being prescriptive, the purpose of this page is to illustrate a particular application of the software to help readers better grasp what running ATESA ia actually like.

Where appropriate, this page will include complete examples of input, template, and configuration files. **Not all of the contents of these files will be appropriate for your application**. Therefore, we strongly discourage you from copying-and-pasting content from this page unless you know what you're doing. The complete files are provided only for illustrative purposes.

Initial Setup and the Model
---------------------------

We assume here that ATESA has already been installed on your system. If that is not the case, you should instead start at the :ref:`Installation` page.

The model we will be working with here is the gas-phase decomposition of ethyl chlorosulfite into chloroethane and sulfur dioxide. The mechanism via S\ :sub:`N`\ i reconfiguration has been studied by `Schreiner et al. 1995 <https://pubs.acs.org/doi/pdf/10.1021/jo00086a041>`_. We chose this reaction because the small number of atoms involved facilitate quick simulations and easy visualizations, but ATESA has also been successfully applied for much larger systems, such as enzyme reaction modeling.

	.. figure:: _images/reaction_pathway.png

	The reaction pathway via S\ :sub:`N`\ i reconfiguration for gas phase ethyl chlorosulfite, per Schreiner *et al.* 1995, who demonstrated that this “frontside” attack (where the chlorine bonds to the same side of the carbon as the oxygen departs from) is energetically favorable compared to attacking from the opposite side. Teal: carbon; white: hydrogen; red: oxygen; yellow: sulfur; green: chlorine.

Setup of ATESA for a new system begins with obtaining initial coordinates through whatever means, such as download from a repository like the Protein Databank, or created bespoke using appropriate software. In either case, as with most molecular simulations, the system should first be minimized, heated to the desired temperature, and equilibrated in the desired relaxed state. Excellent tutorials for these steps (among others) using Amber are available at: https://ambermd.org/tutorials/

In this case, we produced initial coordinates using `Open Babel <http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html>`_ with the SMILES string "CCOS(=O)Cl", and then minimized, heated to 300 K, and briefly equilibrated the system in Amber.

Finding a Transition State
--------------------------

ATESA automates the discovery of suitable initial transition state models using gentle restraints to force reactions to take place, and then identifying the transition state along that forced pathway. The restraints are based only on user-defined definitions of the two stable states that the transition state connects. The complete configuration file used for this job was as follows (remember that many of these options are user- and model-specific!)::

	job_type = 'find_ts'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = False
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_find_ts'
	overwrite = True

	initial_coordinates = ['ethyl_chlorosulfite.inpcrd']

	commit_fwd = ([3,1,1],[5,5,2],[2.2,2.0,3.0],['gt','lt','gt'])
	commit_bwd = ([3,1,1],[5,5,2],[1.5,3.5,2.4],['lt','gt','lt'])

	max_moves = 5

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '04:00:00'
	prod_ppn = 1

As shown below, this job automatically finds a transition state very close to that proposed by Schreiner* et al. [3]

	.. figure:: _images/find_ts.png

	Definitions of stable states and initial and final structures from the example transition state search. The stable state definitions are read by inner index; for example, the first element of the definition of the “bwd” state is read as “the distance between atom 3 and atom 5 is less than (‘lt’) 1.5 Å”. Based on these definitions, the initial coordinates (at left) occupy the “bwd” state, and restraints are automatically constructed to build a transition state (at right) that has roughly equal probabilities of relaxing to either state. The narrow, transparent bonds in the transition state structure show the original topology of the model, for comparison.

Aimless Shooting
----------------

Once a model has been set up near the transition state, aimless shooting can proceed. In this case we used the transition state identified in the previous step as the initial coordinates, with 25 copies to speed up the sampling::

	job_type = 'aimless_shooting'
	topology = 'ethyl_chlorosulfite.prmtop'
	batch_system = 'slurm'
	restart = True
	working_directory = '/scratch/tburgin/ethyl_chlorosulfite_as'
	overwrite = True

	initial_coordinates = ['/scratch/tburgin/ethyl_chlorosulfite_find_ts/as_test/ethyl_chlorosulfite.inpcrd_0_ts_guess_97.rst7']
	degeneracy = 25

	commit_fwd = ([3,1,1],[5,5,2],[2.2,2.0,3.0],['gt','lt','gt'])
	commit_bwd = ([3,1,1],[5,5,2],[1.5,3.5,2.4],['lt','gt','lt'])

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '01:00:00'
	prod_ppn = 1
	
Note that the absence of any specific CVs in this file results in the default behavior, which is building CVs automatically based on the atoms indicated in the commitment definitions. In this case, ATESA derived 156 CVs to sample at each shooting point.

This job collected 13,052 shooting moves before terminating automatically using based on the :ref:`InformationError` termination criterion with the default settings. A visualization of the sampled data projected onto the three dimensions making up the commitment basins is shown here to help readers who may not be familiar with aimless shooting visualize the data, but this plot is not produced automatically by ATESA:

	.. figure:: _images/as_data.gif

Likelihood Maximization and rc_eval.py
--------------------------------------

After aimless shooting terminates, we pass the results to the auxiliary script :ref:`LikelihoodMaximization` in order to obtain a model reaction coordinate that describes the probable fate of a simulation beginning from a given set of initial conditions.

In order to minimize the influence of the initial coordinates chosen to begin aimless shooting with on this result, we use the largest decorrelated aimless shooting output file available, which in this case is named ``as_decorr_13000.out``. Decorrelated output files include only the shooting points after the point where all of the CVs have no correlation with their initial values for that thread with at least 95% confidence, or in other words when the autocorrelation of each CV is less than or equal to 1.96 / sqrt(n) for n shooting moves in the thread. These files are built automatically by ATESA when evaluating the information error termination criterion, but otherwise can be produced manually by running a repeat of the aimless shooting job with ``resample = True``.

To illustrate the flexibility of likelihood maximization, let's include one dimension that we want to be sure makes it into the final reaction coordinate. You might choose to do this if you have a particular interest in the relationship between a given dimension and the fate of a simulation. In this case, we'll choose the distance between the chlorine atom and the carbon to which it bonds in the product state. Consulting the automatically produced ``cvs.txt`` in the working directory (so in this case, located at ``/scratch/tburgin/ethyl_chlorosulfite_as/cvs.txt``) we see that this is CV12, so we call likelihood maximization as::

	lmax.py -i /scratch/tburgin/ethyl_chlorosulfite_as/as_decorr_13000.out --automagic --plots -f 12

The ``--plots`` option produces the sigmoid committor plot and, when automatic is used as is the case here, the two-line test plot shown below. The good relationship between the modeled and ideal committor sigmoids is a necessary, but not a sufficient, condition for a good reaction coordinate:

	.. figure:: _images/lmax.png
	
After selecting a reaction coordinate, we need to call the auxiliary script :ref:`RCEval`. This will build a file named ``rc.out`` in the working directory, which we need in the next step::

	rc_eval.py /scratch/tburgin/ethyl_chlorosulfite_as/ -1.352-2.766*CV12+3.572*CV22+2.985*CV4

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
	rc_threshold = 0.005
	committor_analysis_n = 20

	path_to_input_files = '/home/tburgin/ethyl_chlorosulfite/input_files'
	path_to_templates = '/home/tburgin/ethyl_chlorosulfite/templates'

	prod_walltime = '01:00:00'
	prod_ppn = 1
	
The use of ``as_settings_file`` to point to the ``settings.pkl`` file produced during aimless shooting ensures that the same commitment basin and CV definitions are used. The next block of options specifies how committor analysis will be carried out: all of the shooting points identified in ``/scratch/tburgin/ethyl_chlorosulfite_as/rc.out`` (the file produced just before by ``rc_eval.py``) as having a reaction coordinate absolute value of less than or equal to the threshold value of 0.005 will be used to seed 20 individual committor analysis simulations, the results of which when taken together make up the output of committor analysis.

Plotting the contents of the output file produced by this job (``/scratch/tburgin/ethyl_chlorosulfite_as/committor_analysis/committor_analysis.out``) as a histogram, we see that it is roughly even and centered near 0.5, which affirms that our reaction coordinate is a good one.

	.. figure:: _images/ethyl_chlorosulfite_comana.png

Pathway-restrained Umbrella Sampling
------------------------------------