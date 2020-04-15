.. _Troubleshooting:

Troubleshooting
===============

Although every effort has been made to provide helpful error messages to correspond to every foreseeable issue, unfortunately not every issue is, in fact, foreseeable. This page is designed to help users troubleshoot should they run into issues not accompanied by a sufficiently helpful error message. Just because your transition path sampling shooting is aimless, doesn't mean your troubleshooting has to be!

Issues running ATESA will fall into one of three categories, based on where in the pipeline the issue arises:

#. The local shell or the batch system

#. Amber

#. ATESA itself

A description of how to identify the category of error being encountered, and suggestions to resolve the issue, are provided. If this page does not help you resolve your issue, please raise it on `our GitHub page <https://github.com/team-mayes/atesa>`_.

Issues with the local shell or the batch system
-----------------------------------------------

Errors encountered by these systems will appear in one of two places. If the shell encounters an error, it will be reported directly in the command line after ATESA is called. This likely indicates an issue with the installation of the software, or of Python itself. Command line errors will also arise if the batch system does not recognize the formatting of the batch files submitted by ATESA; double-check that your batch templates are of the correct format. Also, note that if ATESA is called from inside a batch job (as is recommended when performing many simulations), the command line output will be piped to the batch output file for that job.

Errors encountered by the batch system will also appear in batch output files, but in this case they will appear in the individual batch job output files in the working directory. These errors are likely to do with the formatting of the batch template files used by ATESA to create its job batch files. For example, if the user has removed the line in the batch template indicating the nodes to be assigned to the job on a Slurm system, ATESA will not complain (as omitting variables in the templates is supported in general), but the batch system will (because this particular variable is mandatory in Slurm). Other batch system errors include those not specific to ATESA, such as insufficient funds on the user account; to check whether a given error is of this type, try submitting a job that does not involve ATESA.

Issues with Amber
-----------------

Amber has its own output file for each job, usually given by the suffix “.out”  or “.mdout” in the working directory (note: ATESA may produce output files named with the “.out” suffix in this directory, as well.) If an error arises here, it likely indicates an issue with the formatting of the input coordinate or topology files, or else with the input file passed to Amber. Another common source of error when using QM/MM models is the absence of the necessary quantum mechanics files (located in $(AMBERHOME)/dat/slko/) for the chosen level of theory. If these are missing they may need to be installed for you by a user with root access. ATESA has been successfully tested with Amber 14, 16, and 18, but earlier versions are not guaranteed to be compatible.

Issues with ATESA
-----------------

Errors encountered by ATESA itself will be outputted into the command line in the event that the program has been called directly, or else into the batch output file in the event that it has been called in a batch job. Hopefully, most errors of this type will be self-explanatory. There may, however, be unforeseen errors that are not handled in a useful way. Errors of this sort could have to do with improperly installed dependencies, especially pytraj. One possible troubleshooting option is to retry the job with as many default settings as possible to see if the error persists — if it does not, add on custom settings one-by-one (starting with those you are least confident in) until the culprit is found.

If ATESA outputs an error message that confuses you, especially if it happens very early on in the job, double-check that the simulations that have been run so far (if there are any) appear to have behaved properly.

Although the user is of course permitted under the license to attempt to debug the code themselves (and to use the modified code in whatever manner they see fit), no one finds it easy to read someone else’s code. If you encounter an issue with ATESA that you cannot resolve, I encourage you to raise an issue on `our GitHub page <https://github.com/team-mayes/atesa>`_. And if you encounter an error that you did resolve by modifying the code, please submit a pull request so that everyone can benefit!

Not happy with the results?
---------------------------

One final type of issue is the case where the software has functioned as intended, but has not provided results that satisfy you. Perhaps your aimless shooting acceptance ratios are very poor, or even zero, or perhaps ``lmax.py`` was unable to confidently determine a reaction coordinate from the aimless shooting output file. In such cases, the user is encouraged to look to the chemistry. A chronically poor acceptance ratio during aimless shooting, even across many degenerate threads or slightly different input coordinate files, probably indicates that your initial coordinates are not as close to the reaction separatrix as you might have hoped. Even if they are very close, an extremely steep energetic landscape makes aimless shooting difficult, and depending on the context of your study may indicate an incorrect putative reaction mechanism. You should also consider making the *max_dt* and/or *min_dt* settings smaller; although larger values can explore phase space more efficiently, it also increases the risk of a thread straying too far from the separatrix and being unable to climb back up, resulting in poorer acceptance ratios after that step. This effect can be mitigated by setting *always_new* to “True” (which is the default).

If acceptance ratios are good but ``lmax.py`` is struggling to identify a strong reaction coordinate, consider resampling your existing simulation files (using the *resample* option) with more candidate values included; it’s possible that you left something important out. Remember that it is impossible to include too many candidate CVs (insofar as the size of the resulting aimless shooting output file is not prohibitive in terms of I/O or the speed of ``lmax.py``), but it is easy to include too few.

If you have good ``lmax.py`` results (an apparently good fit to the committor probability sigmoid) but a poor committor analysis histogram (peaked near forward commitment probabilities of 1, 0, or both) in spite of including an ample number of committor analysis shooting points (100-200 with 10 simulations per point is appropriate), this likely indicates either that you have underfitted or overfitted your reaction coordinate (too few or too many CVs included), or that you have not performed enough sampling with aimless shooting (your sample of shooting points is not representative of the underlying population). If you suspect the latter issue, try setting the information error threshold lower (see :ref:`AimlessShootingSettings`) and resubmitting the aimless shooting job with *restart = True*.