.. _ExampleStudy:

Example Study
============

This page will detail a complete workflow in ATESA for an example chemical reaction. Not every setting or option applied in this workflow will be appropriate for every application of ATESA; rather than being prescriptive, the purpose of this page is to illustrate a particular application of the software to help readers better grasp what running ATESA ia actually like.

Where appropriate, this page will include complete examples of input, template, and configuration files. **Not all of the contents of these files will be appropriate for your application**. Therefore, we strongly discourage you from copying-and-pasting content from this page unless you know what you're doing. The complete files are provided only for illustrative purposes.

Initial Setup and the Model
---------------------------

We assume here that ATESA has already been installed on your system. If that is not the case, you should instead start at the :ref:`Installation` page.

The model we will be working with here is the gas-phase decomposition of ethyl chlorosulfite into chloroethane and sulfur dioxide. The mechanism via S\ :sub:`N`\ i reconfiguration has been studied by `Schreiner *et al.* 1995 <https://pubs.acs.org/doi/pdf/10.1021/jo00086a041>`_.

	.. figure:: _images/reaction_pathway.png

	The reaction pathway via S\ :sub:`N`\ i reconfiguration for gas phase ethyl chlorosulfite, per Schreiner *et al.* 1995, who demonstrated that this “frontside” attack (where the chlorine bonds to the same side of the carbon as the oxygen departs from) is energetically favorable compared to attacking from the opposite side. Teal: carbon; white: hydrogen; red: oxygen; yellow: sulfur; green: chlorine.

Setup of ATESA for a new system begins with .......

Finding a Transition State
--------------------------

Aimless Shooting
----------------

Likelihood Maximization and rc_eval.py
--------------------------------------

Committor Analysis
------------------

Pathway-restrained Umbrella Sampling
------------------------------------