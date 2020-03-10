.. _OnTerminationCriteria:

On Termination Criteria
=======================

.. _TheProblem:

The Problem
-----------

Most any molecular system that is interesting enough to bother simulating operates within an extremely multidimensional phase space. In the most general case, a complete exploration of this phase space requires not only that the full range of each dimension be explored individually, but also the full range of combinations of dimensions; that is, the complexity of phase spaces scales exponentially with the number of dimensions. Thus, perfect sampling of even small sections of phase space is most often totally impractical.

If we know that our simulation will never fully explore phase space, how do we know whether we've sampled enough to trust the results? This issue has been aptly dubbed "the sampling problem", although given its centrality and importance to the field of molecular simulations, one might simply call it "The Problem".

The only "solution" to The Problem is to circumvent it by relying on free energy to "weight" regions of phase space by their relative importance in determining relevant features of the system. That is, regions of phase space with relatively lower free energy play a proportionally larger role in controlling the system's behavior, so distributing sampling according to free energy is the most efficient means of arriving at a good approximation of the total behavior of the system. Of course, this is not really a solution at all -- that is, it does not address the central issue of identifying when our sampling is sufficient, or even whether we're sampling all of the important regions of phase space in the first place. This can only be done post-hoc, by verifying the putative results of simulations in some other way, and then retrospectively affirming that the sampling had been sufficient.

Relevance to Aimless Shooting
-----------------------------

The process of aimless shooting can be thought of as constrained sampling -- specifically, in the region of phase space near enough to the reaction separatrix that both the product and reactant states are accessible. The aimless shooting algorithm automatically distributes its sampling according to the acceptance probability (that is, the likelihood of forming a complete reactive trajectory from a given initial set of coordinates shot "forward" and "backward") and the free energy (simulations are more likely to proceed towards regions of lower energy, which in aggregate causes successive shooting moves to become centered at local minima along the separatrix), but this is no remedy to The Problem; there remains no means of determining when sufficient sampling has been performed. But we can't go on sampling forever, so how do we choose when to stop?

Choosing When to Stop with Committor Analysis
---------------------------------------------

To date, the most common method of determining when to stop aimless shooting has been :ref:`CommittorAnalysis`. This requires that the aimless shooting data be mined for a reaction coordinate that describes the reaction progress (usually via :ref:`LikelihoodMaximizationTheory`), and then that that reaction coordinate be verified by showing that unbiased simulations beginning from its center are (approximately) equally as likely to proceed to the reactants as to the products. This method is invaluable and should certainly be used to affirm the correctness of the final reaction coordinate, but as a termination criterion it suffers from two issues:

#. It requires an extensive amount of additional simulation, which is expensive and time consuming; and

#. It provides no feedback about whether further sampling should be expected to improve the result.

In practice, this means that the researcher performing aimless shooting must either perform rounds of committor analysis many times throughout the study, or else collect such an excess of data that any unsatisfactory result in committor analysis is assumed to be due to issues elsewhere. Either option wastes time, both computational and real.

.. _InformationError:

Information Error
-----------------

ATESA introduces a new method of addressing the two shortcomings identified above, based on an assessment of the error in the likelihood maximization procedure. Specifically, ATESA measures the mean value of the parametric standard errors from the Godambe information matrix (a more general estimator of the Fisher information matrix appropriate for data that may not be independent and identically distributed, as is generally the case in aimless shooting) for a given model as a function of the amount of data collected, and by default uses the convergence of this parameter as a termination criterion.

This "information error" is a property of likelihood maximization derived from the first and second derivatives of the log likelihood function evaluated at the model optimization solution (that is, at the "maximum likelihood estimator"). It is a measure of the "information" about the optimization parameters that is stored in the dataset, as described by the sensitivity of the optimization result to changes in the values of the parameters. As the information error decreases, the confidence that the optimized model is the "best" possible one for the given dataset increases.

In the context of aimless shooting, the information error in the model can be interpreted as a metric of convergence of sampling *within the explored regions of phase space*. That is, it does not resolve :ref:`TheProblem` in that it remains impossible to determine (except through retrospection) whether the regions of phase space being sampled are sufficient to describe the system. It does, however, address issue (2) above by assessing convergence of the sampling *that does occur*; that is, if information error is converged, one can be assured that further sampling from the same distribution will not improve the result appreciably (to within the chosen tolerance). Furthermore, because information error is derived from likelihood maximization and therefore does not involve any simulation, it is comparatively much faster and more resource-efficient than committor analysis and can be assessed quite frequently without wasting significant resources, addressing issue (1).

To reiterate, information error should not *supplant* committor analysis but should *augment* it. Aimless shooting should be performed until information error has converged, and only then should the more expensive step of committor analysis be run. Besides increasing efficiency, the addition of information error to the aimless shooting workflow adds information to the committor analysis result when convergence has been achieved: if a positive (desirable) result is obtained, one can be assured that it is also a good approximation of the best possible committor analysis result for the given region(s) of sampling and candidate reaction coordinate parameters. Conversely, if a negative (undesirable) result is obtained, it must be due to omission of one or more key reaction coordinate parameters from the reaction coordinate, or else a fundamental flaw in the sampling, but *not* due to insufficient sampling of the sampled region(s).