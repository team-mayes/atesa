.. _OnTerminationCriteria:

On Termination Criteria
=======================

Along with ATESA, I have introduced a novel termination criterion for aimless shooting that I believe should be used in all cases. In summary, the idea is to estimate the error in the reaction coordinate produced with likelihood maximization using the diagonal terms of the Godambe Information matrix. Use of this termination criterion prevents wasted simulations and adds confidence to the final result. For more details and context, see the rest of this page.

.. _TheSamplingProblem:

The Sampling Problem
--------------------

Most any molecular system that is interesting enough to bother simulating operates within an extremely multidimensional state space. In the most general case, a complete exploration of this state space requires not only that the full range of each dimension be explored individually, but also the full range of combinations of dimensions; that is, the complexity of state spaces scales exponentially with the number of dimensions. Thus, perfect sampling of even small sections of state space is most often totally impractical.

If we know that our simulation will never fully explore state space, how do we know whether we've sampled enough to trust the results? This issue has been aptly dubbed "the sampling problem".

The only "solution" to the sampling problem is to circumvent it by relying on free energy to "weight" regions of state space by their relative importance in determining relevant properties of the system. That is, regions of state space with relatively lower free energy play a proportionally larger role in the system's behavior, so distributing sampling according to free energy is the most efficient means of arriving at a good approximation of the total behavior of the system. Of course, this is not really a solution at all -- that is, it does not address the central issue of identifying when our sampling is sufficient, or even whether we're sampling all of the important regions of state space in the first place. This can only be done *post hoc*, by verifying the putative results of simulations in some other way (*e.g.*, experimentally), and then retrospectively affirming that the sampling had been sufficient.

Relevance to Aimless Shooting
-----------------------------

The process of aimless shooting can be thought of as constrained sampling -- specifically, in the region of state space near enough to the reaction separatrix that both the product and reactant states are accessible. The aimless shooting algorithm distributes its sampling according to the acceptance probability (that is, the likelihood of forming a complete reactive trajectory from a given initial set of coordinates shot "forward" and "backward") and the free energy (simulations are more likely to proceed towards regions of lower energy, which in aggregate causes successive shooting moves to become centered at local energy minima along the separatrix). This sort of sampling does favor the most important transition states; however, this is no remedy to the sampling problem, as there remains no means of determining when sufficient sampling has been performed. But we can't go on sampling forever, so how do we choose when to stop?

Choosing When to Stop with Committor Analysis
---------------------------------------------

To date, the most common method of determining when to stop aimless shooting has been committor analysis (:ref:`CommittorAnalysis`). This requires that the aimless shooting data be mined for a reaction coordinate that describes the reaction progress (usually via likelihood maximization (:ref:`LikelihoodMaximizationTheory`)), and then that that reaction coordinate be verified by showing that unbiased simulations beginning from its center are (approximately) equally as likely to proceed to the reactants as to the products. This method is invaluable and should certainly be used to affirm the correctness of the final reaction coordinate, but as a termination criterion it suffers from two major issues:

#. It requires an extensive amount of additional simulation, which is expensive and time consuming; and

#. It provides no feedback about whether further sampling should be expected to improve the result.

In practice, this means that the researcher performing aimless shooting must either perform rounds of committor analysis many times throughout the study, or else collect such an excess of data that any unsatisfactory result in committor analysis is assumed to be due to issues elsewhere. Either option wastes time, both computational and real.

.. _InformationError:

Information Error
-----------------

ATESA introduces a new method of addressing the two shortcomings identified above, based on an assessment of the error in the likelihood maximization procedure. Specifically, ATESA measures the mean value of the parametric standard errors from the Godambe information matrix for a given model as a function of the amount of data collected, and by default uses a threshold on this parameter as its termination criterion during aimless shooting.

This "information error" is a property of likelihood maximization derived from the first and second derivatives of the log likelihood function evaluated at the model optimization solution (the "maximum likelihood estimator"). It is a measure of the "information" about the optimization parameters that is stored in the dataset, as described by the sensitivity of the optimization result to changes in the values of the parameters. As the information error decreases, the confidence that the optimized model is the "best" possible one for the given sets of data/observations and included dimensions increases.

In the context of aimless shooting, the information error in the model can be interpreted as a metric of convergence of sampling *within the explored regions of state space*. That is, it does not resolve :ref:`TheSamplingProblem` in that it remains impossible to determine (except through retrospection) whether the regions of state space being sampled are sufficient to describe the system. It does, however, address issue (2) above by assessing convergence of the sampling *that does occur*. In other words, to the extent that information error is converged, one can be assured that further sampling from the same distribution will not improve the result appreciably (to within the chosen tolerance). Furthermore, because information error is derived from likelihood maximization and therefore does not involve any simulation, it is comparatively much faster and more resource-efficient than committor analysis and can be assessed quite frequently without wasting significant resources, addressing issue (1).

To reiterate, information error should not *supplant* committor analysis but should *augment* it. Aimless shooting should be performed until information error has converged to within some reasonable tolerance, and only then should the more expensive step of committor analysis be attempted. Besides increasing efficiency, the addition of information error to the aimless shooting workflow adds information to the committor analysis result: if a positive (desirable) result is obtained with low information error, it is probably also a good approximation of the best possible committor analysis result for the given region(s) of sampling and candidate reaction coordinate parameters. Conversely, if a negative (undesirable) result is obtained, it is likely due to omission of one or more key parameters from the reaction coordinate, or else a fundamental flaw in the sampling, but *not* due to insufficient sampling of the sampled region(s).
