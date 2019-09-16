"""
This portion of the program is responsible for handling update of the results, checking global termination criteria, and
implementing the calls to JobType methods to control the value of the thread.coordinates attribute for the next step.
"""

try:
    import factory
except ModuleNotFoundError:
    import atesa_v2.factory as factory

def interpret(thread, allthreads, settings):
    """
    The main function of interpret.py. Makes calls to JobType methods to update results, check termination criteria, and
    update thread.coordinates

    Parameters
    ----------
    thread : Thread
        The Thread object on which to act
    allthreads : list
        The list of all extant Thread objects
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    termination: bool
        True if a global termination criterion has been met; False otherwise

    """

    jobtype = factory.jobtype_factory(settings.job_type)

    # todo: need to insert a check that the job step produced the required file(s) BEFORE proceeding into check_termination and update_results
    if jobtype.check_for_successful_step(thread, settings):
        termination = jobtype.check_termination(thread, allthreads, settings)   # check global termination criteria, if any
        jobtype.update_results(thread, allthreads, settings)                    # update results as needed
    jobtype.algorithm(thread, allthreads, settings)                         # set thread parameters for next step

    return termination
