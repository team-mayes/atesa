"""
This portion of the program is responsible for handling setup of the appropriate batch script(s) for the next step in a
Thread, passing them to a task manager to submit them, and updating the list of currently running threads accordingly.
"""

import jinja2
import os
import sys
from jinja2 import Environment, FileSystemLoader
try:
    import factory
except ModuleNotFoundError:
    import atesa_v2.factory as factory

def process(thread, running, settings):
    """
    The main function of process.py. Reads the thread to identify the next step, then builds and submits the batch
    file(s) as needed.

    Parameters
    ----------
    thread : Thread()
        The Thread object on which to act
    running : list
        The list of currently running threads, which will be updated as needed
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    running : list
        The updated list of currently running threads after this process step

    """

    ### Determine next step and, if appropriate, build corresponding list of batch files ###
    thread.current_type = thread.get_next_step(settings)  # todo: test to make sure this line updates the thread
    if thread.current_type == 'terminate':
        thread.terminated = True

    if thread.terminated:
        if thread in running:
            running.remove(thread)
            if running is None:
                running = []       # to keep running as a list, even if empty
            return running
        else:
            return running

    # Set Jinja2 environment
    if os.path.exists(sys.path[0] + '/atesa_v2/data/templates'):
        env = Environment(loader=FileSystemLoader(settings.path_to_templates))
    else:
        sys.exit('Error: could not locate templates folder: ' + settings.path_to_templates)

    batchfiles = []     # initialize list of batch files to fill out
    thread.traj_files = []  # to clear out previous traj_files if any exist
    for job in thread.current_type: # todo: this won't work; I don't want different walltimes for fwd and bwd jobs, for example! I should rework my whole strategy here.
        template = env.get_template(thread.get_batch_template(job, settings))
        filled = template.render(name=thread.name + '_' + job,
                                 nodes=eval('settings.' + job + '_nodes'),
                                 taskspernode=eval('settings.' + job + '_ppn'),
                                 walltime=eval('settings.' + job + '_walltime'),
                                 solver=eval('settings.' + job + '_solver'),
                                 inp=settings.path_to_input_files + '/' + settings.job_type + '_' + job + '_' + settings.md_engine + '.in',
                                 out=thread.name + '_' + job + '.out',
                                 prmtop=thread.topology,
                                 inpcrd=thread.coordinates,
                                 rst=thread.name + '_' + job + '.rst',
                                 nc=thread.name + '_' + job + '.nc',
                                 mem=eval('settings.' + job + '_mem'),
                                 working_directory=settings.working_directory)
        newfilename = thread.name + '_' + job + '.' + settings.batch_system
        with open(newfilename, 'w') as newfile:
            newfile.write(filled)
            newfile.close()

        batchfiles.append(newfilename)
        thread.traj_files.append(thread.name + '_' + job + '.nc')   # todo: what's the best way to support general trajectory formats here?

    ### Submit batch files to task manager ###
    taskmanager = factory.taskmanager_factory(settings.task_manager)
    thread.jobids = []      # to clear out previous jobids if any exist
    for file in batchfiles:
        thread.jobids.append(taskmanager.submit_batch(None, file, settings))

    if thread not in running:
        running.append(thread)
    return running
