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
    thread.current_type, thread.current_name = thread.get_next_step(settings)

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
    if os.path.exists(settings.path_to_templates):
        env = Environment(loader=FileSystemLoader(settings.path_to_templates))
    else:
        sys.exit('Error: could not locate templates folder: ' + settings.path_to_templates)

    batchfiles = []     # initialize list of batch files to fill out
    thread.traj_files = []  # to clear out previous traj_files if any exist
    for job_index in range(len(thread.current_type)):
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]
        template = env.get_template(thread.get_batch_template(type, settings))
        filled = template.render(name=thread.name + '_' + name,
                                 nodes=eval('settings.' + type + '_nodes'),
                                 taskspernode=eval('settings.' + type + '_ppn'),
                                 walltime=eval('settings.' + type + '_walltime'),
                                 mem=eval('settings.' + type + '_mem'),
                                 solver=eval('settings.' + type + '_solver'),
                                 inp=settings.path_to_input_files + '/' + settings.job_type + '_' + type + '_' + settings.md_engine + '.in',
                                 out=thread.name + '_' + name + '.out',
                                 prmtop=thread.topology,
                                 inpcrd=thread.coordinates,
                                 rst=thread.name + '_' + name + '.rst',
                                 nc=thread.name + '_' + name + '.nc',
                                 working_directory=settings.working_directory)
        newfilename = thread.name + '_' + name + '.' + settings.batch_system
        with open(newfilename, 'w') as newfile:
            newfile.write(filled)
            newfile.close()

        batchfiles.append(newfilename)
        thread.traj_files.append(thread.name + '_' + name + '.nc')   # todo: what's the best way to support general trajectory formats here?

    ### Submit batch files to task manager ###
    taskmanager = factory.taskmanager_factory(settings.task_manager)
    thread.jobids = []      # to clear out previous jobids if any exist
    for file in batchfiles:
        thread.jobids.append(taskmanager.submit_batch(None, file, settings))

    if thread not in running:
        running.append(thread)
    return running
