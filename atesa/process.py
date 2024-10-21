"""
This portion of the program is responsible for handling setup of the appropriate batch script(s) for the next step in a
Thread, passing them to a task manager to submit them, and updating the list of currently running threads accordingly.
"""

import django.template
import os
import sys
import warnings
from atesa import factory
from django.conf import settings as django_settings

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

    if not django_settings.configured:  # need to configure just once per thread
        django_settings.configure()

    ### Determine next step and, if appropriate, build corresponding list of batch files ###
    if not thread.skip_update:
        thread.current_type, thread.current_name = thread.get_next_step(settings)
    else:
        thread.skip_update = False

    if thread.terminated:
        if thread in running:
            running.remove(thread)
            if running is None:
                running = []       # to keep running as a list, even if empty
        return running

    # Set file extention for trajectory # todo: figure out a way to do this properly
    traj_format = '.nc'
    if settings.md_engine == 'cp2k':
        traj_format = '.dcd'

    batchfiles = []         # initialize list of batch files to fill out
    jobtype = factory.jobtype_factory(settings.job_type)    # get jobtype for calling jobtype.update_history
    this_inpcrd = jobtype.get_inpcrd(thread)
    for job_index in range(len(thread.current_type)):
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]

        these_kwargs = {'name': thread.name + '_' + name,
                         'nodes': str(eval('settings.' + type + '_nodes')),
                         'taskspernode': str(eval('settings.' + type + '_ppn')),
                         'walltime': str(eval('settings.' + type + '_walltime')),
                         'mem': str(eval('settings.' + type + '_mem')),
                         'solver': str(eval('settings.' + type + '_solver')),
                         'out': thread.name + '_' + name + '.out',
                         'prmtop': thread.topology,
                         'inpcrd': this_inpcrd[job_index],
                         'rst': thread.name + '_' + name + '.rst7',
                         'nc': thread.name + '_' + name + traj_format,
                         'working_directory': settings.working_directory,
                         'thread_name': thread.name,
                         'extra': str(eval('settings.' + type + '_extra'))}

        inp = jobtype.get_input_file(thread, job_index, settings, **these_kwargs)

        template = settings.env.get_template(thread.get_batch_template(type, settings))
        these_kwargs.update({'inp': inp})

        filled = template.render(django.template.Context(these_kwargs))
        newfilename = thread.name + '_' + name + '.' + settings.batch_system
        try:
            with open(newfilename, 'w') as newfile:
                newfile.write(filled)
                newfile.close()
        except OSError as e:
            if 'name too long' in str(e):
                warnings.warn('Encountered too-long filename: ' + newfilename + '. Skipping submission of this thread '
                              'to the task manager, terminating it, and continuing. Consider renaming a sample of the '
                              'initial coordinate files you like best and starting a new ATESA run. If this limitation '
                              'is a problem for you, please raise an issue on GitHub detailing your use-case.')
                thread.terminated = True
                return running

        batchfiles.append(newfilename)
        jobtype.update_history(thread, settings, **these_kwargs)

    ### Submit batch files to task manager ###
    taskmanager = factory.taskmanager_factory(settings.task_manager)
    thread.jobids = []      # to clear out previous jobids if any exist
    for file in batchfiles:
        thread.jobids.append(taskmanager.submit_batch(file, settings))

    if thread not in running:
        running.append(thread)
    return running
