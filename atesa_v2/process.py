"""
This portion of the program is responsible for handling setup of the appropriate batch script(s) for the next step in a
Thread, passing them to a task manager to submit them, and updating the list of currently running threads accordingly.
"""

import jinja2

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

    if thread.terminated:
        if thread in running:
            return running.remove(thread)
        else:
            return running

    # Set jinja2 environment
    if os.path.exists(sys.path[0] + '/atesa_v2/data/templates'):
        env = Environment(loader=FileSystemLoader(sys.path[0] + '/atesa_v2/data/templates'))
    else:
        sys.exit('Error: could not locate templates folder: ' + sys.path[0] + '/atesa_v2/data/templates')

    thread.current_type = thread.get_next_step()    # todo: test to make sure this line updates the thread; if not, need to return it?
    batchfiles = []     # initialize list of batch files to fill out
    for job in thread.current_type:
        template = thread.get_batch_template(job, settings)
        jobstring = ''
        batchfiles.append(filled)
