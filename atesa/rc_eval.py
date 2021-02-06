"""
rc_eval.py
Standalone script to evaluate RC values given an aimless shooting working directory and reaction coordinate definition
"""

import sys
import os
import glob
import argparse
import pickle
import time
import math
from atesa import utilities

def update_progress(progress, message='Progress', eta=0, quiet=False):
    """
    Print a dynamic progress bar to stdout.

    Credit to Brian Khuu from stackoverflow, https://stackoverflow.com/questions/3160699/python-progress-bar

    Parameters
    ----------
    progress : float
        A number between 0 and 1 indicating the fractional completeness of the bar. A value under 0 represents a 'halt'.
        A value at 1 or bigger represents 100%.
    message : str
        The string to precede the progress bar (so as to indicate what is progressing)
    eta : int
        Number of seconds to display as estimated completion time (converted into HH:MM:SS)
    quiet : bool
        If True, suppresses output entirely

    Returns
    -------
    None

    """

    if quiet:
        return None

    barLength = 10  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done!          \r\n"
    block = int(round(barLength * progress))
    if eta:
        # eta is in seconds; convert into HH:MM:SS
        eta_h = str(math.floor(eta/3600))
        eta_m = str(math.floor((eta % 3600) / 60))
        eta_s = str(math.floor((eta % 3600) % 60)) + ' '
        if len(eta_m) == 1:
            eta_m = '0' + eta_m
        if len(eta_s) == 2:
            eta_s = '0' + eta_s
        eta_str = eta_h + ':' + eta_m + ':' + eta_s
        text = "\r" + message + ": [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), round(progress * 100, 2), status) + " ETA: " + eta_str
    else:
        text = "\r" + message + ": [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), round(progress * 100, 2), status)
    sys.stdout.write(text)
    sys.stdout.flush()


def main(working_directory, rc_definition, as_out_file, extrema=False):
    """
    The main function of rc_eval.py. Accepts an aimless shooting working directory and a reaction coordinate definition,
    producing in that directory a new file named 'rc.out' (overwriting if one already exists) identifying each shooting
    point (files in the directory whose names end in "_init.rst7") and its corresponding reaction coordinate value in a
    sorted list.

    If extrema == True, skips producing rc.out and just returns the minimum and maximum RC value for a single accepted
    shooting move, which is useful when preparing umbrella sampling simulations.

    Parameters
    ----------
    working_directory : str
        The path to the aimless shooting working directory in which to act
    rc_definition : str
        A reaction coordinate definition formatted as a string of python-readable code with "CV[X]" standing in for the
        Xth CV value (one-indexed); this RC definition should be in terms of reduced variables (values between 0 and 1)
    as_out_file : str
        Path to the aimless shooting output file used to build the reaction coordinate. Usually this should be a
        decorrelated file (named with "decorr").
    extrema : bool
        If True, skips producing rc.out and just returns the minimum and maximum RC value for a single accepted shooting
        move, which is useful when preparing umbrella sampling simulations.

    Returns
    -------
    None

    """

    # Change to working directory
    os.chdir(working_directory)

    # Unpickle settings object for use in utilities.get_cvs
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:   # replace with more informative error message
        raise FileNotFoundError('the working directory must contain a valid settings.pkl file, which is generated '
                                'automatically when running ATESA, but one was not found in the working directory: '
                                + working_directory)

    if not settings.job_type == 'aimless_shooting':
        raise RuntimeError('rc_eval.py can only be called on an aimless shooting working directory, but the provided '
                           'directory (' + working_directory + ') contains a settings.pkl file with job_type = ' +
                           settings.job_type)

    settings.as_out_file = as_out_file      # for reducing CVs properly
    settings.include_qdot = False           # unnecessary for our purposes

    if extrema:
        from atesa.main import Thread
        print('Evaluating final RC values of forward and backward trajectories from an accepted shooting move...')
        result = []
        allthreads = pickle.load(open('restart.pkl', 'rb'))
        for thread in allthreads:
            if thread.history.last_accepted > -1:   # if accepted move exists in thread
                for job_index in range(2):
                    cvs = utilities.get_cvs(thread.history.prod_trajs[thread.history.last_accepted][job_index], settings, reduce=True).split(' ')
                    result.append(utilities.evaluate_rc(rc_definition, cvs))
                print(' Shooting move name: ' + thread.history.init_coords[thread.history.last_accepted][0])
                print(' extrema: ' + str(result))
                return None   # to exit the script after returning extrema
        raise RuntimeError('none of the shooting moves in the working directory appear to contain any accepted moves.')

    # Obtain list of shooting point coordinate files
    file_list = glob.glob('*_init.rst7')
    if not file_list:
        raise FileNotFoundError('no valid shooting point files (as given by names ending in \'_init.rst7\') were found '
                                'in the working directory \'' + working_directory + '\'. Is this an aimless shooting '
                                'working directory?')

    # Iterate through the list, calling evaluate_rc for each one and storing the result
    results = []
    count = 0
    count_to = len(file_list)
    update_progress(0, 'Evaluating RC values')
    speed_data = [0, 0]
    for file in file_list:
        t = time.time()
        cv_list = utilities.get_cvs(file, settings, reduce=True).split(' ')
        results.append([file + ': ', utilities.evaluate_rc(rc_definition, cv_list)])
        this_speed = time.time() - t
        speed_data = [(speed_data[1] * speed_data[0] + this_speed) / (speed_data[1] + 1), speed_data[1] + 1]
        count += 1
        eta = (count_to - count) * speed_data[0]
        update_progress(count/count_to, 'Evaluating RC values', eta=eta)
    results = sorted(results, key=lambda x: abs(float(x[1])))  # sort results by absolute value of RC

    # Create and write to rc.out file
    open('rc.out', 'w').close()
    with open('rc.out', 'a') as f:
        for result in results:
            f.write(result[0] + str(result[1]) + '\n')
        f.close()


if __name__ == '__main__':
    if len(sys.argv) == 5:
        extrema = sys.argv[4]
    else:
        extrema = False
    if len(sys.argv) <= 3:
        raise RuntimeError('not enough arguments; you gave ' + str(len(sys.argv) - 1) + ' but rc_eval.py takes 3 or 4')
    main(sys.argv[1], sys.argv[2], sys.argv[3], extrema)
