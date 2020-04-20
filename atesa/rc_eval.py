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


def main(working_directory, rc_definition):
    """
    The main function of rc_eval.py. Accepts an aimless shooting working directory and a reaction coordinate definition,
    producing in that directory a new file named 'rc.out' (overwriting if one already exists) identifying each shooting
    point (files in the directory whose names end in "_init.rst7") and its corresponding reaction coordinate value in a
    sorted list.

    Parameters
    ----------
    working_directory : str
        The path to the aimless shooting working directory in which to act
    rc_definition : str
        A reaction coordinate definition formatted as a string of python-readable code with "CV[X]" standing in for the
        Xth CV value (one-indexed); this RC definition should be in terms of reduced variables (values between 0 and 1)

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
        cv_list = utilities.get_cvs(file, settings, reduce=settings.rc_reduced_cvs).split(' ')
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
    main(sys.argv[1], sys.argv[2])
