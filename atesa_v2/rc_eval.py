"""
rc_eval.py
Standalone script to evaluate RC values given an aimless shooting working directory and reaction coordinate definition
"""

import sys
import os
import glob
import argparse
import pickle
from atesa_v2 import utilities

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
    for file in file_list:
        cv_list = utilities.get_cvs(file, settings, reduce=settings.rc_reduced_cvs).split(' ')
        results.append([file + ': ', utilities.evaluate_rc(rc_definition, cv_list)])
    results = sorted(results, key=lambda x: abs(float(x[1])))  # sort results by absolute value of RC

    # Create and write to rc.out file
    open('rc.out', 'w').close()
    with open('rc.out', 'a') as f:
        for result in results:
            f.write(result[0] + str(result[1]) + '\n')
        f.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
