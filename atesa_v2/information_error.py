"""
Script to evaluate and store information error in support of the information error convergence criterion in aimless
shooting. Implemented as a separate script to facilitate multiprocessing.
"""

import sys
import os
import time
import re
import pickle
import subprocess
import shutil
import glob
from atesa_v2.main import Thread
from atesa_v2 import utilities
from statsmodels.tsa.stattools import kpss

def main(as_raw):
    """
    Evaluate the information error during aimless shooting and output results to info_err.out.

    Reads as_raw.out from the present directory to evaluate the information error of the aimless shooting process as
    well as the likelihood that the series of information error values represents a trend stationary process
    (informally, is converged) using the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test.

    This function depends on lmax.py for evaluations of the information error.

    Parameters
    ----------
    as_raw : str
        Name of the raw aimless shooting file to operate on

    Returns
    -------
    None

    """

    # Initialize info_err.out if it does not yet exist
    if not os.path.exists('info_err.out'):
        open('info_err.out', 'w').close()

    # Initialize regular expressions for obtaining strings later
    pattern = re.compile('[0-9.]+')     # pattern to match information error number
    pattern2 = re.compile('[0-9]+')     # pattern to match amount of data in as_raw

    # Get number string from as_raw
    length = pattern2.findall(as_raw)[0]

    # Read settings from pickle file
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:   # replace with more informative error message
        raise FileNotFoundError('the working directory must contain a valid settings.pkl file, which is generated '
                                'automatically when running ATESA, but one was not found in the working directory: '
                                + os.getcwd())

    # Run resampling to make decorrelated data file
    if not settings.DEBUG and not settings.resample_override:
        utilities.resample(settings, suffix='_' + length, write_raw=False)

    # Set resample_override to False (or back to False; it will only be set to True when called by resample)
    settings.resample_override = False
    pickle.dump(settings, open('settings.pkl', 'wb'), protocol=2)

    # Exit if as_decorr.out has no content
    if len(open('as_decorr_' + length + '.out', 'r').readlines()) == 0:
        return None

    # Run likelihood maximization on current as_decorr.out
    if settings.include_qdot:
        q_str = 'present'
    else:
        q_str = 'absent'
    command = 'lmax.py -i as_decorr_' + length + '.out -q ' + q_str + ' --automagic -o ' + as_raw + '_lmax.out'
    subprocess.check_call(command.split(' '), stdout=sys.stdout, preexec_fn=os.setsid)     # check_call waits for completion
    if not os.path.exists(as_raw + '_lmax.out'):
        raise FileNotFoundError('Likelihood maximization did not produce output file:' + as_raw + '_lmax.out')

    # Assemble model dimensions just obtained for calling lmax for other data sets
    model_str = open(as_raw + '_lmax.out', 'r').readlines()[1]
    model_str = model_str[model_str.rindex(': ') + 2:]
    dims = ''
    while 'CV' in model_str:
        model_str = model_str[model_str.index('CV') + 2:]    # chop off everything up to and including the next 'CV'
        if ' ' in model_str:
            dims += model_str[:model_str.index(' ')] + ' '  # CV index is everything up until the next space
        else:
            dims += model_str[:model_str.index('\n')]       # last dimensions has a newline instead of a space
    if dims == '':
        raise RuntimeError('Likelihood maximization output file is improperly formatted: ' + as_raw + '_lmax.out')


    datalengths = [pattern2.findall(item)[-1] for item in glob.glob('as_decorr_*.out') if not 'lmax' in item and not length + '.out' in item]   # todo: sort
    if datalengths:
        open('info_err_temp.out', 'w').close()

    # Call lmax for each further dataset and write new info_err output file
    for datalength in datalengths:
        command = 'lmax.py -i as_decorr_' + datalength + '.out -q ' + q_str + ' -f ' + dims + ' -k ' + str(int(len(dims.split(' ')))) + ' -o ' + datalength + '_redo_lmax.out'
        subprocess.check_call(command.split(' '), stdout=sys.stdout, preexec_fn=os.setsid)
        inf_err = float(pattern.findall(open(datalength + '_redo_lmax.out', 'r').readlines()[-1])[0])
        open('info_err_temp.out', 'a').write(datalength + ' ' + str(inf_err) + '\n')

    # Add information error from previously completed lmax for this length
    inf_err = float(pattern.findall(open(as_raw + '_lmax.out', 'r').readlines()[-1])[0])
    open('info_err_temp.out', 'a').write(length + ' ' + str(inf_err) + '\n')

    # Copy info_err_temp to info_err.out
    shutil.move('info_err_temp.out', 'info_err.out')


if __name__ == "__main__":
    main(sys.argv[1])
