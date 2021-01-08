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
import warnings
from atesa import utilities
# from statsmodels.tsa.stattools import kpss    # deprecated

def main():
    """
    Evaluate the information error during aimless shooting and output results to info_err.out.

    Reads decorrelated aimless shooting output files from the present directory to evaluate the mean parametric variance
    based on the Godambe information error of the aimless shooting process.

    This function depends on lmax.py for evaluations of the information error.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # Initialize info_err.out if it does not yet exist
    if not os.path.exists('info_err.out'):
        open('info_err.out', 'w').close()

    # Initialize regular expressions for obtaining strings later
    pattern = re.compile('[0-9.]+')     # pattern to match information error number in lmax output file
    pattern2 = re.compile('[0-9]+')     # pattern to match amount of data in decorrelated datafile names

    # Assemble list of data lengths to perform information error evaluation for
    datalengths = [int(pattern2.findall(item)[-1]) for item in glob.glob('as_decorr_*.out') if not 'lmax' in item]
    if not datalengths:
        raise RuntimeError('attempted to evaluate information error, but found no files matching the pattern '
                           '\'as_decorr_*.out\' in the directory: ' + os.getcwd() + '\n'
                           'Was utilities.resample called first to produce these files?')
    datalengths = sorted(datalengths, key=lambda x: int(x))     # sort ascending to order chronologically
    if datalengths:
        open('info_err_temp.out', 'w').close()

    # Get length of largest datafile to perform two_line_test model optimization on
    length = max(datalengths)
    datalengths.remove(length)  # so as to skip repeat optimization later

    # Read settings from pickle file
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:   # replace with more informative error message
        raise FileNotFoundError('the working directory must contain a valid settings.pkl file, which is generated '
                                'automatically when running ATESA, but one was not found in the working directory: '
                                + os.getcwd())

    # Exit if the datafile has no content (as in, no decorrelated steps (almost impossible, included for robustness))
    if len(open('as_decorr_' + str(length) + '.out', 'r').readlines()) == 0:
        warnings.warn('skipping information error evaluation because decorrelated datafile as_decorr_' + str(length) +
                      '.out is empty. This may indicate an unusual error, especially if you see this message multiple '
                      'times in the course of sampling.', RuntimeWarning)
        return None

    # Run likelihood maximization on largest decorrelated datafile
    if os.path.exists(str(length) + '_lmax.out'):   # remove pre-existing output file if it exists
        os.remove(str(length) + '_lmax.out')
    if settings.include_qdot:
        q_str = 'present'
    else:
        q_str = 'absent'
    command = 'lmax.py -i as_decorr_' + str(length) + '.out -q ' + q_str + ' ' + settings.information_error_lmax_string + ' -o ' + str(length) + '_lmax.out --quiet'
    subprocess.check_call(command.split(' '))     # check_call waits for completion
    if not os.path.exists(str(length) + '_lmax.out'):
        raise FileNotFoundError('Likelihood maximization did not produce output file:' + str(length) + '_lmax.out')

    # Assemble model dimensions just obtained for calling lmax for other data sets
    model_str = open(str(length) + '_lmax.out', 'r').readlines()[1]
    model_str = model_str[model_str.rindex(': ') + 2:]
    dims = ''
    while 'CV' in model_str:
        model_str = model_str[model_str.index('CV') + 2:]    # chop off everything up to and including the next 'CV'
        if ' ' in model_str:
            dims += model_str[:model_str.index(' ')] + ' '  # CV index is everything up until the next space
        else:
            dims += model_str[:model_str.index('\n')]       # last dimensions has a newline instead of a space
    if dims == '':
        raise RuntimeError('Likelihood maximization output file is improperly formatted: ' + str(length) + '_lmax.out')

    # Call lmax for each further dataset and write new info_err output file
    for datalength in datalengths:
        command = 'lmax.py -i as_decorr_' + str(datalength) + '.out -q ' + q_str + ' -f ' + dims + ' -k ' + str(int(len(dims.split(' ')))) + ' -o ' + str(datalength) + '_lmax.out --quiet'
        subprocess.check_call(command.split(' '))
        inf_err = float(pattern.findall(open(str(datalength) + '_lmax.out', 'r').readlines()[-1])[0])
        open('info_err_temp.out', 'a').write(str(datalength) + ' ' + str(inf_err) + '\n')
        open('info_err_temp.out', 'a').close()

    # Add information error from previously completed lmax for this length
    inf_err = float(pattern.findall(open(str(length) + '_lmax.out', 'r').readlines()[-1])[0])
    open('info_err_temp.out', 'a').write(str(length) + ' ' + str(inf_err) + '\n')
    open('info_err_temp.out', 'a').close()

    # Move info_err_temp.out to info_err.out
    shutil.move('info_err_temp.out', 'info_err.out')

    # Deprecated (no longer using KPSS)
    # # Next, evaluate the KPSS statistic for the dataset in info_err.out
    # # Begin by suppressing warnings that occur frequently during normal operation of kpss
    # warnings.filterwarnings('ignore', message='p-value is greater than the indicated p-value')
    # warnings.filterwarnings('ignore', message='invalid value encountered in double_scalars')
    # warnings.filterwarnings('ignore', message='divide by zero encountered in double_scalars')
    #
    # # Get data back out from info_err.out (weird but robust)
    # info_err_lines = open('info_err.out', 'r').readlines()
    # info_errs = [float(line.split(' ')[1]) for line in info_err_lines]  # cast to float removes '\n'
    # kpssresults = []
    # for cut in range(len(info_errs) - 1):
    #     try:
    #         kpssresult = (kpss(info_errs[0:cut + 1], lags='auto')[1] - 0.01) / (0.1 - 0.01)
    #     except (ValueError, OverflowError):
    #         kpssresult = 1  # 100% certainty of non-convergence!
    #     kpssresults.append(kpssresult)
    #
    # # Write new info_err.out with kpss values
    # open('info_err.out', 'w').write(info_err_lines[0])  # no KPSS statistic for the first line
    # for line_index in range(1, len(info_err_lines)):
    #     line_data = [float(item) for item in info_err_lines[line_index].split(' ')]
    #     new_line = '%.0f' % line_data[0] + ' ' + '%.3f' % line_data[1] + ' ' + '%.3f' % kpssresults[line_index - 1] + '\n'
    #     open('info_err.out', 'a').write(new_line)
    # open('info_err.out', 'a').close()


if __name__ == "__main__":
    main()
