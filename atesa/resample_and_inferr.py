"""
Helper file to call utilities.resample and then information_error.main in sequence in a single process.
"""

import os
import pickle
from atesa import utilities
from atesa import information_error
from atesa import rc_eval
from atesa.main import Thread

def main():
    # Read settings from pickle file
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:   # replace with more informative error message
        raise FileNotFoundError('the working directory must contain a valid settings.pkl file, which is generated '
                                'automatically when running ATESA, but one was not found in the working directory: '
                                + os.getcwd())

    # Call resample and then information_error
    utilities.resample(settings, partial=True)
    information_error.main()

    # If terminating based on information error, run rc_eval.py automatically
    if os.path.exists('info_err.out') and len(open('info_err.out', 'r').readlines()) > 0:
        # Get data from info_err.out
        last_line = open('info_err.out', 'r').readlines()[-1].split(' ')
        last_value = last_line[1]
        last_n_data = last_line[0]
        if float(last_value) <= settings.information_error_threshold:   # this is the criterion for termination
            # Extract RC from appropriate lmax output file
            rc = open(str(last_n_data) + '_lmax.out', 'r').readlines()[1].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '').replace('\n','').replace(' ','')

            # Call rc_eval
            rc_eval.main(settings.working_directory, rc, 'as_decorr_' + str(last_n_data) + '.out', False)

if __name__ == "__main__":
    main()
