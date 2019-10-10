"""
Script to evaluate and store information error in support of the information error convergence criterion in aimless
shooting. Implemented as a separate script to facilitate multiprocessing.
"""

import sys
import os

def main(as_raw):
    """
    Evaluate the information error during aimless shooting and output results to info_err.out.

    Reads as_raw.out from the present directory to evaluate the information error of the aimless shooting process as
    well as the likelihood that the series of information error values represents a trend stationary process
    (informally, is converged) using the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test.

    This function depends on atesa_lmax.py for evaluations of the information error.

    Parameters
    ----------
    as_raw : str
        Name of the raw aimless shooting file to operate on

    Returns
    -------
    None

    """

    # Read settings from pickle file
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:   # replace with more informative error message
        raise FileNotFoundError('the working directory must contain a valid settings.pkl file, which is generated '
                                'automatically when running ATESA, but one was not found in the working directory: '
                                + os.getcwd())

    # Run resampling to get as_decorr.out
    utilities.resample(settings, write_raw=False)

    if len(open('as_decorr.out', 'r').readlines()) > 0:  # proceed only if as_decorr.out has content
        # Run likelihood maximization on current as_decorr.out
        command = sys.executable + ' atesa_lmax.py -i as_decorr.out -q ' + str(settings.include_qdot) + \
                  ' --running ' + str(int(ndims)) + ' --output_file ' + str(len_data) + '_dimchk.out'
        process = subprocess.Popen(['bash', str(len_data) + '_dimchk.sh'], stdout=subprocess.PIPE,
                                   preexec_fn=os.setsid)


if __name__ == "__main__":
    main(sys.argv[1])
