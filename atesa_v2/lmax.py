"""
Likelihood maximization script. This program is designed to be entirely separable from ATESA in that it can be called
manually to perform likelihood maximization to user specifications and with arbitrary input files; however, it is
required by ATESA's aimless shooting information error convergence criterion.
"""

import sys
import os
import numpy
import time
import math
import itertools
import argparse
import warnings
import pickle
from scipy import optimize
from scipy import stats
from scipy.special import erf

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


def objective_function(params, A_data, B_data):
    """
    Evaluate the negative log likelihood function for the given parameters and lists of observations.

    This function evaluates the goodness of fit of the given parameters and data to an error function ansatz, as
    described in Peters, 2012. Chem. Phys. Lett. 554: 248.

    Designed to be called by an optimization routine to obtain the best fitting params.

    Parameters
    ----------
    params : list
        Parameters for the current model to be tested
    A_data : list
        List of observations from aimless shooting that committed to basin "A" (usually the reactants)
    B_data : list
        List of observations from aimless shooting that committed to basin "B" (usually the products)

    Returns
    -------
    negative_log_likelihood : float
        The negative log likelihood of the fit to the error function for the given parameters and observations

    """

    def erflike(arg):
        pl = numpy.ones(len(arg))
        ml = numpy.negative(numpy.ones(len(arg)))
        return numpy.where(arg > 5.7, pl, numpy.where(arg < -5.7, ml, erf(arg)))

    qa = params[0] + numpy.inner(params[1:], A_data)
    qb = params[0] + numpy.inner(params[1:], B_data)
    sum = numpy.sum(numpy.log((1 - erflike(qa)) / 2)) + numpy.sum(numpy.log((1 + erflike(qb)) / 2))

    return -1 * sum


def two_line_test(results):
    """
    Perform a double linear regression on intersecting subsets of the data in results to determine whether to terminate
    and how many dimensions to return in the RC during automagic.

    Can only be called with len(results) >= 5.

    Parameters
    ----------
    results : list
        List of dictionary objects indexed by step of automagic, each possessing attribute 'fun' giving the optimization
        score for that step

    Returns
    -------
    out : int
        Index of selected 'best' RC from two-line test; or, -1 if no best RC could be determined

    """

    if len(results) < 5:
        raise RuntimeError('two_line_test can only be called with at least 5 optimized models')

    best_closest = []  # result for which the intersection is closest to the shared point
    for test_index in range(len(results) - 2):  # - 2 to account for minimum of two points in each line
        first_segment = range(1, 3 + test_index)
        second_segment = range(first_segment[-1], len(results) + 1)
        opt1 = stats.linregress(first_segment, [results[i - 1].fun for i in first_segment])
        opt2 = stats.linregress(second_segment, [results[i - 1].fun for i in second_segment])

        # Now evaluate closest point in results to the intersection of the two lines
        x_intersect = (opt1.intercept - opt2.intercept) / (opt2.slope - opt1.slope)
        y_intersect = (opt1.slope * x_intersect) + opt1.intercept
        x_val = 0       # initialize index for keeping track of x values
        min_diff = -1   # initialize smallest distance between intersection and point
        closest = 0     # initialize index of closest point to intersection
        for result in results:
            y_val = result.fun
            x_val += 1
            y_diff = y_val - y_intersect
            x_diff = x_val - x_intersect
            diff = numpy.sqrt(y_diff**2 + x_diff**2)
            if min_diff < 0:
                min_diff = diff
                closest = [x_val, diff]
            elif diff < min_diff:
                min_diff = diff
                closest = [x_val, diff]

        # if the closest point to the intersection is the shared point of the lines;
        if closest[0] == test_index + 2:
            if not best_closest:                    # for the first time
                best_closest = [closest, opt1, opt2]
            elif closest[1] < best_closest[0][1]:   # update the closest yet
                best_closest = [closest, opt1, opt2]

    if not best_closest:    # no pairs of lines whose intersection was closest to their shared point
        return -1

    slope_fract = best_closest[2].slope / best_closest[1].slope
    if slope_fract > 0.35:  # best point does not meet threshold for relative difference in slopes
        return -1
    else:                   # DOES meet threshold; return the index of the passing result
        return best_closest[0][0] - 1   # - 2 because of different indexing standards


def main(**kwargs):
    """
    Main runtime function of lmax.py.

    Assembles lists of models to optimize in the form of lists of CVs, passes them to optimize, interprets results, and
    repeats or terminates in accordance with argument-dependent termination criteria.

    Parameters
    ----------
    kwargs : dict
        Dictionary object containing arguments

    Returns
    -------
    None

    """

    # Ensure existence and validity of input file
    input_file = kwargs['i'][0]
    if not os.path.exists(input_file):
        raise FileNotFoundError('could not find input file: ' + input_file)
    input_file_lines = open(input_file, 'r').readlines()
    open(input_file, 'r').close()
    if False in [char == 'A' or char == 'B' for char in [line[0] for line in input_file_lines]]:
        raise RuntimeError('input file ' + input_file + ' does not have \'A\' or \'B\' as the first character in each '
                           'line. Is this the correct file? Be sure to remove any blank lines.')

    # Bring in other arguments, just for neatness
    dims = kwargs['k'][0]
    fixed = kwargs['f']      # we actually want this one to stay a list
    qdot = kwargs['q'][0]
    running = kwargs['r'][0]
    output_file = kwargs['o'][0]
    automagic = kwargs['automagic']
    quiet = kwargs['quiet']

    # Ignore arguments as described in documentation
    if running:
        fixed = []
        dims = running
    if automagic:
        fixed = []
        dims = -1
        running = 0

    # Load settings object from .pkl file if present, to check for information criterion override
    try:
        settings = pickle.load(open('settings.pkl', 'rb'))
    except FileNotFoundError:
        pass

    # Get data from input file, and determine minimum and maximum values for each CV, reduce data
    input_data = [[float(item) for item in
                   line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')]
                  for line in input_file_lines]     # [[obs1cv1, obs1cv2], [obs2cv1, obs2cv2]]
    A_data = [[float(item) for item in line.replace('A <- ', '').replace(' \n', '').replace('\n', '').split(' ')] for
              line in input_file_lines if line[0] == 'A']
    B_data = [[float(item) for item in line.replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')] for
              line in input_file_lines if line[0] == 'B']
    mapped = list(map(list, zip(*input_data)))      # [[obs1cv1, obs2cv1], [obs1cv2, obs2cv2]]
    minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]    # [[mincv1, mincv2], [maxcv1, maxcv2]]
    N = len(input_file_lines)   # number of observations
    NA = len(A_data)            # number of observations that committed to A...
    NB = len(B_data)            # ... and to B
    num_cvs = len(minmax[0])    # number of CVs recorded in each observation
    reduced_A = [[(A_data[jj][ii] - minmax[0][ii]) / (minmax[1][ii] - minmax[0][ii]) for ii in range(num_cvs)] for jj in range(NA)]
    reduced_B = [[(B_data[jj][ii] - minmax[0][ii]) / (minmax[1][ii] - minmax[0][ii]) for ii in range(num_cvs)] for jj in range(NB)]

    if qdot == 'present' or qdot == 'ignore':
        if not num_cvs % 2 == 0:
            raise RuntimeError('likelihood maximization was attempted with input file: ' + input_file + ' and '
                               'include_qdot (q) = True, but this input file has an odd number of entries per line. Are'
                               ' you sure it includes rate-of-change data?')
        num_cvs = int(num_cvs / 2)

    # Prepare for and then enter optimization loop
    termination = False     # initialize primary termination criterion flag
    termination_2 = False   # additional termination flag for use with qdot = 'present', to perform final optimization
    two_line_result = -1    # initialize current model dimensionality for automagic
    cv_combs = [[]]         # initialize list of CV combinations to iterate through
    results = []  # initialize for automagic
    while not termination and len(cv_combs[0]) <= N:
        # Initialize current best result
        current_best = [argparse.Namespace(), [0]]
        current_best[0].fun = math.inf

        # Assemble list of RCs to optimize
        if len(fixed) == dims:
            cv_combs = [fixed]
        elif running or automagic:
            cv_combs = [fixed + [new] for new in range(1, num_cvs + 1) if not new in fixed]
        else:
            cv_combs = [comb for comb in itertools.combinations(range(1, num_cvs + 1), dims) if (not fixed or set(fixed).issubset(comb))]
        if qdot == 'present' and not termination_2:
            cv_combs_temp = cv_combs
            cv_combs = []
            for comb in cv_combs_temp:
                cv_combs.append([])
                for item in comb:
                    cv_combs[-1].append(item)
                    cv_combs[-1].append(item + num_cvs)

        # Perform optimization
        start_params = [0 for null in range(len(cv_combs[0]) + 1)]  # + 1 for constant term
        count = 0
        count_to = len(cv_combs)
        update_progress(0, 'Optimizing ' + str(count_to) + ' combination(s) of CVs', quiet=quiet)
        speed_data = [0,0]
        for comb in cv_combs:
            t = time.time()
            this_A = []
            this_B = []
            for index in comb:  # produce k-by-len(A_data) matrices (list of lists) for the selected CVs
                this_A.append([obs[index - 1] for obs in reduced_A])
                this_B.append([obs[index - 1] for obs in reduced_B])
            this_A = list(map(list, zip(*this_A)))  # transpose the matrices to get desired format
            this_B = list(map(list, zip(*this_B)))
            this_result = optimize.minimize(objective_function, numpy.asarray(start_params), (this_A, this_B),
                                            method='BFGS', options={"disp": False, "maxiter": 20000 * (len(comb) + 1)}) # try SR1?
            if this_result.fun < current_best[0].fun:
                current_best = [this_result, comb]
            this_speed = time.time() - t
            speed_data = [(speed_data[1] * speed_data[0] + this_speed) / (speed_data[1] + 1), speed_data[1] + 1]
            count += 1
            eta = (count_to - count) * speed_data[0]
            update_progress(count / count_to, 'Optimizing ' + str(count_to) + ' combination(s) of CVs', eta, quiet=quiet)

        # Update fixed and results parameters as needed
        if automagic:
            results.append(current_best)
        if running or automagic:
            fixed = current_best[1]
            if qdot == 'present':
                for item in fixed:
                    if item > num_cvs:
                        fixed.remove(item)

        # Check termination criteria
        if not running and not automagic:
            termination = True
        elif running and not automagic:
            if int(len(current_best[1])) == running:
                termination = True
        elif automagic and not termination_2:
            if len(results) >= 5:   # can only confidently check for convergence with at least 5 points
                two_line_result = two_line_test([result[0] for result in results])
                if two_line_result >= 0:
                    termination = True
                    current_best = results[two_line_result]
        if termination_2:
            termination = True
        if qdot == 'present' and termination and not termination_2:
            termination = False
            termination_2 = True
            fixed = current_best[1]
            for item in fixed:
                if item > num_cvs:
                    fixed.remove(item)
            dims = len(fixed)

    if automagic and two_line_result < 0:   # ran out of CVs to append and two_line_test never passed
        err = RuntimeError('The automagic convergence criterion was never satisfied even after including every '
                           'candidate CV in the model reaction coordinate.\nThis almost certainly indicates that either'
                           ' one or more key CVs are absent from the aimless shooting output file supplied, or that not'
                           ' enough unimportant CVs were included to give context to the important ones. Either way you'
                           ' should add more CVs to the list.\nThis error can by bypassed by running lmax.py in a '
                           'directory containing a settings.pkl file with the line "information_error_override = True" '
                           '(without quotes). If you did supply this setting, then you are seeing this message because '
                           'the settings.pkl file could not be found.')
        try:
            if settings.information_error_override:
                pass
            else:
                raise err
        except NameError:
            raise err

    # Return output in desired format
    rc_string = str('%.3f' % current_best[0].x[0]) + ' + ' + ' + '.join(['%.3f' % current_best[0].x[i+1] + '*CV' +
                        str(current_best[1][i]) for i in range(len(current_best[1]))])
    output_string = 'Likelihood maximization complete!\n' \
                    'The optimized reaction coordinate (with CVs indexed from 1) is: ' + rc_string + '\n' \
                    'The negative log likelihood of this model is: ' + '%.3f' % current_best[0].fun + '\n' \
                    'The mean information error for this model is: ' + '%.3f' % numpy.mean([float(numpy.sqrt(item)) for item in numpy.diag(current_best[0].hess_inv)])

    if output_file:
        open(output_file, 'w').write(output_string)
    else:
        print(output_string)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform LMAX on the given input data')
    parser.add_argument('-i', metavar='input_file', type=str, nargs=1, default=['as_decorr.out'],
                        help='input filename (output from aimless shooting). Default=as_decorr.out')
    parser.add_argument('-k', metavar='dimensionality', type=int, nargs=1, default=[0],
                        help='number of CVs to include in RC. Default=none')
    parser.add_argument('-f', metavar='fixed', type=int, nargs='*', default=[None],
                        help='CVs to require inside the RC. Default=none')
    parser.add_argument('-q', metavar='include_qdot', type=str, nargs=1, default=['absent'],
                        help='valid options are: "present", "absent", and "ignore" (quotes excluded). If "present" or '
                             '"ignore", the input file is assumed to include rate-of-change ("q") data for each CV '
                             '(formatted as in "A <- CV0 CV1 q0 q1"); in the former case, a q term will be included in '
                             'the RC for each CV, implementing inertial likelihood maximization. In the latter, no such'
                             ' terms will be included. Finally, if "absent", the q data will be assumed not to be '
                             'present in the input file. Default=absent')
    parser.add_argument('-r', metavar='running', type=int, nargs=1, default=[0],
                        help='if > 0, runs from k = 1 to "running" using the previously obtained k - 1 results as the '
                             'argument for f, ignoring the arguments passed for k and f. Default=0')
    parser.add_argument('-o', metavar='output_file', type=str, nargs=1, default=[''],
                        help='Prints output to a new file whose name is given with this argument, instead of directly '
                             'to the terminal. The file will be overwritten if it exists. Default=none')
    parser.add_argument('--quiet', action='store_true',
                        help='If this option is given, output to the terminal is suppressed and only the output file is'
                             ' written.')
    parser.add_argument('--automagic', action='store_true',
                        help='If this option is given, arguments passed for k, f, and r are ignored, and the best RC is'
                             ' determined based on the two-line method (see documentation).')

    arguments = vars(parser.parse_args())  # Retrieves arguments as a dictionary object

    # Suppress numpy.log warnings that occur frequently during normal operation
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in log')
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in double_scalars')
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in double_scalars')

    main(**arguments)
