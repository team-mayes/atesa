"""
Standalone script for converting equilibrium path sampling output files into free energy profiles via Boltzmann
weighting.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import os
import copy
from scipy import optimize

def objective_function(slope, probs, RC_values, kT):
    # Determines and returns a score measuring goodness of fit between the expected distribution of data from the given
    # slope, and the actual data in probs located at the RC_values
    energies = [slope * rc for rc in RC_values]
    putative_probs_temp = [float(np.exp(-1 * item / kT)) for item in energies]   # unnormalized
    putative_probs = [item/sum(putative_probs_temp) for item in putative_probs_temp]
    score = sum([abs(putative_probs[i] - probs[i]) for i in range(len(probs))])
    return score

def main(**kwargs):
    """
    Process equilibrium path sampling output file kwargs['i'] to obtain free energy profile.

    This is the main function of boltzmann_weight.py, which converts each EPS window into a stretch of the energy
    profile and stitches these stretches together into a continuous plot. Simply put, this function descretizes and
    reweights the data according to the Boltzmann weight (E = kT*ln(p), where p is the probability of a given state).

    Parameters
    ----------
    kwargs : dict
        Dictionary object containing arguments passed in from command line in "if __name__ == '__main__'" section.

    Returns
    -------
    None

    """

    # Change format of items in kwargs from one-length lists of inputs to just the inputs
    for item in kwargs.keys():
        if isinstance(kwargs[item], list):
            kwargs[item] = kwargs[item][0]

    # First, check that each input argument is valid
    if not os.path.exists(kwargs['i']):
        raise FileNotFoundError('could not find input file: ' + kwargs['i'])
    if not kwargs['t'] >= 0:
        raise RuntimeError('temperature argument \'t\' must be greater than or equal to zero')
    if not kwargs['n'] >= 2:
        raise RuntimeError('bins per window argument \'n\' must be greater than or equal to two')
    if not kwargs['c'] >= 0:
        raise RuntimeError('bootstrap cycles argument \'c\' must be greater than or equal to zero')
    if kwargs['b'] == 0:
        raise RuntimeError('bootstrap samples per window argument \'b\' cannot be equal to zero')
    if not 0 <= kwargs['e'] <= 1:
        raise RuntimeError('exclusion fraction argument \'e\' must be between zero and one')

    # Set Bolztmann factor for this temperature
    kT = kwargs['t'] * 0.001987  # kcal/mol-K

    # Load in data from input file
    file = open(kwargs['i'], 'r').readlines()
    open(kwargs['i'], 'r').close()

    windows = []    # nested list of format [[lower0, upper0], [lower1, upper1], ...]
    data = []       # nested list with indices corresponding windows, format [[x00, x01, ...], [x10, x11,...], ...]
    # alldata = []    # simple list [x00, x01, ... x0N, x10, x11, ...]

    # Determine the window boundaries
    for line in file:
        line = line.strip('\n')
        split = line.split()
        if [float('%.3f' % float(split[0])), float('%.3f' % (float(split[1])))] not in windows:
            windows.append([float('%.3f' % (float(split[0]))), float('%.3f' % (float(split[1])))])
            data.append([])

    windows.sort(key=lambda x: x[0])  # need to be sorted for building the PMF

    # Build the data nested list by sorting each line from the input file into the appropriate sublist of data
    for line in file:
        line = line.strip('\n')
        split = line.split(' ')
        if float('%.3f' % float(split[0])) <= float(split[2]) <= float('%.3f' % float(split[1])):
            data[windows.index([float('%.3f' % float(split[0])), float('%.3f' % float(split[1]))])].append(float(split[2]))
            # alldata.append(float(split[2]))
        else:
            raise RuntimeError('Impossible line in input file ' + kwargs['i'] + ': sampled value is not within window '
                               'boundaries.\n The offending line: ' + line)

    # Implement exclusion of the first 'e' of the data in each window
    for window_index in range(len(windows)):
        data[window_index] = data[window_index][int(len(data[window_index]) * kwargs['e']):]

    std_err = []                    # initialize standard error of the mean within each window
    subsampled_fullPMFs = []        # initialize list of PMFs from each bootstrapping cycle
    error_index = 0                 # initialize index to keep track of index within std_err when plotting

    if not kwargs['noplot']:
        fig0 = plt.figure()             # initialize first matplotlib figure
        ax0 = fig0.add_subplot(111)     # add axes to the first figure
        plt.ylabel('Free Energy (kcal/mol)', weight='bold')
        plt.xlabel('Reaction Coordinate', weight='bold')
        fig1 = plt.figure()             # initialize second matplotlib figure
        ax1 = fig1.add_subplot(111)     # add axes to the second figure
        plt.ylabel('Relative Frequency', weight='bold')
        plt.xlabel('Reaction Coordinate', weight='bold')

    cycle_index =  0    # initialize index of cycle for progress bar

    # 'c' cycles, plus one for final complete run
    for cycle in range(kwargs['c'] + 1):

        # Initialize list of last and second-to-last energy and RC values for stitching together adjacent windows
        boundary_values = [0, 0, 0, 0]

        # Also initialize list of full PMF energy values and RC values
        fullPMF = []
        fullRCs = []

        # Subsample data, if appropriate
        if cycle < kwargs['c']:
            temp_data = copy.deepcopy(data)  # necessary to avoid overwriting data in next line
            this_data = temp_data    # initialize
            for i in range(len(windows)):
                if kwargs['b'] > 0:
                    bootstrapN = kwargs['b']
                else:
                    bootstrapN = len(data[i])
                this_data[i] = np.random.choice(temp_data[i], bootstrapN)
        else:
            this_data = data    # full data set for final cycle

            # Calculate std_err from subsampled_fullPMFs
            std_err = [np.std(list) for list in np.transpose(subsampled_fullPMFs)]
            if std_err == []:  # no data to bootstrap, so we'll use "zero" for every error value
                std_err = [0 for null in range(kwargs['n'] * len(windows))]

        # Main processing step; for each window, evaluate the shape of the energy profile and adjust it vertically based on
        # the previous window to obtain the full free energy profile
        for window_index in range(len(windows)):
            cycle_index += 1
            update_progress(cycle_index / ((kwargs['c'] + 1) * len(windows)), 'Evaluating PMF')

            bin_edges = np.linspace(windows[window_index][0], windows[window_index][1], kwargs['n'] + 1)  # list RC value of each bin
            RC_values = [np.mean([bin_edges[i], bin_edges[i+1]]) for i in range(len(bin_edges) - 1)]
            probs = [0 for null in range(kwargs['n'])]    # initialize probabilities by bin

            # Check that we have at least two data points in this window
            if min(this_data[window_index]) == max(this_data[window_index]):
                raise RuntimeError('Found a window containing only a single sampled value. The boundaries of the offending '
                                   'window are: ' + str(windows[window_index]) + '. Sample more in this window or remove '
                                   'it from the input file.')

            # Build "probs" for this window
            for value in this_data[window_index]:
                reduced = (value - min(this_data[window_index]))/(max(this_data[window_index]) - min(this_data[window_index]))     # reduce to between 0 and 1
                local_index = int(np.floor(reduced * kwargs['n']))        # appropriate index within probs
                if local_index == kwargs['n']:
                    local_index -= 1                # handle case where reduced == 1
                probs[local_index] += 1             # increment probability count in the appropriate window

            # Scale probability counts to get fractions
            for i in range(len(probs)):
                probs[i] = probs[i]/len(this_data[window_index])
                if probs[i] == 0:
                    if not kwargs['slope_only']:
                        if cycle == kwargs['c']:
                            raise RuntimeError('at least one window contains a bin with zero samples. Either sample more in'
                                               ' this window or use fewer bins.\n The offending window is: ' +
                                               str(windows[window_index][0]) + ' to ' + str(windows[window_index][1]))
                        else:
                            raise RuntimeError('bootstrapping encountered a window containing a bin with zero samples. This'
                                               ' is most commonly caused by a window that has not been sampled across its '
                                               'full range of possible values. The offending window is: ' +
                                               str(windows[window_index][0]) + ' to ' + str(windows[window_index][1]) + '\n'
                                               'If you do not encounter a similar error when running without bootstrapping,'
                                               ' then you probably set the number of bootstrapping samples per window too '
                                               'low; increase it or set it to -1. Otherwise, please remove all of the data '
                                               'from this window from the input file.')



            # Calculate energy corresponding to each probability in this window
            U = [0 for null in range(kwargs['n'])]    # initialize energy values for this window
            local_index = 0
            offset = 0      # offset is used to anchor the first value in each window to zero before adjusting the whole window in the next step

            if not kwargs['slope_only']:
                for prob in probs:
                    if local_index == 0:
                        offset = -1 * kT * np.log(prob)
                    U[local_index] = -1 * kT * np.log(prob) - offset
                    local_index += 1

            else:
                this_result = optimize.minimize(objective_function, 0, (probs, RC_values, kT))
                for rc in RC_values:
                    if local_index == 0:
                        offset = this_result.x[0] * rc
                    U[local_index] = this_result.x[0] * rc - offset
                    local_index += 1

            # Adjust energy values uniformly up or down to stitch adjacent windows together smoothly
            if window_index == 0 or cycle < kwargs['c']:  # turn off boundary value matching during bootstrapping to avoid propagating errors in this step into PMF error in final step
                left_boundary = 0
            else:  # calculate the adjustment; left_boundary1 == left_boundary2 iff the slopes between the first two points of this window and the last two points of the previous window are equal; otherwise, we compromise to get the best fit
                f1 = (RC_values[0] - boundary_values[3]) / (boundary_values[2] - boundary_values[3])    # fraction of distance from last two points of previous window at which the first point of this window falls (e.g., 0.5 is exactly between, 1.2 is 20% further, etc.)
                left_boundary1 = ((boundary_values[0] - boundary_values[1]) * f1) + boundary_values[1]  # vertical shift to lower this window so that its first point intersects the line connecting the last two of the previous window
                f2 = (boundary_values[2] - RC_values[0]) / (RC_values[1] - RC_values[0])                # fraction of distance from last first two points of this window at which the last point of the previous window falls
                left_boundary2 = boundary_values[0] - ((U[1] - U[0]) * f2) + U[0]                       # vertical shift to lower this window so that the last point of the previous window intersects the line connecting the first two points of this window
                left_boundary = np.mean([left_boundary1, left_boundary2])                               # average of shift amounts (the amount we'll actually shift by)

            U = [item + left_boundary for item in U]    # apply the adjustment

            # Set boundary_values for next window
            boundary_values = [U[-1], U[-2], RC_values[-1], RC_values[-2]]  # [x2, x1, r2, r1]; store data for this step to calculate left_boundary for next step

            fullPMF += list(U)            # append the just-calculated data to the full energy profile
            fullRCs += list(RC_values)    # append list of RC values corresponding to fullPMF energy values

            if cycle == kwargs['c'] and not kwargs['noplot']:    # final cycle, so put together the first plot but don't show it yet
                ax0.errorbar(list(RC_values), list(U), std_err[error_index:error_index + len(U)])
                fig0.canvas.draw()
                nextcolor = list(colors.to_rgb(next(ax1._get_patches_for_fill.prop_cycler).get('color'))) + [0.75]
                ax1.bar(np.linspace(windows[window_index][0], windows[window_index][1], len(probs)), probs,
                        width=(RC_values[1] - RC_values[0]), color=nextcolor)
                fig1.canvas.draw()
                error_index += len(U)

        if cycle < kwargs['c']:
            subsampled_fullPMFs.append(fullPMF)     # this is the last line of this bootstrapping cycle
        elif not kwargs['noplot']:
            plt.show()

    # For smoothing the PMF into a single continuous line
    smoothPMF = []
    smoothRC = []
    smoothErr = []
    i = 0
    while i < len(fullPMF):
        if (i + 1) % kwargs['n'] == 0 and i + 1 < len(fullPMF):
            i += 1  # skip next point
        else:
            smoothPMF.append(fullPMF[i])
            smoothRC.append(fullRCs[i])
            smoothErr.append(std_err[i])
        i += 1

    with open(kwargs['o'], 'w') as f:
        for i in range(len(smoothRC)):
            f.write(str(smoothRC[i]) + ' ' + str(smoothPMF[i]) + ' ' + str(smoothErr[i]) + '\n')

    if not kwargs['noplot']:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(smoothRC,smoothPMF,color='#0072BD',lw=2)
        plt.fill_between(np.asarray(smoothRC), np.asarray(smoothPMF) - np.asarray(smoothErr), np.asarray(smoothPMF) + np.asarray(smoothErr),
                       alpha=0.5, facecolor='#0072BD')
        plt.ylabel('Free Energy (kcal/mol)', weight='bold')
        plt.xlabel('Reaction Coordinate', weight='bold')
        fig.canvas.draw()
        plt.show()


def update_progress(progress, message='Progress'):
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

    Returns
    -------
    None

    """
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
        status = "Done!\r\n"
    block = int(round(barLength * progress))
    text = "\r" + message + ": [{0}] {1}% {2}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 2), status)
    sys.stdout.write(text)
    sys.stdout.flush()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate free energy profile from the given equilibrium path sampling'
                                                 ' data')
    parser.add_argument('-i', metavar='input_file', type=str, nargs=1, default='eps.out',
                        help='input filename (output from equilibrium path sampling). Default=eps.out')
    parser.add_argument('-o', metavar='output_file', type=str, nargs=1, default='fep.out',
                        help='output filename. Default=fep.out')
    parser.add_argument('-t', metavar='temp', type=int, nargs=1, default=300,
                        help='temperature in Kelvin to evaluate energy at. Default=300')
    parser.add_argument('-n', metavar='nbins', type=int, nargs=1, default=5,
                        help='number of bins to divide each window into (at least 2). Default=5')
    parser.add_argument('-b', metavar='bootstrapN', type=int, nargs=1, default=-1,
                        help='number of bootstrapping samples to include in each window; if negative, then the length '
                             'of the full data set in that window is used instead. Default=-1')
    parser.add_argument('-c', metavar='bootstrapCyc', type=int, nargs=1, default=100,
                        help='number of bootstrapping cycles to average over. Default=100')
    parser.add_argument('-e', metavar='exclude_fract', type=float, nargs=1, default=0,
                        help='fraction of data from each window to exclude from calculations for decorrelation. Must be'
                             ' between 0 and 1. Default=0')
    parser.add_argument('--slope_only', action='store_true', default=False,
                        help='Evaluate energy profile using only a single line of slope that best explains the data '
                             'in that window (instead of using the Boltzmann equation). This option tolerates bins with'
                             ' no data in them, but ignores any changes in slope over the course of a window.')
    parser.add_argument('--noplot', action='store_true', default=False,
                        help='suppress free energy profile and window histogram plots')

    arguments = vars(parser.parse_args())  # Retrieves arguments as a dictionary object

    main(**arguments)
