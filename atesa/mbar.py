"""
Standalone script for converting umbrella sampling output files into free energy profiles via pymbar
"""

import os
import sys
import numpy as np, numpy
import scipy
import pymbar
import glob
import re
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


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


def main(**kwargs):
    """
    Main function for mbar.py. Collect all files in the present directory whose names begin with "rcwin_" and end with
    "us.dat", interpret them as umbrella sampling data centered on the RC value indicated by the first
    underline-separated float in the name, using the pymbar package and user-provided settings.

    Parameters
    ----------
    kwargs : dict
        Dictionary object containing arguments passed in from command line in "if __name__ == '__main__'" section.

    Returns
    -------
    None

    """

    # Compile and sort data files
    data_files = [name for name in glob.glob('rcwin_*.dat') if len(open(name, 'r').readlines()) > kwargs['min_data'][0]]
    pattern = re.compile('[-0-9.]+')
    data_files = sorted(data_files, key=lambda x: float(pattern.findall(x)[0]))

    if len(data_files) == 0:
        raise RuntimeError('did not find any data files matching the expected naming convention and with length of at '
                           'least min_data (' + str(kwargs['min_data'][0]) + ')')

    # Set thermodynamic parameters
    R = 1.987e-3                        # gas constant in kcal/mol-K
    beta = 1 / (R * kwargs['t'][0])     # inverse temperature of simulations (in 1/(kcal/mol))
    kconst = kwargs['k'][0]             # spring constant in kcal/mol

    # PMF range
    try:                    # try for user-defined rc_min
        rc_min = float(kwargs['rc_min'][0])
    except ValueError:      # default; use smallest window center
        rc_min = float(pattern.findall(data_files[0])[0])
    try:                    # try for user-defined rc_max
        rc_max = float(kwargs['rc_max'][0])
    except ValueError:      # default; use largest window center
        rc_max = float(pattern.findall(data_files[-1])[0])

    K = len(data_files)             # number of umbrellas   # todo: is there a use-case for making this user-defined?
    N_max = max(len(open(data_files[k], 'r').readlines()) for k in range(K))    # maximum number of snapshots/simulation

    rc_kn = numpy.zeros([K,N_max])  # rc_kn[k,n] is the RC for snapshot n from umbrella simulation k

    N_k = []    # number of samples from simulation k
    K_k = []    # spring constants
    rc0_k = []  # window centers

    # pymbar wants data in the form of a np.ndarray u_kn of shape [k, n] where there are k states and n samples
    # n.b. that u_kn has every sample in every row, regardless of the "source" window
    centers = []
    alldata = []
    if not kwargs['quiet']:
        count = 0
        update_progress(0, 'Collecting data')
    for k in range(K):
        lines = open(data_files[k], 'r').readlines()
        data = np.asarray([float(line.split()[1]) for line in lines][:])[kwargs['ignore'][0]:]

        if not kwargs['decorr'][0]:
            A_n = data
        else:
            # Use pymbar.timeseries to cut out pre-equilibrated data and decorrelate
            [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(data)  # compute indices of uncorrelated timeseries
            data_equil = data[t0:]
            indices = pymbar.timeseries.subsampleCorrelatedData(data_equil, g=g)
            A_n = data[indices]

        center = float(pattern.findall(data_files[k])[0])
        if not center in centers:
            centers.append(center)
            alldata.append(list(data))
        else:
            center_index = centers.index(center)
            alldata[center_index] += list(data)

        N_k.append(len(A_n))        # number of decorrelated samples in this window
        K_k.append(kconst)          # constant 20 kcal/mol spring constant
        rc0_k.append(float(pattern.findall(data_files[k])[0]))  # window center from file name
        for n in range(len(A_n)):
            rc_kn[k,n] = A_n[n]        # list of all decorrelated samples in this window

        if not kwargs['quiet']:
            count += 1
            update_progress(count / K, 'Collecting data')

    # print([rc0_k[k] for k in range(K)])
    # print([numpy.mean([item for item in rc_kn[k,:] if not item == 0]) - rc0_k[k] for k in range(K)])
    # print([scipy.stats.sem([item for item in rc_kn[k,:] if not item == 0]) for k in range(K)])

    print('checkpoint 0')

    # todo: figure out how to clean this particular plot up
    if not kwargs['quiet']:
        fig, ax = plt.subplots()
        ax.errorbar([rc0_k[k] for k in range(K)], [numpy.mean([item for item in rc_kn[k,:] if not item == 0]) - rc0_k[k] for k in range(K)], yerr=[scipy.stats.sem([item for item in rc_kn[k,:] if not item == 0]) for k in range(K)], color='#0072BD', lw=1)
        plt.ylabel('Mean Value', weight='bold')
        plt.xlabel('Reaction Coordinate', weight='bold')
        plt.title('THIS IS NOT YOUR FREE ENERGY PROFILE\nThis plot should probably be smooth.\nOtherwise, you may have some'
                  ' sampling abnormalities.')
        plt.show()

    open(kwargs['o'][0], 'w').close()   # initialize output file (and overwrite it if it exists)
    with open(kwargs['o'][0], 'a') as f:
        f.write('ATESA mbar.py output file. Skip to the end for the free energy profile if desired.\n')
        f.write('\n~~Mean Value Plot~~\n')
        f.write('THIS IS NOT YOUR FREE ENERGY PROFILE\nThis plot should probably be smooth.'
                '\nOtherwise, you may have some sampling abnormalities.\nSee ATESA\'s documentation for further '
                'information.\n')
        f.write('Reaction coordinate    Mean value\n')
        for k in range(K):
            f.write(str(rc0_k[k]) + '   ' + str(numpy.mean([item for item in rc_kn[k,:] if not item == 0])) + '\n')

    # Deprecated first-value plot
    # fig, ax = plt.subplots()
    # ax.plot([rc0_k[k] for k in range(K)], [rc_kn[k,0] for k in range(K)], color='#0072BD', lw=1)
    # plt.ylabel('First Value', weight='bold')
    # plt.xlabel('Reaction Coordinate', weight='bold')
    # plt.show()

    nbins = int(1.5 * len(set(rc0_k)))      # number of bins for 1D PMF (equal to 1.5 times the number of unique window centers)

    # Deprecated all-in-one histogram
    # fig, ax = plt.subplots()
    # n, bins, patches = ax.hist([item for item in np.reshape([rc_kn[k,:] for k in range(K)],-1) if not item == 0], int(10 * nbins))
    # plt.show()

    print('checkpoint 1')

    with open(kwargs['o'][0], 'a') as f:
        #f.write('\n~~Sampling Histogram~~\nEach line represents the samples in a single window\n')

        if not kwargs['quiet']:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for this_data in alldata:
            if not kwargs['quiet']:
                nextcolor = list(colors.to_rgb(next(ax._get_patches_for_fill.prop_cycler).get('color'))) + [0.75]
                ax.hist(this_data, list(numpy.arange(min(this_data), max(this_data) + 0.05, 0.05)), color=nextcolor)
                fig.canvas.draw()
            #f.write(str(this_data) + '\n')

        if not kwargs['quiet']:
            plt.xlabel('Reaction Coordinate', weight='bold')
            plt.ylabel('Number of Samples', weight='bold')
            plt.show()

    print('checkpoint 2')

    N_max = numpy.max(N_k)  # largest amount of data in a single simulation
    u_kn = numpy.zeros([K,N_max], numpy.float64)
    u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

    # Bin data
    delta = (rc_max - rc_min) / float(nbins)
    bin_center_i = numpy.zeros([nbins], numpy.float64)
    for i in range(nbins):
        bin_center_i[i] = rc_min + delta/2 + delta * i
    bin_kn = numpy.zeros([K,N_max], numpy.int32)
    for k in range(K):
        for n in range(N_k[k]):
            bin_kn[k,n] = int((rc_kn[k][n] - rc_min) / delta)   # compute bin assignment

    print('checkpoint 3')

    # Evaluate reduced energies in all umbrellas
    count = 0
    update_progress(0, 'Processing data')
    for k in range(K):
        for n in range(N_max):
            for l in range(K):
                # Compute deltaRC value of each point to this window center
                drc = rc_kn[k,n] - rc0_k[l]

                # Compute energy of point n interpreted as from window k with umbrella restraint l
                u_kln[k,l,n] = u_kn[k,n] + beta * K_k[k] * drc**2

        if not kwargs['quiet']:
            count += 1
            update_progress(count / K, 'Processing data')

    print('checkpoint 4')

    if not kwargs['quiet']:
        verbose = True
    else:
        verbose = False
    mbar = pymbar.mbar.MBAR(u_kln, N_k, verbose=verbose, maximum_iterations=1000)
    (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

    # Normalize by beta to convert from kT's of energy to kcal/mol of energy
    f_i /= beta
    df_i /= beta

    print('checkpoint 5')

    if not kwargs['quiet']:
        print("PMF (kcal/mol)")
        print("%8s %8s %8s" % ('bin', 'f', 'df'))
        for i in range(nbins):
            print("%8.2f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i]))

        fig, ax = plt.subplots()
        ax.plot([bin_center_i[i] for i in range(nbins)], [f_i[i] for i in range(nbins)], color='#0072BD', lw=2)
        plt.fill_between([bin_center_i[i] for i in range(nbins)], np.asarray([f_i[i] for i in range(nbins)]) - np.asarray([df_i[i] for i in range(nbins)]),
                         np.asarray([f_i[i] for i in range(nbins)]) + np.asarray([df_i[i] for i in range(nbins)]),
                         alpha=0.5, facecolor='#0072BD')
        plt.ylabel('Free Energy (kcal/mol)', weight='bold')
        plt.xlabel('Reaction Coordinate', weight='bold')
        plt.show()
    with open(kwargs['o'][0], 'a') as f:
        f.write('\n~~Free Energy Profile~~\nReaction coordinate bin center    Free energy (kcal/mol)    Error (kcal/mol)')
        for i in range(nbins):
            f.write(str(bin_center_i[i]) + '    ' + str(f_i[i]) + '    ' + str(df_i[i]) + '\n')


if __name__ == '__main__':

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(description='Evaluate free energy profile from the given umbrella sampling data')

    parser.add_argument('-t', metavar='temp', type=int, nargs=1, default=[300],
                        help='temperature in Kelvin to evaluate energy at. Default=300')
    parser.add_argument('-k', metavar='kconst', type=float, nargs=1, default=[50.0],
                        help='restraint weight in each window. The default matches the default value for umbrella '
                             'sampling in ATESA. Default=50.0')
    parser.add_argument('-o', metavar='output', type=str, nargs=1, default=['mbar.out'],
                        help='the name of the output file containing the data for the histogram, mean value, and free '
                             'energy plots. This file will be overwritten if it exists. Default=mbar.out')
    parser.add_argument('--min_data', metavar='min_data', type=int, nargs=1, default=[0],
                        help='minimum number of lines in a given data file for it to be eligible for inclusion in mbar.'
                             ' Default=0')
    parser.add_argument('--ignore', metavar='ignore', type=int, nargs=1, default=[1],
                        help='Number of samples from beginning of each data file to ignore in the analysis, as time to '
                             'decorrelate from initial coordinates. Generally should be used in a mutually exclusive '
                             'manner with the --decorr option. Default=1')
    parser.add_argument('--decorr', metavar='decorr', type=str2bool, nargs=1, default=[True],
                        help='if True, use pymbar.timeseries.detectEquilibration and '
                             'pymbar.timeseries.subsampleCorrelatedData to attempt to automatically use only '
                             'equilibrated and decorrelated data in the analysis. Default=True')
    parser.add_argument('--rc_min', metavar='rc_min', nargs=1, default=[''],
                        help='lower bound for range of RC values to include in the energy profile. The default setting '
                             'automatically uses the smallest window center value available, so only set this option if'
                             ' the default isn\'t working or if you want less than the full profile.')
    parser.add_argument('--rc_max', metavar='rc_max', nargs=1, default=[''],
                        help='upper bound for range of RC values to include in the energy profile. The default setting '
                             'automatically uses the largest window center value available, so only set this option if'
                             ' the default isn\'t working or if you want less than the full profile.')
    parser.add_argument('--quiet', action='store_true', default=False,
                        help='suppress all plots and progress bars. The only output will be the output file.')

    arguments = vars(parser.parse_args())  # Retrieves arguments as a dictionary object

    main(**arguments)
