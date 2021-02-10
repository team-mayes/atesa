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
import copy
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

    # Format input directory string to ensure it has a trailing '/'
    input_path = kwargs['i'][0]
    if not input_path[-1] == '/':
        input_path += '/'

    # Compile and sort data files
    data_files = [name.replace(input_path, '') for name in glob.glob(input_path + 'rcwin_*_us.dat') if (len(open(name, 'r').readlines()) - 1) > kwargs['min_data'][0]]
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

    # Refresh data files to reflect new min and max, if necessary
    if not rc_max == float(pattern.findall(data_files[-1])[0]) or not rc_min == float(pattern.findall(data_files[0])[0]):
        data_files = [name.replace(input_path, '') for name in glob.glob(input_path + 'rcwin_*_us.dat') if rc_min <= float(pattern.findall(name)[0]) <= rc_max]
        data_files = sorted(data_files, key=lambda x: float(pattern.findall(x)[0]))

    # Combine data from identical states into temp files, implementing decorr and ignore in the process
    if not kwargs['quiet']:
        count = 0
        update_progress(0, 'Collecting data')
    data_files_temp = []
    for data_file in data_files:
        center = str(pattern.findall(data_file)[0])
        if not os.path.exists('mbar_temp_' + center + '_0.dat'):
            open('mbar_temp_' + center + '_0.dat', 'w').close()
        lines = open(data_file, 'r').readlines()[kwargs['ignore'][0] + 1:]
        data = np.asarray([float(line.split()[1]) for line in lines if len(line.split()) > 1])

        if kwargs['decorr']:
            [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(data)  # compute indices of uncorrelated timeseries
            data_equil = data[t0:]
            indices = pymbar.timeseries.subsampleCorrelatedData(data_equil, g=g)
            A_n = data[indices]
        else:
            A_n = data
        with open('mbar_temp_' + center + '_0.dat', 'a') as f:
            for item in A_n:
                f.write(str(item) + '\n')
        if not 'mbar_temp_' + center + '_0.dat' in data_files_temp:
            data_files_temp.append('mbar_temp_' + center + '_0.dat')
        if not kwargs['quiet']:
            count += 1
            update_progress(count / len(data_files), 'Collecting data')
    data_files = copy.deepcopy(data_files_temp)

    K = len(data_files)             # number of umbrellas
    N_max = max(len(open(data_files[k], 'r').readlines()) for k in range(K))    # maximum number of snapshots/simulation

    rc_kn = numpy.zeros([K,N_max])  # rc_kn[k,n] is the RC for snapshot n from umbrella simulation k

    N_k = []    # number of samples from simulation k
    K_k = []    # spring constants
    rc0_k = []  # window centers

    # pymbar wants data in the form of a np.ndarray u_kn of shape [k, n] where there are k states and n samples
    # n.b. that u_kn has every sample in every row, regardless of the "source" window
    centers = []
    alldata = []
    for k in range(K):
        center = float(pattern.findall(data_files[k])[0])

        if rc_min <= center <= rc_max:
            lines = open(data_files[k], 'r').readlines()
            A_n = np.asarray([float(line) for line in lines])

            if not center in centers:
                centers.append(center)
                alldata.append(list(A_n))
            else:
                center_index = centers.index(center)
                alldata[center_index] += list(A_n)

            N_k.append(len(A_n))        # number of decorrelated samples in this window
            K_k.append(kconst)          # constant spring constant
            rc0_k.append(float(pattern.findall(data_files[k])[0]))  # window center from file name
            for n in range(len(A_n)):
                rc_kn[k,n] = A_n[n]        # list of all decorrelated samples in this window

    # Now that all of the necessary data is loaded into memory, delete the temp files
    for filename in data_files_temp:
        os.remove(filename)

    # Deprecated printing of data for debugging
    # print([rc0_k[k] for k in range(K)])
    # print([numpy.mean([item for item in rc_kn[k,:] if not item == 0]) - rc0_k[k] for k in range(K)])
    # print([scipy.stats.sem([item for item in rc_kn[k,:] if not item == 0]) for k in range(K)])

    if not kwargs['quiet']:
        fig, ax = plt.subplots()
        ax.errorbar([rc0_k[k] for k in range(K)], [numpy.mean([item for item in rc_kn[k,:] if not item == 0]) - rc0_k[k] for k in range(K)], yerr=[scipy.stats.tstd([item for item in rc_kn[k,:] if not item == 0]) for k in range(K)], color='#0072BD', lw=1)
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
            f.write(str(rc0_k[k]) + '   ' + '%.3f' % numpy.mean([item for item in rc_kn[k,:] if not item == 0]) + '\n')

    # Deprecated first-value plot
    # fig, ax = plt.subplots()
    # ax.plot([rc0_k[k] for k in range(K)], [rc_kn[k,0] for k in range(K)], color='#0072BD', lw=1)
    # plt.ylabel('First Value', weight='bold')
    # plt.xlabel('Reaction Coordinate', weight='bold')
    # plt.show()

    # Set number of bins for 1D PMF (equal to 1.5 times the number of unique window centers)
    # This is somewhat arbitrary and can be changed if desired, but there should be no need
    nbins = int(1.5 * len(set(rc0_k)))

    # Deprecated all-in-one histogram
    # fig, ax = plt.subplots()
    # n, bins, patches = ax.hist([item for item in np.reshape([rc_kn[k,:] for k in range(K)],-1) if not item == 0], int(10 * nbins))
    # plt.show()

    with open(kwargs['o'][0], 'a') as f:
        f.write('\n~~Sampling Histogram~~\nEach pair of lines is bin edges followed by bin heights\n')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for this_data in alldata:
            nextcolor = list(colors.to_rgb(next(ax._get_patches_for_fill.prop_cycler).get('color'))) + [0.75]
            n, bins, null = ax.hist(this_data, list(numpy.arange(min(this_data), max(this_data) + 0.05, 0.05)), color=nextcolor)
            if not kwargs['quiet']:
                fig.canvas.draw()
            f.write(' '.join(['%.3f' % item for item in bins]) + '\n')
            f.write(' '.join([str(int(item)) for item in n]) + '\n\n')

        if not kwargs['quiet']:
            plt.xlabel('Reaction Coordinate', weight='bold')
            plt.ylabel('Number of Samples', weight='bold')
            plt.show()

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

    if not kwargs['quiet']:
        verbose = True
    else:
        verbose = False
    mbar = pymbar.mbar.MBAR(u_kln, N_k, verbose=verbose, maximum_iterations=1000)
    (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

    # Normalize by beta to convert from kT's of energy to kcal/mol of energy
    f_i /= beta
    df_i /= beta

    if not kwargs['quiet']:
        print('\nFinal results. This analysis was based on the pyMBAR package. If you publish work based on this '
              'output, in addition to ATESA please cite:\n')
        print('Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states. J.'
              ' Chem. Phys. 129:124105 (2008). DOI: 10.1063/1.2978177\n')
        if kwargs['decorr']:
            print('Since this run also included the "decorr" option, you also need to cite:\n')
            print('Chodera JD. A simple method for automated equilibration detection in molecular simulations. J. Chem.'
                  ' Theor. Comput. 12:1799, 2016. DOI: 10.1021/acs.jctc.5b00784')
        print("\nPMF (kcal/mol)")
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
        f.write('\nFinal results. This analysis was based on the pyMBAR package. If you publish work based on this '
                'output, in addition to ATESA please cite:\n')
        f.write('Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states.'
                ' J. Chem. Phys. 129:124105 (2008). DOI: 10.1063/1.2978177\n')
        if kwargs['decorr']:
            f.write('Since this run also included the "decorr" option, you also need to cite:\n')
            f.write('Chodera JD. A simple method for automated equilibration detection in molecular simulations. J. '
                    'Chem. Theor. Comput. 12:1799, 2016. DOI: 10.1021/acs.jctc.5b00784')
        f.write('\n~~Free Energy Profile~~\nReaction coordinate    Free energy (kcal/mol)    Error (kcal/mol)\n')
        for i in range(nbins):
            f.write('%.3f' % bin_center_i[i] + '    ' + '%.3f' % f_i[i] + '    ' + '%.3f' % df_i[i] + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate free energy profile from the given umbrella sampling data')

    parser.add_argument('-t', metavar='temp', type=int, nargs=1, default=[300],
                        help='temperature in Kelvin to evaluate energy at. Default=300')
    parser.add_argument('-k', metavar='kconst', type=float, nargs=1, default=[50.0],
                        help='restraint weight in each window. The default matches the default value for umbrella '
                             'sampling in ATESA. Default=50.0')
    parser.add_argument('-o', metavar='output', type=str, nargs=1, default=['mbar.out'],
                        help='the name of the output file containing the data for the histogram, mean value, and free '
                             'energy plots. This file will be overwritten if it exists. Default=mbar.out')
    parser.add_argument('-i', metavar='input_path', type=str, nargs=1, default=['./'],
                        help='the path to the directory containing the data files. Both relative and absolute paths are'
                             ' accepted. Default=./')
    parser.add_argument('--min_data', metavar='min_data', type=int, nargs=1, default=[0],
                        help='minimum number of lines in a given data file for it to be eligible for inclusion in mbar.'
                             ' Default=0')
    parser.add_argument('--ignore', metavar='ignore', type=int, nargs=1, default=[1],
                        help='Number of samples from beginning of each data file to ignore in the analysis, as time to '
                             'decorrelate from initial coordinates. Generally should be used in a mutually exclusive '
                             'manner with the --decorr option. Default=1')
    parser.add_argument('--decorr', action='store_true', default=False,
                        help='use pymbar.timeseries.detectEquilibration and pymbar.timeseries.subsampleCorrelatedData '
                             'to attempt to automatically use only equilibrated and decorrelated data in the analysis.')    # todo: make decorr the default, change this to --no-decorr?
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
