# Temporary script to automate collection of data for and execution of partial committor analysis runs.
# There should be no real use case for this for users; it is merely useful to me as a method developer to illustrate my
# new method!

from atesa import main
from atesa import utilities
import subprocess
import pickle
import numpy
import re
import sys
import os
import argparse
from jinja2 import Environment, FileSystemLoader

settings = pickle.load(open('settings.pkl', 'rb'))
allthreads = pickle.load(open('restart.pkl', 'rb'))

# Set Jinja2 environment (not stored in settings.pkl)
if os.path.exists(settings.path_to_templates):
    settings.env = Environment(loader=FileSystemLoader(settings.path_to_templates))
else:
    sys.exit('Error: could not locate templates folder: ' + settings.path_to_templates)

with open('comana_results.out', 'w') as f:
    f.write('data_length information_error comana_n comana_mean comana_std\n')

# Entire script will be written as a for loop over the points of interest (let's do every 750 frames, since I have 7500)
for data_length in range(750, 7500 + 750, 750):
    # First, perform LMAX and information error
    pattern = re.compile('[0-9.]+')     # pattern to match information error number in lmax output file
    settings.as_out_file = 'as_decorr_' + str(data_length) + '.out'

    command = 'lmax.py -i as_decorr_' + str(data_length) + '.out -q present --automagic -o ' + str(data_length) + '_comana_lmax.out'
    subprocess.check_call(command.split(' '), stdout=sys.stdout, preexec_fn=os.setsid)

    inf_err = float(pattern.findall(open(str(data_length) + '_comana_lmax.out', 'r').readlines()[-1])[0])
    open('comana_results.out', 'a').write(str(data_length) + ' ' + str(inf_err))

    rc = open(str(data_length) + '_comana_lmax.out', 'r').readlines()[1].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ','')
    cutoff_timestamp = int(open('as_raw_timestamped.out', 'r').readlines()[data_length - 1].split(' ')[0])

    # Then, loop over threads to build a list of every shooting point with timestamp less than or equal to the timestamp
    # for the data_length'th move and absolute RC value less than 0.05
    comana_init_coords = []
    for thread in allthreads:
        for move_index in range(len(thread.history.init_coords)):
            if thread.history.timestamps[move_index] <= cutoff_timestamp:
                cvs = utilities.get_cvs(thread.history.init_coords[move_index][0], settings, reduce=True).split(' ')
                rc_value = utilities.evaluate_rc(rc, cvs)
                if abs(rc_value) <= 0.05 and len(comana_init_coords) < 200:
                    comana_init_coords.append(thread.history.init_coords[move_index][0])

    # Finally, call committor analysis on the collected shooting points and collect the results into an output file
    comana_settings = argparse.Namespace()
    comana_settings.__dict__.update(settings.__dict__)
    comana_settings.job_type = 'committor_analysis'
    comana_settings.committor_analysis_n = 10
    comana_settings.committor_analysis_use_rc_out = False
    comana_settings.rc_threshold = 0.05
    comana_settings.initial_coordinates = comana_init_coords
    comana_settings.rc_definition = rc
    comana_settings.as_out_file = 'as_decorr_' + str(data_length) + '.out'
    comana_settings.rc_reduced_cvs = True
    comana_settings.restart = False
    comana_settings.overwrite = True
    comana_settings.working_directory = settings.working_directory + '/comana_' + str(data_length)
    comana_settings.dont_dump = True

    main.main(comana_settings)

    os.chdir(settings.working_directory)     # cd back to original working directory

    results_lines = open(comana_settings.working_directory + '/committor_analysis.out', 'r').readlines()[1:]
    total_committed = sum([int(string[string.index('/') + 1:]) for string in results_lines])
    evald = [eval(string) for string in results_lines]
    mean = numpy.mean(evald)
    std = numpy.std(evald)

    open('comana_results.out', 'a').write(' ' + str(total_committed) + ' ' + str(mean) + ' ' + str(std) + '\n')
