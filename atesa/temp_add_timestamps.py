# Quick, temporary script for adding timestamps to existing aimless shooting data.
# This is just going to load in restart.pkl, go through each thread appending timestamps based on the modified time of
# each step's output file, and write a new restart.pkl (keep the old one for safety)

import os
import shutil
import sys
import pickle
from atesa.main import Thread

os.chdir('/oasis/scratch/comet/tburgin/temp_project/atesa_v2_working/')
shutil.copy('restart.pkl', 'restart_original.pkl')

allthreads = pickle.load(open('restart.pkl', 'rb'))

for thread in allthreads:
    thread.history.timestamps = []
    for trajs in thread.history.prod_trajs:
        outname = trajs[0].replace('.nc', '.out')
        if not os.path.exists(outname):
            raise FileNotFoundError('missing output file: ' + outname + '\naborting')
        thread.history.timestamps.append(int(os.path.getmtime(outname)))

# quick sanity check
for thread in allthreads:
    if not len(thread.history.prod_trajs) == len(thread.history.timestamps):
        raise RuntimeError('error in thread ' + thread.name + ': length of prod_trajs = ' + str(len(thread.history.prod_trajs)) + ' but lngth of timestamps = ' + str(len(thread.history.timestamps)))

pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)
