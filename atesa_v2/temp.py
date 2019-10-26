from atesa_v2 import utilities
from atesa_v2 import information_error
from main import Thread
import os
import shutil

os.chdir('tests/test_data')
shutil.copy('restart.pkl', '../test_temp')

information_error.main('as_big_123.out')
