import subprocess
import os

command = 'python temp.py'
process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, preexec_fn=os.setsid)
