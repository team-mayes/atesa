import psutil
import pytraj

print('RAM memory % used:', psutil.virtual_memory()[2])

traj = pytraj.load('atesa/tests/test_data/test.rst7', 'atesa/tests/test_data/test.prmtop')

print('RAM memory % used:', psutil.virtual_memory()[2])

def dothing():
    distance = eval('pytraj.distance(traj, \'@105 @102\')[0]')
    return distance

for _ in range(100000):
    a = dothing()
    print('RAM memory % used:', psutil.virtual_memory()[2])