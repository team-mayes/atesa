# A test script to generate spoofed aimless shooting data according to a fixed process, to help me investigate the
# relationship between the sampling procedure and the convergence behavior of the LMAX estimator variance.

import random

# The distribution is of a single variable x, where x can be between zero and one.
# The separatrix is orthogonal to the variable, located at x = 0.5.
# The distribution of samples is defined here:

# Unbiased, independent sampling
def unbiased():
    return random.random()

# Biased (gaussian), independent sampling
def biased():
    x = random.gauss(0.75,0.07)
    if x < 0:
        x = 0
    elif x > 1:
        x = 1
    return x

# Unbiased, correlated sampling
def unbias_corr(last_x):
    scale = 0.01
    delta = scale*random.random()
    sign = random.randint(0,1)
    if sign == 0:
        x = last_x + delta
    elif sign == 1:
        x = last_x - delta
    if x < 0:
        x = 0
    elif x > 1:
        x = 1
    return x

# Unbiased, anti-correlated sampling

data = []

x = 0.6 # initial condition, if appropriate

for i in range(25):
    for j in range(250):

        x = unbias_corr(x)

        basin = random.random()
        if basin > x:
            basin_name = 'A'
        else:
            basin_name = 'B'

        data.append(basin_name + ' <- ' + str(x) + '\n')

    open('spoofed_as_' + str(int((i + 1) * 250)) + '.out', 'w').close()
    with open('spoofed_as_' + str(int((i + 1) * 250)) + '.out', 'a') as f:
        for item in data:
            f.write(item)
