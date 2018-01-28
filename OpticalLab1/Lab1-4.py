# AY120 Optical Lab #1, wk 1: Photon counting with an oscilliscope


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sys import exit

rcParams['font.size'] = 15
rcParams['axes.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['legend.frameon'] = False


# LOAD DATA (time in clock ticks)
data = np.loadtxt('pmt_170901_1508_40.csv', delimiter=',', dtype='int32', usecols=[1])  # if dtype is left to default, numbers are converted to floats in scientific notation: numpy.float64
dt = data[1:] - data[0:-1]  # could also make this a float here right off the bat

mu = np.mean(dt)
s = np.std(dt)


stepsizes = np.arange(50, 1000, 50)
stds = np.array([])
inv_sqrt_N = np.array([])

for j in range(10, 1000):
    nstep = j*1
    means = []
    for i in np.arange(0, len(dt), nstep): 
        mean_i = np.mean(dt[i:i+nstep])
        means.append(mean_i)
    means = np.array(means)
    stds = np.append(stds, np.std(means))
    inv_sqrt_N = np.append(inv_sqrt_N, 1 / np.sqrt(nstep))

plt.figure(figsize=(10,7))
plt.plot(stds, '.')
plt.xlabel('Number events averaged (N)')
plt.ylabel('Standard deviation of the mean (ticks)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig('std_dev_mean_corrected.png')
plt.show()
plt.close()

plt.figure(figsize=(10,7))
plt.plot(inv_sqrt_N, stds, '.')
plt.plot(inv_sqrt_N, s*inv_sqrt_N)
plt.xlabel(r'$ 1 / \sqrt{N}$')
plt.ylabel('Standard deviation of the mean (ticks)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(['Data', r'$ s / \sqrt{N} $'])
plt.savefig('std_rtN_corrected.png')
plt.show()
plt.close()