
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.special import factorial

#load the data in clock ticks
b = np.loadtxt('pmt_170908_1838_40.csv', delimiter = ',', dtype=np.int32)
t = b[:,1]
dt = t[1:] - t[0:-1]
dt = dt[dt > 3000]
dt_float = dt.astype(float) #convert to float for cumsum
#wa is used to eliminate anything less than 3000 ticks within each time interval dt.
t1 = np.cumsum(dt_float)

plt.plot(t1[:2500])
plt.ylabel('Time[Clock Ticks]',fontsize = 17)
plt.xlabel('Event number',fontsize = 17)

plt.show()

#Counts per bin vs time: random distribution
n = 1500 #number of time bins
time_intervals = t1[-1] / n * np.arange(n)
counts = []
for i in range(len(time_intervals)-1):
	t_bin = [x for x in t1 if time_intervals[i] < x < time_intervals[i+1]]
	count = len(t_bin)
	counts.append(count)
counts = np.array(counts)
plt.plot(time_intervals[1:], counts)
plt.xlabel('Time[Ticks]',fontsize = 17)
plt.ylabel('Counts per bin',fontsize = 17)
plt.show()


#Histogram of counts per bin and frequency of occurance --poisson distribution
mu = np.mean(counts)
x = np.arange(15)

nbins = 16
binw = np.max(counts)/nbins
hist = np.histogram(counts, bins=nbins)
norm = np.sum(hist[0] * binw)

#poisson distribution
def poisson(x,mu):
	prob = np.exp(-1 * mu) * mu**x / factorial(x)
	return prob

plt.hist(counts, bins=nbins, histtype='step')
plt.plot(poisson(x,mu)*norm, 'k-')
plt.xlabel('Counts per bin', fontsize = 17)
plt.ylabel('Frequency', fontsize = 17)
plt.show()