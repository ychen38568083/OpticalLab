import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sys import exit

rcParams['font.size'] = 15
rcParams['axes.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['legend.frameon'] = False

x = np.loadtxt('pmt_170901_1508_40.csv', delimiter = ',', dtype='int32')
x.shape
x[0,0]
x[0,1]
x[3:6,1]

t = x[:,1]
print(t)

plt.plot(t)
plt.xlabel('Event Number',fontsize = 17)
plt.ylabel('Clock Tick',fontsize = 17)

plt.show()

dt = t[1:] - t[0:-1]
print('dt is', dt)
plt.figure()
plt.plot(dt)
plt.xlabel('Event Number', fontsize = 17)
plt.ylabel('Interval [Clock Ticks]',fontsize = 17)

plt.show()

dt.size
plt.plot(dt, ',')
plt.xlabel('Event Number',fontsize = 17)
plt.ylabel('Clock Tick',fontsize = 17)
plt.plot()

plt.show()

i = np.arange(dt.size)
nstep = 1000
for j in i[0::nstep]:
	print (j, j+nstep-1, np.mean(dt[j:j+nstep]))

#figure 4
marr = np.array([])
i = np.arange(dt.size)
nstep = 1000
for j in i[0::nstep]:
	m = np.mean(dt[j:j+nstep])
	marr = np.append(marr,m)

plt.plot(marr, 'o')	
plt.ylabel('Mean Interval', fontsize = 17)
plt.xlabel('Start Index', fontsize = 17)
plt.show()

#figure 5
i = np.arange(dt.size)
nstep = 100
for j in i[0::nstep]:
	m = np.mean(dt[0:j+nstep])
	marr = np.append(marr,m)

plt.plot(marr, 'o')	
plt.ylabel('Mean Interval', fontsize = 17)
plt.xlabel('Number of Intervals averaged', fontsize =17)
plt.show()

#figure 6
marr2 = np.array([])
i = np.arange(dt.size)
nstep = 100
for j in i[0::nstep]:
	m = np.mean(dt[j:j+nstep])
	marr2 = np.append(marr2,m)

plt.plot(marr2, 'o')
plt.ylabel('Mean Interval',fontsize = 17)
plt.xlabel('Start Index',fontsize = 17)
plt.show()

#figure 7
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
#plt.show()
plt.close()

plt.figure(figsize=(10,7))
plt.plot(inv_sqrt_N, stds, '.')
plt.plot(inv_sqrt_N, s*inv_sqrt_N)
plt.xlabel(r'$ 1 / \sqrt{N}$')
plt.ylabel('Standard deviation of the mean (ticks)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(['Data', r'$ s / \sqrt{N} $'])
plt.savefig('std_rtN_corrected.png')
#plt.show()
plt.close()

#Histograms
N = 500
#define the lower and upper bin edges and bin width
bw = (dt.max() - dt.min())/(N-1.)
binl = dt.min() + bw * np.arange(N)

#define the array to hold the occurrence count
bincount = np.array([])
#loop through the bins
for bin in binl:
	count = np.where((dt >= bin) & (dt < bin+bw))[0].size
	bincount = np.append(bincount,count)

#compute the bin centers for plotting
#Figure 9
binc = binl + 0.5*bw
plt.figure()
plt.xlim(0, .7e8)
plt.ylabel('Frequency',fontsize = 17)
plt.xlabel('Interval [Ticks]',fontsize = 17)
plt.plot(binc, bincount, drawstyle='steps-mid')
plt.show()

#Figure 11-A
N = 100
#define the lower and upper bin edges and bin width
bw = (dt.max() - dt.min())/(N-1.)
binl = dt.min() + bw * np.arange(N)

#define the array to hold the occurrence count
bincount = np.array([])
#loop through the bins
for bin in binl:
	count = np.where((dt >= bin))[0].size
	bincount = np.append(bincount,count)

#compute the bin centers for plotting
binc = binl + 0.5*bw
p = (1/np.mean(binc))**(-1/np.mean(binc))
plt.figure()
plt.xlim(0, 1e8)
plt.ylabel('Frequency',fontsize = 17)
plt.xlabel('Interval [Ticks]',fontsize = 17)
plt.plot(binc, bincount, drawstyle='steps-mid')
plt.plot(p, '.')
plt.show()
print(p)

#Figure 11-B
plt.semilogy(binc, bincount, drawstyle='steps-mid')
plt.xlim(0,1e8)
plt.ylabel('Frequency',fontsize = 17)
plt.xlabel('Interval [Ticks]',fontsize = 17)
plt.show()


#3.3 figure 12 
x = np.loadtxt('pmt_170908_1827_40.csv', delimiter = ',', dtype='int32')
print(x.shape)
