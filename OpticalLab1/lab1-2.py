import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#3.3 figure 12 
b6 = np.loadtxt('pmt_170908_1827_40.csv', delimiter = ',', dtype='int32')
b5 = np.loadtxt('pmt_170908_1829_40.csv', delimiter = ',', dtype='int32')
b4 = np.loadtxt('pmt_170908_1831_40.csv', delimiter = ',', dtype='int32')
b3 = np.loadtxt('pmt_170908_1835_40.csv', delimiter = ',', dtype='int32')
b2 = np.loadtxt('pmt_170908_1837_40.csv', delimiter = ',', dtype='int32')
b1 = np.loadtxt('pmt_170908_1838_40.csv', delimiter = ',', dtype='int32')

t1 = b1[:,1]
t2 = b2[:,1]
t3 = b3[:,1]
t4 = b4[:,1]
t5 = b5[:,1]
t6 = b6[:,1]

dt1 = t1[1:] - t1[0:-1]

dt2 = t2[1:] - t2[0:-1]

dt3 = t3[1:] - t3[0:-1]

dt4 = t4[1:] - t4[0:-1]

dt5 = t5[1:] - t5[0:-1]

dt6 = t6[1:] - t6[0:-1]


mu6 = np.sum(dt6)/np.float(dt6.size)
std6 = np.sqrt(np.sum((dt6 - mu6)**2.)/(np.float(dt6.size)-1.))

mu5 = np.sum(dt5)/np.float(dt5.size)
std5 = np.sqrt(np.sum((dt5 - mu5)**2.)/(np.float(dt5.size)-1.))

mu4 = np.sum(dt4)/np.float(dt4.size)
std4 = np.sqrt(np.sum((dt4 - mu4)**2.)/(np.float(dt4.size)-1.))

mu3 = np.sum(dt3)/np.float(dt3.size)
std3 = np.sqrt(np.sum((dt3 - mu3)**2.)/(np.float(dt3.size)-1.))

mu2 = np.sum(dt2)/np.float(dt2.size)
std2 = np.sqrt(np.sum((dt2 - mu2)**2.)/(np.float(dt2.size)-1.))

mu1 = np.sum(dt1)/np.float(dt1.size)
std1 = np.sqrt(np.sum((dt1 - mu1)**2.)/(np.float(dt1.size)-1.))

t1m = np.mean(dt1)
t2m = np.mean(dt2)
t3m = np.mean(dt3)
t4m = np.mean(dt4)
t5m = np.mean(dt5)
t6m = np.mean(dt6)

print(std1, t1m, std6, t6m)
plt.plot([0, std6], [0, t6m], 'k-')
plt.plot(std6, t6m, 'o', std5, t5m, 'o', std4, t4m, 'o', std3, t3m, 'o', std2, t2m, 'o', std1, t1m, 'o')
plt.xlabel('Interval Sample Mean [Ticks]',fontsize = 17)
plt.ylabel('Interval Standard Deviation [Ticks]',fontsize = 17)
print(t1m)
print(b1.shape)
print(t1.shape)
plt.show()

x = np.loadtxt('pmt_170901_1508_40.csv', delimiter = ',', dtype='int32')
t = x[:,1]
dt = t[1:] - t[0:-1]
wa = np.where(dt > 3000)[0]
dt = dt[wa]

t1 = np.cumsum(dt)
dt1 = t1[1:] - t[0:-1]
print(dt1)
plt.plot(dt1)
plt.show()