import numpy as np
import astropy.io.fits as pf
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image
import pdb
import math

sun22 = pf.getdata("Nov-17-2017/solar_scan-035-2.fit")
sun35 = pf.getdata("Nov-17-2017/solar_scan-035-2.fit")
sun36 = pf.getdata("Nov-17-2017/solar_scan-036-2.fit")
sun37 = pf.getdata("Nov-17-2017/solar_scan-037-2.fit")
sun38 = pf.getdata("Nov-17-2017/solar_scan-038-2.fit")
sun39 = pf.getdata("Nov-17-2017/solar_scan-039-2.fit")
sun40 = pf.getdata("Nov-17-2017/solar_scan-040-2.fit")


sliced_sun37 = plt.imshow(sun37[660:703,:])
plt.show()
sun37_array = np.sum(sun37[660:703,:], axis = 0)
mean_sun37 = sun37_array

sliced_sum = np.sum(mean_sun37/ len(mean_sun37))
print(sliced_sum)

sub_sliced = mean_sun37 - sliced_sum

zeros1 = np.zeros(200)
zeros2 = np.zeros(200)

zero_sun1 = np.append(zeros1, sub_sliced)
zero_sun37 = np.append(zero_sun1, zeros2)
win_ham =np.hamming(len(zero_sun37))

##############################################################
sun22_array = np.sum(sun22[660:703,:], axis = 0)
mean_sun22 = sun22_array

sliced_sum22 = np.sum(mean_sun22/ len(mean_sun22))

sub_sliced22 = mean_sun37 - sliced_sum22

zeros22_1 = np.zeros(200)
zeros22_2 = np.zeros(200)

zero_sun1 = np.append(zeros22_1, sub_sliced)
zero_sun22 = np.append(zero_sun1, zeros22_2)

plt.plot(zero_sun37 * win_ham)
plt.plot(zero_sun22 * win_ham)
plt.show()