import numpy as np
import astropy.io.fits as pf
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image
import pdb
import math

#'sun/solar_scan-021_12.fit'
sun22 = pf.getdata("Nov-17-2017/solar_scan-022-2.fit")
sun30 = pf.getdata("Nov-17-2017/solar_scan-030-2.fit")
sun31 = pf.getdata("Nov-17-2017/solar_scan-031-2.fit")
sun32 = pf.getdata("Nov-17-2017/solar_scan-032-2.fit")
sun33 = pf.getdata("Nov-17-2017/solar_scan-033-2.fit")
sun34 = pf.getdata("Nov-17-2017/solar_scan-034-2.fit")
sun35 = pf.getdata("Nov-17-2017/solar_scan-035-2.fit")
sun36 = pf.getdata("Nov-17-2017/solar_scan-036-2.fit")
sun37 = pf.getdata("Nov-17-2017/solar_scan-037-2.fit")
sun38 = pf.getdata("Nov-17-2017/solar_scan-038-2.fit")
sun39 = pf.getdata("Nov-17-2017/solar_scan-039-2.fit")
sun40 = pf.getdata("Nov-17-2017/solar_scan-040-2.fit")
sun41 = pf.getdata("Nov-17-2017/solar_scan-041-2.fit")
sun42 = pf.getdata("Nov-17-2017/solar_scan-042-2.fit")
sun43 = pf.getdata("Nov-17-2017/solar_scan-043-2.fit")
sun44 = pf.getdata("Nov-17-2017/solar_scan-044-2.fit")
sun45 = pf.getdata("Nov-17-2017/solar_scan-045-2.fit")
sun46 = pf.getdata("Nov-17-2017/solar_scan-046-2.fit")
sun47 = pf.getdata("Nov-17-2017/solar_scan-047-2.fit")
sun48 = pf.getdata("Nov-17-2017/solar_scan-048-2.fit")
sun49 = pf.getdata("Nov-17-2017/solar_scan-049-2.fit")
sun50 = pf.getdata("Nov-17-2017/solar_scan-050-2.fit")
sun51 = pf.getdata("Nov-17-2017/solar_scan-051-2.fit")
sun52 = pf.getdata("Nov-17-2017/solar_scan-052-2.fit")
sun53 = pf.getdata("Nov-17-2017/solar_scan-053-2.fit")
sun54 = pf.getdata("Nov-17-2017/solar_scan-054-2.fit")
sun55 = pf.getdata("Nov-17-2017/solar_scan-055-2.fit")
sun56 = pf.getdata("Nov-17-2017/solar_scan-056-2.fit")
sun57 = pf.getdata("Nov-17-2017/solar_scan-057-2.fit")
sun58 = pf.getdata("Nov-17-2017/solar_scan-058-2.fit")
dark_file = 'Scans/Dark.fit'
dark_data = pf.getdata(dark_file)

plt.figure()
plt.imshow(dark_data)
plt.show()

#########################################################################

def peak_finder(signal_array):
    threshold = 300
    peaks = []                            #x positions of the peaks, or rather, their index
    for i in range(2,len(signal_array)-2): 
        if signal_array[i] > signal_array[i+2]  and signal_array[i] > signal_array[i+1] and signal_array[i] > signal_array[i-1] and signal_array[i] > signal_array[i-2]:  #four conditions to be a peak (see description)
            if signal_array[i] > threshold:                      #is the value of the spectrum at i higher than our threshold?
                peaks.append(i)
    return peaks
    
#print('peaks at:',peak_finder(data1y))
#print('peak intensities are:', peak_finder(data1y))

def centroid(x_range,y_range):
    '''A function to return the centroid given equally sized x and y ranges over which to perform the calculation'''
    x_range = np.array(x_range) #make sure these are arrays if they aren't already
    y_range = np.array(y_range) #make sure these are arrays if they aren't already
    ... #convert the math formula for a centroid into code in these lines
    x_centroid = sum(x_range*y_range)/(sum(y_range))
    return x_centroid

def find_all_centroids(x_range,y_range):
	peaks = peak_finder(y_range) #define the peak positions in x indicies
	multicen = [] #empty array to append
	for i in peaks: #for loops for indicies in peaks
		y = y_range[i] #define the y which uses the y-axis indicies
		halfmax = y/2 #half of each peaks
		#print(halfmax)
		#multicen.append(centroid(x_range[i-4:i+4],y_range[i-4:i+4]))
		#The following codes are for more general way:
		dr = np.where(y_range[i:] < halfmax)[0][0] # everything to the right after half of each peaks
		dl = np.where(y_range[:i] < halfmax)[0][-1] # everything to the left
		multicen.append(centroid(x_range[dl:i+dr], y_range[dl:i+dr])) #append centroid back
	return multicen #returns multicen = [] with each newl updated centroid

neon_file = 'Scans/Neon.fit'
neon_header = pf.getheader(neon_file) #Get the header from the fits file
neon_data = pf.getdata(neon_file) #Get the CCD data from the fits file


for key in neon_header:
    print(key + ": ", neon_header[key])

plt.figure()

ds_neon = plt.imshow(neon_data - dark_data)

clb = plt.colorbar(ds_neon, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]', fontsize = 17)

plt.xlabel('X [pixels]', fontsize = 17)
plt.ylabel('Y [pixels]', fontsize = 17)
plt.title('Neon (Dark-subtracted)', fontsize = 17)
#600:647
plt.show()
sliced_neon = plt.imshow(neon_data[600:647,:] - dark_data[600:647,:])
plt.xlabel('X [pixels]', fontsize = 17)
plt.ylabel('Y [pixels]', fontsize = 17)
plt.title('Neon spectra slice from row 600 to 647', fontsize = 17)
plt.show()

print('neon_data =', neon_data)
print('Length of neon_data from 600 to 647 row =', len(neon_data[600:647,0]))

########################################################################
sum_array = np.sum(neon_data[600:647,:], axis = 0)
mean_sum = sum_array/47

print('sum of rows array = ', sum_array)
print('average of those 47 rows = ', mean_sum)
print(len(mean_sum))

plt.plot(sum_array)
#plt.show()
plt.title('1D Image of Echelle Order 6 w/ Centroids', fontsize = 17)
plt.xlabel('Pixels', fontsize = 17)
plt.ylabel('Intensity[Counts]', fontsize = 17)
centroids_neon = find_all_centroids(np.arange(0,len(mean_sum)),mean_sum)
plt.plot(centroids_neon,np.zeros(len(centroids_neon))+70000, 'r.')
plt.show()

print(centroids_neon)

wl_neon = np.array([609.616, 614.316, 626.649, 633.443, 640.225,650.654])
#wl_neon = np.array([638.299, 640.225, 650.653, 653.288, 659.895, 667.288])
nfit, ncov = np.polyfit(centroids_neon, wl_neon, 1, full=False, cov = True)
pn = np.poly1d(nfit)
nys = pn(centroids_neon)


plt.plot(centroids_neon, wl_neon, 'r.')
plt.plot(centroids_neon, nys, 'k-')
plt.title('Wavelength-Pixel Ratio w/fit line', fontsize = 17)
plt.xlabel('Pixels', fontsize = 17)
plt.ylabel('Wavelength (NM)', fontsize = 17)
plt.show()

plt.imshow(sun37)
plt.xlabel('Pixels', fontsize = 17)
plt.ylabel('Pixels', fontsize = 17)
plt.title('Solar Spectra 37', fontsize = 17)
plt.show()


sun_data1 = pf.getdata('Nov-17-2017/solar_scan-021-2.fit')
sum_sun1 = np.sum(sun_data1)
header = pf.getheader('Nov-17-2017/solar_scan-059-2.fit')
for key in header:
    print(key + ": ", header[key])
    
    
header2 = pf.getheader('Nov-17-2017/solar_scan-074-2.fit')
for key in header2:
    print(key + ": ", header2[key])
    
    
    
print(sun_data1)
print(sum_sun1)
#names = np.arange(1,70)
#strnames = str(names)

sun_intense1 = [] #sun scans from 1-9
for i in np.arange(1, 10):
    sun_datas = pf.getdata("Nov-17-2017/solar_scan-00" + str(i) + "-2.fit")
    sum_data = np.sum(np.mean(sun_datas))
    sun_intense1.append(sum_data)

sun_intense2 = [] #sun scans from 10 - 75
for i in np.arange(10, 75):
    sun_datas2 = pf.getdata("Nov-17-2017/solar_scan-0" + str(i) + "-2.fit")
    sum_data2 = np.sum(np.mean(sun_datas2))
    sun_intense2.append(sum_data2)
    
total_sun = sun_intense1 + sun_intense2
#print('total = ', total_sun)
plt.plot(total_sun, 'r.')
plt.xlabel('Time', fontsize = 17)
plt.ylabel('Average Intensity[Counts]', fontsize = 17)
plt.title('Solar Transition Data',fontsize=17)

sunfit = np.poly1d(total_sun)
plt.plot(sunfit)
plt.show()

#for i in np.arange(21, 60):
    #sun_datas2 = pf.getdata("Nov-17-2017/solar_scan-0" + str(i) + "-2.fit")
    #sum_data2 = np.sum(np.mean(sun_datas2))
    #print('number', i , ' = ', sum_data2)


sun37_array = np.sum(sun37[660:703,:], axis = 0)/43
sun38_array = np.sum(sun38[660:703,:], axis = 0)/43
sliced_sum37= np.sum(sun37_array/ len(sun37_array))

sub_sliced37 = sun37_array - sliced_sum37

plt.imshow(sun37[660:703,:])
plt.ylabel('Pixels', fontsize =17)
plt.xlabel('Pixels', fontsize =17)
plt.title('Slice of the Solar Spectra 37 from row 660 to 703', fontsize =17)
plt.show()

plt.plot(sun38_array)
plt.plot(sun37_array)
plt.ylabel('Intensity[Counts]', fontsize =17)
plt.xlabel('Pixels', fontsize =17)
plt.title('1D Image of Spectra 37 and 38', fontsize =17)
plt.show()



laser_file = 'Scans/650nmlaser.fit'
laser_header = pf.getheader(laser_file) #Get the header from the fits file
laser_data = pf.getdata(laser_file) #Get the CCD data from the fits file



#for key in laser_header:
    #print(key + ": ", laser_header[key])

plt.figure()

plt.imshow(laser_data - dark_data)

plt.xlabel('X [pixels]',fontsize=17)
plt.ylabel('Y [pixels]', fontsize=17)
plt.title('650 NM Laser (Dark-subtracted)',fontsize=17)

plt.show()
########################################################################

halogen_file = 'Scans/halogen.fit'
halogen_header = pf.getheader(halogen_file) #Get the header from the fits file
halogen_data = pf.getdata(halogen_file) #Get the CCD data from the fits file

halogen_file = 'Scans/halogen.fit'
halogen_header = pf.getheader(halogen_file) #Get the header from the fits file
halogen_data = pf.getdata(halogen_file) #Get the CCD data from the fits file



#for key in halogen_header:
    #print(key + ": ", halogen_header[key])

plt.figure()

plt.imshow(halogen_data - dark_data)

plt.xlabel('X [pixels]', fontsize=17)
plt.ylabel('Y [pixels]', fontsize=17)
plt.title('Quartz-Halogen (Dark-subtracted)',fontsize=17)

plt.show()