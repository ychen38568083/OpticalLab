import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt

def peak_finder(signal_array):
    threshold = 200
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

##################################################################################################################################
##################################################################################################################################
sumbg1 = np.zeros((1336, 2004))
for i in np.arange(1, 9):
	bg1 = pf.getdata("./10-10-2017/000" + str(i) + ".fts")
	sumbg1 = sumbg1 + bg1

sumbg2 = np.zeros((1336,2004))
for i in np.arange(10, 14, 1):
	bg2 = pf.getdata("./10-10-2017/00" + str(i) + ".fts")
	sumbg2 = sumbg2 + bg2

sumfocus1 = np.zeros((1336,2004))
for i in np.arange(78,97):
	focus = pf.getdata("./10-10-2017/00" + str(i) + ".fts")
	sumbfocus1 = sumfocus1 + focus
	
sumbg = (sumbg1 + sumbg2)/14
sumfocus = sumfocus1/10

print('sumbg = ', sumbg)
print('sumfocus =', sumfocus)
for i in np.arange(22,25):
	epsilonpeg = pf.getdata("./10-10-2017/00" + str(i) + ".fts")

epsilonpegmean = epsilonpeg/4 - sumfocus - sumbg

plt.imshow(epsilonpegmean, cmap='gray_r')
plt.show()


#data = pf.getdata('./10-10-2017/0017.fts')
#plt.imshow(data, cmap = 'gray')
#plt.colorbar()
#plt.show()
sumpho = np.zeros((1336, 2004))
for i in np.arange(32,37):
	pho = pf.getdata("./10-10-2017/00" + str(i) + ".fts") - sumbg - sumfocus
	sumpho = sumpho + pho
totalpho = sumpho/6

pixelx = np.arange(0, len(totalpho[:,0]))
pixely = np.arange(0, len(totalpho[0,:]))

plt.imshow(totalpho, cmap = 'gray_r')
plt.colorbar()
plt.show()

brightcol1 = totalpho[:,0]
brightrow1 = totalpho[0,:]


#print('First col values = ', brightcol1)
#print('Number of rows in 1st coloum = ', brightcol1.size)
#print('First row values = ', brightrow1)
#print('Number of cols in 1st row = ', brightrow1.size)

i = 0
cols = []
while i < len(totalpho[0,:]):
	colnumber = totalpho[:,i]
	#print('Coloum number', i,'=', colnumber)
	#cols.append(colnumber)
	#colcent = find_all_centroids(pixelx, colnumber)
	stackcol = np.vstack((colnumber, colnumber))
	cols.append(stackcol)
	#print(colcent)
	i = i + 1


i = 0
rows = []
while i < len(totalpho[:,0]):
	rownumber = totalpho[i,:]
	#print('Row number', i, '=', rownumber)
	rows.append(rownumber)
	#rowcent = find_all_centroids(pixely, rownumber)
	#print(rowcent)
	i = i + 1

#print('brightcols = ', cols)
#print('brightrows = ', rows)


###############################################################################################################################
brightcol = np.array(totalpho[:,1366])
brightrow = np.array(totalpho[1193,:])

peakbrightcol = np.array(peak_finder(brightcol))
peakbrightrow = np.array(peak_finder(brightrow))

centroidcol = find_all_centroids(pixelx,brightcol)
centroidrow = find_all_centroids(pixely,brightrow)



halfpeak = peakbrightrow/2
#print("Brightcol = ", brightcol)
#print("Peak Brightcol = ", peakbrightcol)
#print("Brightrow = ", brightrow)
#print("Peak Brightrow = ", peakbrightrow)

#print("Centroid col = ", centroidcol)
#print("Centroid row = ", centroidrow)

plt.plot(pixelx, brightcol)
plt.xlabel('Position in Rows')
plt.ylabel('Intensity')
plt.show()

plt.plot(pixely, brightrow)
plt.xlabel('Position in Coloum')
plt.ylabel('Intensity')
plt.show()





##plt.plot(brightrow)
#plt.title("birghtrow")
#plt.show()
#maxpho = totalpho.max()

#print(brightcol)
#print("Peak of Brightcol= ", peak_finder(brightcol))
#plt.plot(brightcol)
#plt.title("brightcol")
#plt.show()


