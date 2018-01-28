import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.signal import argrelextrema
import math

fluorpixelx = np.loadtxt('flour.txt')[:,0]
fluorpixely = np.loadtxt('flour.txt')[:,1]

#Background Noises in pixel and wavelength for calibration
##############################################################################################


bias1 = np.genfromtxt('./data/bias_3ms/bias00.txt',skip_header=17, skip_footer=1)[:,1]
bias2 = np.genfromtxt('./data/bias_3ms/bias01.txt',skip_header=17, skip_footer=1)[:,1]
bias3 = np.genfromtxt('./data/bias_3ms/bias02.txt',skip_header=17, skip_footer=1)[:,1]
bias4 = np.genfromtxt('./data/bias_3ms/bias03.txt',skip_header=17, skip_footer=1)[:,1]
bias5 = np.genfromtxt('./data/bias_3ms/bias04.txt',skip_header=17, skip_footer=1)[:,1]
bias6 = np.genfromtxt('./data/bias_3ms/bias05.txt',skip_header=17, skip_footer=1)[:,1]
bias7 = np.genfromtxt('./data/bias_3ms/bias06.txt',skip_header=17, skip_footer=1)[:,1]
bias8 = np.genfromtxt('./data/bias_3ms/bias07.txt',skip_header=17, skip_footer=1)[:,1]
bias9 = np.genfromtxt('./data/bias_3ms/bias08.txt',skip_header=17, skip_footer=1)[:,1]
bias10 = np.genfromtxt('./data/bias_3ms/bias09.txt',skip_header=17, skip_footer=1)[:,1]

meanbias = (bias1+bias2+bias3+bias4+bias5+bias6+bias7+bias8+bias9+bias10)/10
plt.title('Bias', fontsize = 17)
plt.ylabel('Intensity (Counts)', fontsize = 17)
plt.xlabel('Wavelength (NM)',fontsize = 17)
plt.ylim(80, 220)
plt.plot(meanbias)	
plt.show()

############################################################################################

#max1 = argrelextrema(data1y, np.greater,order = 50) #lists all the maximas in x indicies 
#min1 = argrelextrema(data1y, np.less, order = 50) #lists all the minimas in x indicies
#print('maximas:', max1)
#print('minimas:', min1)

#Centroids method
###############################################################################################################


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

##############################################################################################################
pixf = np.loadtxt('flour.txt')
wavef = np.loadtxt('Waveflour.txt')

#Fluor in pixel
wfluor1 = np.loadtxt("./data/fluorescent_22ms/fluorescent000.txt")
wfluor2 = np.loadtxt("./data/fluorescent_22ms/fluorescent001.txt")
wfluor3 = np.loadtxt("./data/fluorescent_22ms/fluorescent002.txt")
wfluor4 = np.loadtxt("./data/fluorescent_22ms/fluorescent003.txt")
wfluor5 = np.loadtxt("./data/fluorescent_22ms/fluorescent004.txt")
wfluor6 = np.loadtxt("./data/fluorescent_22ms/fluorescent005.txt")
wfluor7 = np.loadtxt("./data/fluorescent_22ms/fluorescent006.txt")
wfluor8 = np.loadtxt("./data/fluorescent_22ms/fluorescent007.txt")
wfluor9 = np.loadtxt("./data/fluorescent_22ms/fluorescent008.txt")
wfluor10 = np.loadtxt("./data/fluorescent_22ms/fluorescent009.txt")

wfluor11 = np.loadtxt("./data/fluorescent_22ms/fluorescent090.txt")
wfluor12 = np.loadtxt("./data/fluorescent_22ms/fluorescent091.txt")
wfluor13 = np.loadtxt("./data/fluorescent_22ms/fluorescent092.txt")
wfluor14 = np.loadtxt("./data/fluorescent_22ms/fluorescent093.txt")
wfluor15 = np.loadtxt("./data/fluorescent_22ms/fluorescent094.txt")
wfluor16 = np.loadtxt("./data/fluorescent_22ms/fluorescent095.txt")
wfluor17 = np.loadtxt("./data/fluorescent_22ms/fluorescent096.txt")
wfluor18 = np.loadtxt("./data/fluorescent_22ms/fluorescent097.txt")
wfluor19 = np.loadtxt("./data/fluorescent_22ms/fluorescent098.txt")
wfluor20 = np.loadtxt("./data/fluorescent_22ms/fluorescent099.txt")

neon0 = np.genfromtxt('./data/neon_16ms/red00.txt',skip_header=17, skip_footer=1)
neon1 = np.genfromtxt('./data/neon_16ms/red01.txt',skip_header=17, skip_footer=1)
neon2 = np.genfromtxt('./data/neon_16ms/red02.txt',skip_header=17, skip_footer=1)
neon3 = np.genfromtxt('./data/neon_16ms/red03.txt',skip_header=17, skip_footer=1)
neon4 = np.genfromtxt('./data/neon_16ms/red04.txt',skip_header=17, skip_footer=1)
neon5 = np.genfromtxt('./data/neon_16ms/red05.txt',skip_header=17, skip_footer=1)
neon6 = np.genfromtxt('./data/neon_16ms/red06.txt',skip_header=17, skip_footer=1)
neon7 = np.genfromtxt('./data/neon_16ms/red07.txt',skip_header=17, skip_footer=1)
neon8 = np.genfromtxt('./data/neon_16ms/red08.txt',skip_header=17, skip_footer=1)
neon9 = np.genfromtxt('./data/neon_16ms/red09.txt',skip_header=17, skip_footer=1)

neon10 = np.genfromtxt('./data/neon_16ms/red21.txt',skip_header=17, skip_footer=1)
neon11 = np.genfromtxt('./data/neon_16ms/red22.txt',skip_header=17, skip_footer=1)
neon12 = np.genfromtxt('./data/neon_16ms/red23.txt',skip_header=17, skip_footer=1)
neon13 = np.genfromtxt('./data/neon_16ms/red24.txt',skip_header=17, skip_footer=1)
neon14 = np.genfromtxt('./data/neon_16ms/red25.txt',skip_header=17, skip_footer=1)
neon15 = np.genfromtxt('./data/neon_16ms/red26.txt',skip_header=17, skip_footer=1)
neon16 = np.genfromtxt('./data/neon_16ms/red27.txt',skip_header=17, skip_footer=1)
neon17 = np.genfromtxt('./data/neon_16ms/red28.txt',skip_header=17, skip_footer=1)
neon18 = np.genfromtxt('./data/neon_16ms/red29.txt',skip_header=17, skip_footer=1)
neon19 = np.genfromtxt('./data/neon_16ms/red30.txt',skip_header=17, skip_footer=1)
neon20 = np.genfromtxt('./data/neon_16ms/red31.txt',skip_header=17, skip_footer=1)

neonpx = np.genfromtxt('Neon.txt')[:,0]

meanwf = (wfluor1 + wfluor2 + wfluor3 + wfluor4 + wfluor5 + wfluor6 + wfluor7 + wfluor8 + wfluor9 +wfluor10)/10
meanwfx = meanwf[:,0]
meanwfy = meanwf[:,1] - meanbias

plt.plot(meanwfx, meanwfy)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Intensity (Counts)', fontsize = 17)
plt.title('Bias Subtracted Fluorescent Spectra in Wavelength', fontsize = 17)
plt.show()
meanwf2 = (wfluor11 + wfluor12 + wfluor13 + wfluor14 + wfluor15 + wfluor16 + wfluor17 + wfluor18 + wfluor19 + wfluor20)/10
meanwfx2 = meanwf2[:,0]
meanwfy2 = meanwf2[:,1] - meanbias

plt.plot(fluorpixelx, meanwfy)
plt.xlabel('Pixel Number', fontsize = 17)
plt.ylabel('Intensity (Counts)', fontsize = 17)
plt.title('Bias Subtracted Fluorescent In Pixel Number', fontsize = 17)
plt.show()


#plt.plot(wavef[:,1])
#plt.plot(pixf[:,1])
#plt.title('2')
#plt.show()

#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#ax1.plot(x, y)
#ax1.set_title('Sharing Y axis')
#ax2.scatter(x, y)

#plt.plot(fluorpixelx, meanwfy2) #in pixel number
#plt.plot(meanwfx, meanwfy)
#plt.title('3')
#plt.show()
#print('Peaks for meanwfy:', peak_finder(meanwfy), 'Peaks for meanwfy2:', peak_finder(meanwfy2))

# centroids in wavelength
centroidswf = find_all_centroids(meanwfx, meanwfy)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Intensity(Counts)', fontsize = 17)
plt.title('Bias Subtracted Fluorescent Light in Wavelength w/centroids', fontsize = 17)
plt.plot(meanwfx, meanwfy)
for i in centroidswf:
	plt.axvline(i, ls='--')
plt.plot(300,0, ls='--', label='Centroid Location')
plt.legend(loc =2)
plt.show()


# centroids in pixel number
centroidsf = find_all_centroids(fluorpixelx,meanwfy2)
plt.xlabel('Pixel Number', fontsize = 17)
plt.ylabel('Intensity(Counts)', fontsize = 17)
plt.title('Bias Subtracted Fluorescent Light in Pixel Number w/centroids', fontsize = 17)
plt.plot(meanwfy2)

for i in centroidsf:
	plt.axvline(i, ls='--')
plt.plot(300,0, ls='--', label='Centroid Location')
plt.legend(loc =2)
plt.show()

#print(np.size(centroidsf), np.size(centroidsf))
#print('centroids in pixel:', centroidsf)
#print('centroids in wavelength:', centroidswf)

print(centroidswf)
print(centroidsf)
realCentp = []
realCentp.append(centroidsf[0])
realCentp.append(centroidsf[1])
realCentp.append(centroidsf[2])
realCentp.append(centroidsf[6])
realCentp.append(centroidsf[8])
realCentp.append(centroidsf[11])
realCentp.append(centroidsf[18])
print(realCentp)

realCentw = []
realCentw.append(centroidswf[0])
realCentw.append(centroidswf[1])
realCentw.append(centroidswf[2])
realCentw.append(centroidswf[9])
realCentw.append(centroidswf[10])
realCentw.append(centroidswf[12])
realCentw.append(centroidswf[19])
print(realCentw)

#coefficients = np.polyfit(realCentw, realCentp, 1)
fit, cov = np.polyfit(realCentw, realCentp, 1, full=False, cov = True)
p = np.poly1d(fit)
ys = p(realCentw)
plt.plot(realCentw, realCentp, 'o')
plt.title('Hg I Locations (Wavelength vs Pixel Number)', fontsize = 17)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Pixel Number', fontsize = 17)
plt.show()

plt.plot(realCentw, realCentp, 'o')
plt.plot(realCentw, ys, 'k-')
plt.title('Hg I Locations (Wavelength vs Pixel Number) w/Fit Line', fontsize = 17)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Pixel Number', fontsize = 17)
plt.show()

print('Covarience of Fluor:', cov)
print('Variance of Slope of Fluor: ', cov[0][0])
print('Standard Deviation of Slope of Fluor:', np.sqrt(cov[0][0]))

print('Variance of y-intercept of Fluor:', cov[1][1])
print('Standard Deviation of y-intercept of Fluor:', np.sqrt(cov[1][1]))	
print('slope (Gain):', p[1])
print('Y-intercept: ', p[0])
print('Coeff:', fit)	

###################################################################################NEON################################################################################


#var = G*(meanfy) + b

q= 1.6*10**(-16) #in culombs
c= 1*10**(-12) #in farady
g= 1/(160*10**(-9)) # in volts 
var = (q*g/c)*meanwfy 

print(var)
plt.plot(meanwfy, var, 'o')
plt.show()

vfit, vcov = np.polyfit(meanwfy, var, 1, full=False, cov = True)
vp = np.poly1d(vfit)
vys = vp(meanwfy)
plt.plot(meanwfy, var, 'o')
plt.plot(meanwfy, vys, 'k-')
plt.title('6')
plt.show()

print('Covarience:', vcov)
print('Variance of Slope : ', vcov[0][0])
print('Standard Deviation of Slope:', np.sqrt(vcov[0][0]))

print('Variance of y-intercept (Read Noise):', vcov[1][1])
print('Standard Deviation of y-intercept:', np.sqrt(vcov[1][1]))
print('slope (Gain):', vp[1])
print('Y-intercept: ', vp[0])
print('Coeff:', vfit)	

###################################################################3 #NEON ##############################################################################################
#########################################################################################################################################################################


mneonw = (neon0+neon1+neon2+neon3+neon4+neon5+neon6+neon7+neon8+neon9)/10
mneonwx = mneonw[:,0]
mneonwy = mneonw[:,1] - meanbias
mneonp = (neon10+neon11+neon12+neon13+neon14+neon15+neon16+neon17+neon18+neon19+neon20)/11
mneonpx = mneonp[:,0]
mneonpy = mneonp[:,1] - meanbias

wneoncents = find_all_centroids(mneonwx, mneonwy)
pneoncents = find_all_centroids(neonpx, mneonpy)

#print('Wavelength centroids: ', wneoncents)
plt.plot(mneonwx, mneonwy)
plt.xlabel('Wavelength(NM)')
plt.ylabel('Intensity(Counts)')
plt.title('Neon Light in Wavelength w/ centroids')
for i in wneoncents:
	plt.axvline(i, ls='--')
plt.show()

#print('Pixel Number centroids:', pneoncents)
plt.xlabel('Pixel Number')
plt.ylabel('Intensity(Counts)')
plt.title('Neon Light in Pixel Number w/ centroids')
plt.plot(neonpx, mneonpy)
for i in pneoncents:
	plt.axvline(i, ls='--')
plt.show()

neonrealCentw = []
neonrealCentw.append(wneoncents[0])
neonrealCentw.append(wneoncents[2])
neonrealCentw.append(wneoncents[4])
neonrealCentw.append(wneoncents[5])
neonrealCentw.append(wneoncents[8])
neonrealCentw.append(wneoncents[10])
neonrealCentw.append(wneoncents[11])
neonrealCentw.append(wneoncents[12])
neonrealCentw.append(wneoncents[13])
neonrealCentw.append(wneoncents[16])
print('Real Centroids in wavelength = ', neonrealCentw)

neonrealCentp = []
neonrealCentp.append(pneoncents[0])
neonrealCentp.append(pneoncents[2])
neonrealCentp.append(pneoncents[5])
neonrealCentp.append(pneoncents[6])#
neonrealCentp.append(pneoncents[9])
neonrealCentp.append(pneoncents[11])
neonrealCentp.append(pneoncents[12])
neonrealCentp.append(pneoncents[13])#
neonrealCentp.append(pneoncents[14])#
neonrealCentp.append(pneoncents[17])
print('Real Centroids in pixel number = ', neonrealCentp)

nfit, ncov = np.polyfit(neonrealCentw, neonrealCentp, 1, full=False, cov = True)
pn = np.poly1d(nfit)
nys = pn(neonrealCentw)
plt.plot(neonrealCentw, neonrealCentp, 'o')
plt.title('Ne I Locations (Wavelength vs Pixel Number)', fontsize = 17)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Pixel Number', fontsize = 17)
plt.show()

plt.plot(neonrealCentw, neonrealCentp, 'o')
plt.plot(neonrealCentw, nys, 'k-')
plt.title('Ne I Locations (Wavelength vs Pixel Number) w/Fit Line', fontsize = 17)
plt.xlabel('Wavelength (NM)', fontsize = 17)
plt.ylabel('Pixel Number', fontsize = 17)
plt.show()
print('Covarience for neon:', ncov)
print('Variance of Slope for neon: ', ncov[0][0])
print('Standard Deviation of Slope for neon:', np.sqrt(ncov[0][0]))

print('Variance of y-intercept for neon:', ncov[1][1])
print('Standard Deviation of y-intercept for neon:', np.sqrt(ncov[1][1]))
print('slope for neon:', pn[1])
print('Y-intercept: ', pn[0])
print('Coeff for neon:', nfit)	

q= 1.6*10**(-16) #in culombs
c= 1*10**(-12) #in farady
g= 1/(160*10**(-9)) # in volts 
nvar = (q*g/c)*mneonwy 

print(nvar)
plt.plot(mneonwy, nvar, 'o')
plt.show()
######################################################################################################
vnfit, vncov = np.polyfit(mneonwy, nvar, 1, full=False, cov = True)
vpn = np.poly1d(nfit)
vnys = pn(mneonwy)
plt.plot(mneonwy, nvar, 'o')
plt.plot(mneonwy, vnys, 'k-')
plt.title('6')
plt.show()

print('Covarience for neon:', vncov)
print('Variance of Slope for neon: ', vncov[0][0])
print('Standard Deviation of Slope for neon:', np.sqrt(vncov[0][0]))

print('Variance of y-intercept for neon:', vncov[1][1])
print('Standard Deviation of y-intercept for neon:', np.sqrt(vncov[1][1]))
print('slope for neon:', vpn[1])
print('Y-intercept: ', vpn[0])
print('Coeff for neon:', vnfit)

ptotalcentroids = []
ptotalcentroids = np.append(realCentp, neonrealCentp)
wtotalcentroids = []
wtotalcentroids = np.append(realCentw, neonrealCentw)
print('Total Centroids in pixel =', ptotalcentroids)
print('Total Centroids in wavelength =', wtotalcentroids)

tfit, cov = np.polyfit(wtotalcentroids, ptotalcentroids, 1, full=False, cov = True)
pt = np.poly1d(tfit)
tys = pt(wtotalcentroids)

######################################################################## Centroid line of Fluor
fit, cov = np.polyfit(realCentw, realCentp, 1, full=False, cov = True)
p = np.poly1d(fit)
ys = p(realCentw)
plt.plot(realCentw, ys, 'g', label = 'Hg I Linear Fit')
######################################################################## Centroid line of Neon
nfit, cov = np.polyfit(neonrealCentw, neonrealCentp, 1, full=False, cov = True)
pn = np.poly1d(nfit)
nys = pn(neonrealCentw)
plt.plot(neonrealCentw, nys,'g', label = 'Ne I Linear Fit')
######################################################################## Centroid line of Neon + Fluor
#plt.plot(wtotalcentroids, tys, 'k-')
plt.plot(wtotalcentroids, ptotalcentroids, 'o')
plt.legend(loc = 4)
plt.title('Total Centroids of Hg I and Ne I With Linear Fit', fontsize=17)
plt.xlabel('Centroids in Pixel Number', fontsize=17)
plt.ylabel('Centroids in Wavelength', fontsize=17)
plt.show()

tfit2, cov = np.polyfit(wtotalcentroids, ptotalcentroids, 2, full=False, cov = True)
pt2 = np.poly1d(tfit2)
tys2 = pt2(wtotalcentroids)

plt.plot(wtotalcentroids, tys2, 'k-')
plt.plot(wtotalcentroids, ptotalcentroids, 'o')
plt.title('Total Centroids of Fluorescent and Neon Lights With Quadratic Fit', fontsize=17)
plt.xlabel('Centroids in Pixel Number', fontsize=17)
plt.ylabel('Centroids in Wavelength', fontsize=17)
plt.show()

totalw = []
totalw = np.append(meanwfy, mneonwy)


q= 1.6*10**(-16) #in culombs
c= 1*10**(-12) #in farady
g= 1/(160*10**(-9)) # in volts 
tvar = (q*g/c)*totalw

tfit, tcov = np.polyfit(totalw, tvar, 1, full=False, cov = True)
tp = np.poly1d(tfit)
tys = tp(totalw)
plt.plot(totalw, tvar, 'o')
plt.title('Variance Versus Mean Wavelength', fontsize = 17)
plt.xlabel('Combined Mean Wavelength', fontsize = 17)
plt.ylabel('Variance', fontsize = 17)
plt.plot(totalw, tys, 'k-')
plt.show()

print('Covarience for neon:', tcov)
print('Variance of Slope for neon: ', tcov[0][0])
print('Standard Deviation of Slope for neon:', np.sqrt(tcov[0][0]))

print('Variance of y-intercept for neon:', tcov[1][1])
print('Standard Deviation of y-intercept for neon:', np.sqrt(tcov[1][1]))
print('slope for neon:', tp[1])
print('Y-intercept: ', tp[0])
print('Coeff for neon:', tfit)