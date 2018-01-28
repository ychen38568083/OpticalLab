
import numpy as np
import astropy.io.fits as pf
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image
import pdb
import math

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (8,6)

def getstars(image, max_stars, star_radius, gaussian_sigma):
    #
    # James R. Graham 11/7/2011 University of Toronto
    #
    # image - image
    # max_stars - number of stars to find
    # star_radius - radius of exclusion zone
    # gaussian_sigma - sigma in the Gaussian smoothing function 
    #
    # Step 1) Find the stars in an image starting with the brightest peak.
    # Step 2) Remove that star and the search for the next brightest peak.

    filtered_image = ndi.filters.gaussian_filter(image, gaussian_sigma) # filter the image
    
    starxy = np.zeros([max_stars,2])

    x = np.arange(np.float(image.shape[1]))  # row/column order
    y = np.arange(np.float(image.shape[0]))
    xx, yy = np.meshgrid(x, y)

    for star in np.arange(max_stars):
        coordinate = np.unravel_index(np.argmax(filtered_image), image.shape)  # only pick out one value

        starxy[star,0] = coordinate[1]   # x - row/column order
        starxy[star,1] = coordinate[0]   # y 

        r2 = (xx-coordinate[1])**2. + (yy-coordinate[0])**2.
       
        filtered_image[np.where(r2 <= star_radius**2.)] = -999.0 # zero out that part of the array 

    return starxy, filtered_image

###################################################################################################################
#Dark Current
sumdark1 = np.zeros((1336, 2004))
for i in np.arange(16, 23):
	dark1 = pf.getdata("./10-11-2017/00" + str(i) + ".fts")
	sumdark1 = sumdark1 + dark1
avgdark1 = sumdark1/8

sumdark2 = np.zeros((1336, 2004))
for i in np.arange(136, 143):
	dark2 = pf.getdata("./10-18-2017/0" + str(i) + ".fts")
	sumdark2 = sumdark2 + dark2
avgdark2 = sumdark2/8
    
sumdark3 = np.zeros((1336, 2004))
for i in np.arange(17, 24):
	dark3 = pf.getdata("./10-26-2017/00" + str(i) + ".fts")
	sumdark3 = sumdark3 + dark3
avgdark3 = sumdark3/8

file1 = '10-11-2017/0013.fts'
file2 = '10-18-2017/0128.fts'
file3 = '10-26-2017/0008.fts'

header = pf.getheader(file1) #Get the header from the fits file
ccd_data1 = pf.getdata(file1) #Get the CCD data from the fits file
ccd_data2 = pf.getdata(file2)
ccd_data3 = pf.getdata(file3)

#for key in header:
    #print(key + ": ", header[key])
############################################################################################
#plt.figure()

img1 = plt.imshow(np.flipud(np.rot90(ccd_data1)), vmin = 700, vmax=1000)
clb = plt.colorbar(img1, fraction = 0.03, pad = 0.07)


ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
#plt.show()

#plt.imshow(np.flipud(np.rot90(ccd_data2)), vmin = 600, vmax=900)
clb = plt.colorbar(img1, fraction = 0.03, pad = 0.07)
ax = plt.gca()
ax.invert_yaxis()
plt.title('25 Phocaea Starfield (10-18-2017)')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
#plt.show()

#plt.imshow(np.flipud(np.rot90(ccd_data2 - avgdark2)), vmin = 200, vmax=500)
clb = plt.colorbar(img1, fraction = 0.03, pad = 0.07)
ax = plt.gca()
ax.invert_yaxis()
plt.title('Starfield Dark-Subtracted (10-18-2017)')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
#plt.show()


#plt.imshow(np.flipud(np.rot90(ccd_data3)), vmin = 500, vmax=800)
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
#plt.show()
###########################################################################################
darksubccd = ccd_data2 - avgdark2
ccd_data_1D = darksubccd.flatten()

#plt.figure()

plt.hist(ccd_data_1D, bins=50)

plt.xlabel('Pixel Value [ADU]')
plt.ylabel('Number of Pixels')
plt.title('Number of Pixels vs Pixel Values Histogram')
plt.yscale("log", nonposy='clip')

#plt.show()
#############################################################################################
img = plt.imshow(np.flipud(np.rot90(ccd_data2 - avgdark2)), vmin=500, vmax=800, cmap='gray_r')

clb = plt.colorbar(img, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')

ax = plt.gca()
ax.invert_yaxis()
plt.title('25 Phocaea Starfield (Gray)')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.show()

plt.figure()
####################################################################################3
dark_subtracted_ccd_data2 = ccd_data2 - avgdark2
image2 = dark_subtracted_ccd_data2 - np.median(dark_subtracted_ccd_data2)

star_radius = 18. # radius of region to blank out when removing star from image
gaussian_sigma = 4. # smoothing rms to use when search for stars
max_stars = 16 # maximum number of stars to locate

# find the brightest stars and return smoothed image

starxy2, filtered_image2 = getstars(image2, max_stars, star_radius, gaussian_sigma)

plt.figure()

img2 = plt.imshow(np.transpose(image2), vmin=200, vmax=500, cmap="gray_r") #Set colormap

plt.plot(starxy2[:,1], starxy2[:,0],'ro', markerfacecolor='none', markersize=10)

clb = plt.colorbar(img2, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')
#invert y-axis
ax = plt.gca()
ax.invert_yaxis()
plt.title('25 Phocaea with Star Locations (Circled)')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.ylim(0, 2000)
plt.show()

def star_centroids(image, starxy, star_radius, annulus_width):
    
    # James R. Graham 11/7/2011 University of Toronto
    #
    # Measure star centroids
    #
    # starrad  radius of sky aperture 
    # nsky     number of pixels insky annulus relative to star 

    x = np.arange(0, image.shape[1])  # row/column order
    y = np.arange(0, image.shape[0])

    sky_radius = np.sqrt(annulus_width+1)*star_radius

    xx, yy = np.meshgrid(x, y)

    x_centroids = np.array([]) # x-centroid
    y_centroids = np.array([]) # y-centroid
    starflux = np.array([]) # star counts in aperture
    rms_x = np.array([]) # rms width in x
    rms_y = np.array([]) # rms width in y

    i = 1
    for star in starxy: 

        r2 = (xx-star[0])**2. + (yy-star[1])**2.

        wstar = np.where(  r2 <= star_radius**2.)
        wsky  = np.where( (r2 >  star_radius**2.) & (r2 < sky_radius**2.) )

        # measure the centroid 

        medsky = np.median(image[wsky])

        # print 'Star %d'%i,star,' Median sky = ',medsky

        # compute the moments

        si   = np.sum((image[wstar] - medsky))
        six  = np.sum((image[wstar] - medsky)*xx[wstar])/si
        siy  = np.sum((image[wstar] - medsky)*yy[wstar])/si
        six2 = np.sum((image[wstar] - medsky)*xx[wstar]**2.)/si
        siy2 = np.sum((image[wstar] - medsky)*yy[wstar]**2.)/si

        rms_x     = np.append(rms_x, np.sqrt(six2 - six**2. )) 
        rms_y     = np.append(rms_y, np.sqrt(siy2 - siy**2. )) 
        x_centroids = np.append(x_centroids, six)
        y_centroids = np.append(y_centroids, siy)
        starflux = np.append(starflux,si) 
        i += 1

    return x_centroids, y_centroids, rms_x, rms_y, starflux

x_centroids, y_centroids, rms_x, rms_y, starflux = star_centroids(image2, starxy2, star_radius, annulus_width)

filtered_image = ndi.filters.gaussian_filter(image2, gaussian_sigma) # filter the image

plt.figure()
img2 = plt.imshow(np.transpose(filtered_image), vmin=200, vmax=500, cmap="gray_r") #Set colormap

plt.plot(y_centroids, x_centroids, 'r+')

clb = plt.colorbar(img2, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')
#invert y-axis
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.show()