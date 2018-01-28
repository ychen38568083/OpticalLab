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


#Dark Current
sumdark = np.zeros((1336, 2004))
for i in np.arange(136, 143):
	dark = pf.getdata("./10-18-2017/0" + str(i) + ".fts")
	sumdark = sumdark + dark

totaldark = sumdark/8


file = '10-18-2017/0128.fts'
#file = '10-11-2017/0009.fts'
header = pf.getheader(file) #Get the header from the fits file
ccd_data = pf.getdata(file) #Get the CCD data from the fits file



for key in header:
    print(key + ": ", header[key])
#########################################################################

plt.figure()

plt.imshow(np.flipud(np.rot90(ccd_data)), vmin = 500, vmax=800)

ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')

plt.show()
#####################################################################################################
ccd_data_1D = ccd_data.flatten()

plt.figure()

plt.hist(ccd_data_1D, bins=50)

plt.xlabel('Pixel Value [ADU]')
plt.ylabel('Number of Pixels')

plt.yscale("log", nonposy='clip')

plt.show()
####################################################################################################
plt.figure()

img = plt.imshow(np.flipud(np.rot90(ccd_data)), vmin=700, vmax=1000, cmap='gray_r')

clb = plt.colorbar(img, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')

ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.title('1')
plt.show()
###################################################################################################
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

dark_subtracted_ccd_data = ccd_data - totaldark
image = dark_subtracted_ccd_data - np.median(dark_subtracted_ccd_data)

star_radius = 18. # radius of region to blank out when removing star from image
gaussian_sigma = 4. # smoothing rms to use when search for stars
max_stars = 16 # maximum number of stars to locate

# find the brightest stars and return smoothed image

starxy, filtered_image = getstars(image, max_stars, star_radius, gaussian_sigma)

plt.figure()

img = plt.imshow(np.transpose(image), vmin=200, vmax=500, cmap="gray_r") #Set colormap

plt.plot(starxy[:,1], starxy[:,0],'ro', markerfacecolor='none', markersize=10)

clb = plt.colorbar(img, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')
#invert y-axis
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.title('2')
plt.show()
###################################################################################################################33

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

star_radius = 15. # radius of circular region to measure centroid
annulus_width = 4. # size of sky region Area(sky) = nsky x Area(star)

x_centroids, y_centroids, rms_x, rms_y, starflux = star_centroids(image, starxy, star_radius, annulus_width)

filtered_image = ndi.filters.gaussian_filter(image, gaussian_sigma) # filter the image

plt.figure()
img = plt.imshow(np.transpose(filtered_image), vmin=200, vmax=500, cmap="gray_r") #Set colormap

plt.plot(y_centroids, x_centroids, 'r+')

clb = plt.colorbar(img, fraction=0.03, pad=0.07)
clb.ax.set_ylabel('Counts [ADU]')
#invert y-axis
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.ylim(0, 2000)
plt.title('Centroids')
plt.show()

#######################################################################################################################

def usno(radeg,decdeg,fovam,epoch):
	import string as str
	import urllib.request as url

	str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1'
	str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)

	str = str1 + str2
	print('Calling Vizier', str)
	f= url.urlopen(str)

	# read from the obkect, sotring the page's contents in 's'.
	s = f.read()
	f.close()

	sl = s.splitlines()
	sl = sl[45:-1] # get rid of header 
	name = np.array([])
	rad = np.array([]) #RA in degrees
	ded = np.array([]) # DEC in degrees
	rmag = np.array([]) # rmage

	for k in sl:
		kw =k.split(b'\t')
		ded0 = float(kw[2])
		pmrad = float(kw[6])/3600e3/np.cos(np.deg2rad(ded0)) #convert from mas/yr to deg/year
		pmded = float(kw[7])/3600e3

		name = np.append(name,kw[0])
		rad = np.append(rad, float(kw[1]) + pmrad*(epoch-2000.0))
		ded = np.append(ded,float(kw[2]) + pmded*(epoch-2000.0))

		if kw[12] !=b'     ':
			#print('This is RMag,', kw[12])
			rmag = np.append(rmag,float(kw[12]))
		else:
			rmag = np.append(rmag,np.nan)
	return name,rad,ded,rmag




#######################################################################################################################


#ras = '20:59:54.03'
#des = '+07:16:50.6'

ras = '21:00:25.00'
des = '+07:19:35.0'

radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
dsgn = np.sign(float(des[0:3]))
dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.

fovam = 22.0
name, rad, ded, rmag = usno(radeg, dedeg, fovam, 2017)

#for i in np.arange(len(rmag)):
	#print('rmag= ',rmag[i])
	#print(rmag[i] < 12)

w = np.where(rmag < 12.5)[0]
#w = np.array([])
#for i in np.arange(len(rmag)):
	#elements = rmag[i]
	#w = np.append(w, elements < 12)


	#w = rmag[rmag < 12]	


plt.plot(rad[w],ded[w],'g.')
plt.locator_params(axis='x',nbins=4)
plt.locator_params(axis='y',nbins=4)
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]')
plt.title('USNO Starfield in RA and DEC')
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.show()

###########################################################

f = 6300# focal length in mm
p = 0.018 # pixel in mm
alpha = rad[w]*(np.pi/180) #degree to rad conversions
delta = ded[w]*(np.pi/180)
alpha_0= 315.104*(np.pi/180)
delta_0= 7.326*(np.pi/180)

X = -(np.cos(delta)*np.sin(alpha - alpha_0))/(np.cos(delta_0)*np.cos(delta)*np.cos(alpha - alpha_0) + np.sin(delta)*np.sin(delta_0))
Y = -(np.sin(delta_0)*np.cos(delta)*np.cos(alpha - alpha_0) - np.cos(delta_0)*np.sin(delta))/(np.cos(delta_0)*np.cos(delta)*np.cos(alpha - alpha_0) + np.sin(delta)*np.sin(delta_0))

y_0 = 1150
x_0 = 500	
print('X = ', X)
print('Y = ', Y)

converted_x = f*(X/p) + x_0
converted_y = f*(Y/p) + y_0

plt.plot(converted_x, converted_y, 'b.')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.title('USNO Starfield in Pixel Position')
plt.show()

###########################################################



plt.plot(y_centroids, x_centroids, 'r+', label = 'Centroids')
#invert y-axis
plt.plot(converted_x, converted_y, 'b.', label = 'USNO')
plt.xlabel('X [pixels]')
plt.ylabel('Y [pixels]')
plt.title('USNO Data and Centroids Comparison')
plt.legend()
plt.show()


def mymatch(x1,y1,x2,y2,dthres):

    # James R. Graham 11/7/2011 University of Toronto
    #
    # Match up two lists of stars
    #
    # (x1, y1) are the positions of stars in list 1
    # (x2, y2) are the positions of stars in list 2
    #
    # dthres is the star matching threshold distance in pixels
    #
    # id1 and id2 are the indices that match up lists 1 and 2
    # so that x1[id1] should match with x2[id2] etc.
    #

    import numpy as np

    dthres2=dthres**2 # distance squared

    # Step through each index of list 1, look for closest match in 2.

    nstar = np.size(x1)
    print (nstar)

    id1 = np.array([],dtype=int)
    id2 = np.array([],dtype=int)

    for i in np.arange(nstar):
        d2=(x1[i]-x2)**2+(y1[i]-y2)**2 # compute the distance squared

        m      = np.argmin(d2)
        dmatch = d2[m]

        if dmatch <= dthres2:
            # print np.sqrt(dmatch),i,m
            id1 = np.append(id1,i)
            id2 = np.append(id2,m)

    return id1,id2

a_xx = []
a_yy = []

a_xx.append(x_centroids[0])
a_xx.append(x_centroids[6])
a_xx.append(x_centroids[9])
a_xx.append(x_centroids[7])
a_xx.append(x_centroids[12])
a_xx.append(x_centroids[8])
a_xx.append(x_centroids[10])
a_xx.append(x_centroids[1])
a_xx.append(x_centroids[4])
a_xx.append(x_centroids[5])
a_xx.append(x_centroids[2])

a_yy.append(y_centroids[0])
a_yy.append(y_centroids[6])
a_yy.append(y_centroids[9])
a_yy.append(y_centroids[7])
a_yy.append(y_centroids[12])
a_yy.append(y_centroids[8])
a_yy.append(y_centroids[10])
a_yy.append(y_centroids[1])
a_yy.append(y_centroids[4])
a_yy.append(y_centroids[5])
a_yy.append(y_centroids[2])

a_x = np.array(a_xx)
a_y = np.array(a_yy)

aXX = []
aXX.append(X[18])
aXX.append(X[17])
aXX.append(X[7])
aXX.append(X[11])
aXX.append(X[6])
aXX.append(X[5])
aXX.append(X[4])
aXX.append(X[3])
aXX.append(X[16])
aXX.append(X[2])
aXX.append(X[10])

aYY = []
aYY.append(Y[18])
aYY.append(Y[17])
aYY.append(Y[7])
aYY.append(Y[11])
aYY.append(Y[6])
aYY.append(Y[5])
aYY.append(Y[4])
aYY.append(Y[3])
aYY.append(Y[16])
aYY.append(Y[2])
aYY.append(Y[10])

aX = np.array(aXX)
aY = np.array(aYY)

b1 = f*(np.array(aX)/p)
b2 = f*(np.array(aY)/p)
b3 = np.ones(len(b1))

B = np.column_stack((b1, b2, b3))

tb = np.transpose(B)
btb = np.dot(tb, B)
invbtb = np.linalg.inv(btb)
invbtbtb = np.dot(invbtb, tb)

cx = np.dot(invbtbtb, a_x)
cy = np.dot(invbtbtb, a_y)

#print(cy)
#print(cx)

t1 = []
t1.append((f/p)*cy[0])
t1.append((f/p)*cx[0])
t1.append(0)

t2 = []
t2.append((f/p)*cy[1])
t2.append((f/p)*cx[1])
t2.append(0)

t3 = []
t3.append(cy[2])
t3.append(cx[2])
t3.append(1)

T = np.column_stack((t1, t2, t3))
print(T)

Tdet = np.linalg.det(T)
sqrtdet = np.sqrt(Tdet)
print(sqrtdet)
print(sqrtdet/350000)

print(cy)
print(cx)
invT = np.linalg.inv(T)
xmatrix = []
xmatrix.append(converted_x)
xmatrix.append(converted_y)
xmatrix.append(1)


x_pix_usno, y_pix_usno, dum = np.dot(T, np.array([aX, aY,1]))
x_rms = np.sqrt( np.sum( (x_pix_usno - a_y)**2 ) / len(x_pix_usno))
y_rms = np.sqrt( np.sum( (y_pix_usno - a_x)**2 ) / len(y_pix_usno))
print('x rms = ', x_rms)
print('y rms = ', y_rms)

fx = x_pix_usno - a_y
fy = y_pix_usno - a_x
index = [2,5,6]

plt.plot(fx, 'r^', label = 'X')
plt.plot(fy, 'b^', label = 'Y')
plt.xlabel('Number of Stars Indices')
plt.ylabel('Error (Pixel)')
plt.ylim(-2.0, 2.0)
plt.legend()
plt.show()

nx_pix = np.delete(x_pix_usno, index)
ny_pix = np.delete(y_pix_usno, index)

nalx = np.delete(a_y, index)
naly = np.delete(a_x, index)

x_rms = np.sqrt( np.sum( (nx_pix - nalx)**2 ) / len(nx_pix))
y_rms = np.sqrt( np.sum( (ny_pix - naly)**2 ) / len(ny_pix))
print('x rms - outliers = ', x_rms)
print('y rms - outliers = ', y_rms)


nfx = np.delete(fx, index)
nfy = np.delete(fy, index)

plt.plot(nfx, 'r^', label = 'X')
plt.plot(nfy, 'b^', label = 'Y')
plt.xlabel('Number of Stars Indices')
plt.ylabel('Error (Pixel)')
plt.title('Linear Least-Squares Fit of Eight Stars')
plt.ylim(-1.0, 1.0)
plt.legend()
plt.show()


