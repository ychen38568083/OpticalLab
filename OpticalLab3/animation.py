import numpy as np
import astropy.io.fits as pf
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image
import pdb

#ffmpeg -f image2 -r 1/5 -i image%05d.png -vcodec mpeg4 -y movie.mp4

#Dark Current
sumdark = np.zeros((1336, 2004))
for i in np.arange(136, 143):
	dark = pf.getdata("./10-18-2017/0" + str(i) + ".fts")
	sumdark = sumdark + dark

totaldark = sumdark/8

for i in np.arange(128,136):

	file = '10-18-2017/0' + str(i) + '.fts'
	header = pf.getheader(file) #Get the header from the fits file
	ccd_data = pf.getdata(file) - totaldark #Get the CCD data from the fits file
	plt.figure(figsize=(20,10))
	
	plt.imshow(ccd_data, vmin = 200, vmax=500)

	plt.xlabel('X [pixels]')
	plt.ylabel('Y [pixels]')

	plt.show()