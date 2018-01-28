import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pf

################################################################################################################################################################################################################################
#Find the peaks 
def peak_finder(signal_array):
    threshold = 300
    peaks = []                                                    
    for i in range(2,len(signal_array)-2): 
        if signal_array[i] > signal_array[i+1] and signal_array[i] > signal_array[i-1] and signal_array[i] > signal_array[i-2] and signal_array[i] > signal_array[i+2]:  
            if signal_array[i] > threshold:                     
                peaks.append(i)
    return peaks
#print('peaks = ',peak_finder(FlorIntents))

def centroid(x_range,y_range):
    '''A function to return the centroid given equally sized x and y ranges over which to perform the calculation'''
    x_range = np.array(x_range)
    y_range = np.array(y_range) 
    x = np.sum(x_range*y_range) 
    y = np.sum(y_range)
    x_centroid = x/y
    return x_centroid


def find_all_centroids(x_range,y_range):
	peaks = peak_finder(y_range) #define the peak positions in x indicies
	#print(peaks)
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
################################################################################################################################################################################################################################



z1 = np.zeros(100)
z2 = np.zeros(100)
alcents = []


file37 = './Nov-17-2017/solar_scan-037-2.fit'

header37 = pf.getheader(file37) #Get the header from the fits file
data37 = pf.getdata(file37) #Get the CCD data from the fits file

sliced37 = data37[660:703,:]
od37 = np.sum(sliced37, axis = 0)
mean37 = od37/43

sliced_av37 = np.sum(mean37/(len(mean37)))
Nsum37 = mean37 - sliced_av37

z1 = np.zeros(100)
z2 = np.zeros(100)

asum37 = np.append(z1, Nsum37)
NewSum37 = np.append(asum37, z2)

ham37 = np.hamming(len(NewSum37))
win = ham37 * NewSum37


################################################################################################################
file21 = './Nov-17-2017/solar_scan-021-2.fit'
header21 = pf.getheader(file21) #Get the header from the fits file
data21 = pf.getdata(file21) #Get the CCD data from the fits file

sliced21 = data21[660:703,:]
od21 = np.sum(sliced21, axis = 0)
mean21 = od21/43

sliced_av21 = np.sum(mean21/(len(mean21)))
Nsum21 = mean21 - sliced_av21

asum21 = np.append(z1, Nsum21)
NewSum21 = np.append(asum21, z2)

ham21 = np.hamming(len(NewSum21))
win21 = ham21 * NewSum21

al = []
for i in np.arange(-15, 15):
    rwin21 = np.roll(win21, i)
    totalsum = np.sum(win * rwin21)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents21 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents21)

################################################################################################################
file22 = './Nov-17-2017/solar_scan-022-2.fit'
header22 = pf.getheader(file22) #Get the header from the fits file
data22 = pf.getdata(file22) #Get the CCD data from the fits file

sliced22 = data22[660:703,:]
od22 = np.sum(sliced22, axis = 0)
mean22 = od22/43

sliced_av22 = np.sum(mean22/(len(mean22)))
Nsum22 = mean22 - sliced_av22

asum22 = np.append(z1, Nsum22)
NewSum22 = np.append(asum22, z2)

ham22 = np.hamming(len(NewSum22))
win22 = ham22 * NewSum22

al = []
for i in np.arange(-15, 15):
    rwin22 = np.roll(win22, i)
    totalsum = np.sum(win * rwin22)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents22 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents22)
################################################################################################################
file23 = './Nov-17-2017/solar_scan-023-2.fit'
header23 = pf.getheader(file23) #Get the header from the fits file
data23 = pf.getdata(file23) #Get the CCD data from the fits file

sliced23 = data23[660:703,:]
od23 = np.sum(sliced23, axis = 0)
mean23 = od23/43

sliced_av23 = np.sum(mean23/(len(mean23)))
Nsum23 = mean23 - sliced_av23

asum23 = np.append(z1, Nsum23)
NewSum23 = np.append(asum23, z2)

ham23 = np.hamming(len(NewSum23))
win23 = ham23 * NewSum23

al = []
for i in np.arange(-15, 15):
    rwin23 = np.roll(win23, i)
    totalsum = np.sum(win * rwin23)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents23 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents23)
################################################################################################################
file24 = './Nov-17-2017/solar_scan-024-2.fit'
header24 = pf.getheader(file24) #Get the header from the fits file
data24 = pf.getdata(file24) #Get the CCD data from the fits file

sliced24 = data24[660:703,:]
od24 = np.sum(sliced24, axis = 0)
mean24 = od24/43

sliced_av24 = np.sum(mean24/(len(mean24)))
Nsum24 = mean24 - sliced_av24

asum24 = np.append(z1, Nsum24)
NewSum24 = np.append(asum24, z2)

ham24 = np.hamming(len(NewSum24))
win24 = ham24 * NewSum24

al = []
for i in np.arange(-15, 15):
    rwin24 = np.roll(win24, i)
    totalsum = np.sum(win * rwin24)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents24 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents24)
################################################################################################################
file25 = './Nov-17-2017/solar_scan-025-2.fit'
header25 = pf.getheader(file25) #Get the header from the fits file
data25 = pf.getdata(file25) #Get the CCD data from the fits file

sliced25 = data25[660:703,:]
od25 = np.sum(sliced25, axis = 0)
mean25 = od25/43

sliced_av25 = np.sum(mean25/(len(mean25)))
Nsum25 = mean25 - sliced_av25

asum25 = np.append(z1, Nsum25)
NewSum25 = np.append(asum25, z2)

ham25 = np.hamming(len(NewSum25))
win25 = ham25 * NewSum25

al = []
for i in np.arange(-15, 15):
    rwin25 = np.roll(win25, i)
    totalsum = np.sum(win * rwin25)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents25 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents25)
################################################################################################################
file26 = './Nov-17-2017/solar_scan-026-2.fit'
header26 = pf.getheader(file26) #Get the header from the fits file
data26 = pf.getdata(file26) #Get the CCD data from the fits file

sliced26 = data26[660:703,:]
od26 = np.sum(sliced26, axis = 0)
mean26 = od26/43

sliced_av26 = np.sum(mean26/(len(mean26)))
Nsum26 = mean26 - sliced_av26

asum26 = np.append(z1, Nsum26)
NewSum26 = np.append(asum26, z2)

ham26 = np.hamming(len(NewSum26))
win26 = ham26 * NewSum26

al = []
for i in np.arange(-15, 15):
    rwin26 = np.roll(win26, i)
    totalsum = np.sum(win * rwin26)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents26 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents26)
################################################################################################################
file27 = './Nov-17-2017/solar_scan-027-2.fit'
header27 = pf.getheader(file27) #Get the header from the fits file
data27 = pf.getdata(file27) #Get the CCD data from the fits file

sliced27 = data27[660:703,:]
od27 = np.sum(sliced27, axis = 0)
mean27 = od27/43

sliced_av27 = np.sum(mean27/(len(mean27)))
Nsum27 = mean27 - sliced_av27

asum27 = np.append(z1, Nsum27)
NewSum27 = np.append(asum27, z2)

ham27 = np.hamming(len(NewSum27))
win27 = ham27 * NewSum27

al = []
for i in np.arange(-15, 15):
    rwin27 = np.roll(win27, i)
    totalsum = np.sum(win * rwin27)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents27 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents27)
################################################################################################################
file28 = './Nov-17-2017/solar_scan-028-2.fit'
header28 = pf.getheader(file28) #Get the header from the fits file
data28 = pf.getdata(file28) #Get the CCD data from the fits file

sliced28 = data28[660:703,:]
od28 = np.sum(sliced28, axis = 0)
mean28 = od28/43

sliced_av28 = np.sum(mean28/(len(mean28)))
Nsum28 = mean28 - sliced_av28

asum28 = np.append(z1, Nsum28)
NewSum28 = np.append(asum28, z2)

ham28 = np.hamming(len(NewSum28))
win28 = ham28 * NewSum28

al = []
for i in np.arange(-15, 15):
    rwin28 = np.roll(win28, i)
    totalsum = np.sum(win * rwin28)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents28 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents28)
################################################################################################################
file29 = './Nov-17-2017/solar_scan-029-2.fit'
header29 = pf.getheader(file29) #Get the header from the fits file
data29 = pf.getdata(file29) #Get the CCD data from the fits file

sliced29 = data29[660:703,:]
od29 = np.sum(sliced29, axis = 0)
mean29 = od29/43

sliced_av29 = np.sum(mean29/(len(mean29)))
Nsum29 = mean29 - sliced_av29

asum29 = np.append(z1, Nsum29)
NewSum29 = np.append(asum29, z2)

ham29 = np.hamming(len(NewSum29))
win29 = ham29 * NewSum29

al = []
for i in np.arange(-15, 15):
    rwin29 = np.roll(win29, i)
    totalsum = np.sum(win * rwin29)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents29 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents29)
################################################################################################################
file30 = './Nov-17-2017/solar_scan-030-2.fit'
header30 = pf.getheader(file30) #Get the header from the fits file
data30 = pf.getdata(file30) #Get the CCD data from the fits file

sliced30 = data30[660:703,:]
od30 = np.sum(sliced30, axis = 0)
mean30 = od30/43

sliced_av30 = np.sum(mean30/(len(mean30)))
Nsum30 = mean30 - sliced_av30

asum30 = np.append(z1, Nsum30)
NewSum30 = np.append(asum30, z2)

ham30 = np.hamming(len(NewSum30))
win30 = ham30 * NewSum30

al = []
for i in np.arange(-15, 15):
    rwin30 = np.roll(win30, i)
    totalsum = np.sum(win * rwin30)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents30 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents30)
################################################################################################################
file31 = './Nov-17-2017/solar_scan-031-2.fit'
header31 = pf.getheader(file31) #Get the header from the fits file
data31 = pf.getdata(file31) #Get the CCD data from the fits file

sliced31 = data31[660:703,:]
od31 = np.sum(sliced31, axis = 0)
mean31 = od31/43

sliced_av31 = np.sum(mean31/(len(mean31)))
Nsum31 = mean31 - sliced_av31

asum31 = np.append(z1, Nsum31)
NewSum31 = np.append(asum31, z2)

ham31 = np.hamming(len(NewSum31))
win31 = ham31 * NewSum31

al = []
for i in np.arange(-15, 15):
    rwin31 = np.roll(win31, i)
    totalsum = np.sum(win * rwin31)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents31 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents31)
################################################################################################################
file32 = './Nov-17-2017/solar_scan-032-2.fit'
header32 = pf.getheader(file32) #Get the header from the fits file
data32 = pf.getdata(file32) #Get the CCD data from the fits file

sliced32 = data32[660:703,:]
od32 = np.sum(sliced32, axis = 0)
mean32 = od32/43

sliced_av32 = np.sum(mean32/(len(mean32)))
Nsum32 = mean32 - sliced_av32

asum32 = np.append(z1, Nsum32)
NewSum32 = np.append(asum32, z2)

ham32 = np.hamming(len(NewSum32))
win32 = ham32 * NewSum32

al = []
for i in np.arange(-15, 15):
    rwin32 = np.roll(win32, i)
    totalsum = np.sum(win * rwin32)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents32 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents32)
################################################################################################################
file33 = './Nov-17-2017/solar_scan-033-2.fit'
header33 = pf.getheader(file33) #Get the header from the fits file
data33 = pf.getdata(file33) #Get the CCD data from the fits file

sliced33 = data33[660:703,:]
od33 = np.sum(sliced33, axis = 0)
mean33 = od33/43

sliced_av33 = np.sum(mean33/(len(mean33)))
Nsum33 = mean33 - sliced_av33

asum33 = np.append(z1, Nsum33)
NewSum33 = np.append(asum33, z2)

ham33 = np.hamming(len(NewSum33))
win33 = ham33 * NewSum33

al = []
for i in np.arange(-15, 15):
    rwin33 = np.roll(win33, i)
    totalsum = np.sum(win * rwin33)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents33 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents33)
################################################################################################################
file34 = './Nov-17-2017/solar_scan-034-2.fit'
header34 = pf.getheader(file34) #Get the header from the fits file
data34 = pf.getdata(file34) #Get the CCD data from the fits file

sliced34 = data34[660:703,:]
od34 = np.sum(sliced34, axis = 0)
mean34 = od34/43

sliced_av34 = np.sum(mean34/(len(mean34)))
Nsum34 = mean34 - sliced_av34

asum34 = np.append(z1, Nsum34)
NewSum34 = np.append(asum34, z2)

ham34 = np.hamming(len(NewSum34))
win34 = ham34 * NewSum34

al = []
for i in np.arange(-15, 15):
    rwin34 = np.roll(win34, i)
    totalsum = np.sum(win * rwin34)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents34 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents34)
################################################################################################################
file35 = './Nov-17-2017/solar_scan-035-2.fit'
header35 = pf.getheader(file35) #Get the header from the fits file
data35 = pf.getdata(file35) #Get the CCD data from the fits file

sliced35 = data35[660:703,:]
od35 = np.sum(sliced35, axis = 0)
mean35 = od35/43

sliced_av35 = np.sum(mean35/(len(mean35)))
Nsum35 = mean35 - sliced_av35

asum35 = np.append(z1, Nsum35)
NewSum35 = np.append(asum35, z2)

ham35 = np.hamming(len(NewSum35))
win35 = ham35 * NewSum35

al = []
for i in np.arange(-15, 15):
    rwin35 = np.roll(win35, i)
    totalsum = np.sum(win * rwin35)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents35 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents35)
################################################################################################################
file36 = './Nov-17-2017/solar_scan-036-2.fit'
header36 = pf.getheader(file36) #Get the header from the fits file
data36 = pf.getdata(file36) #Get the CCD data from the fits file

sliced36 = data36[660:703,:]
od36 = np.sum(sliced36, axis = 0)
mean36 = od36/43

sliced_av36 = np.sum(mean36/(len(mean36)))
Nsum36 = mean36 - sliced_av36

asum36 = np.append(z1, Nsum36)
NewSum36 = np.append(asum36, z2)

ham36 = np.hamming(len(NewSum36))
win36 = ham36 * NewSum36

al = []
for i in np.arange(-15, 15):
    rwin36 = np.roll(win36, i)
    totalsum = np.sum(win * rwin36)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents36 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents36)
################################################################################################################
file37 = './Nov-17-2017/solar_scan-037-2.fit'
header37 = pf.getheader(file37) #Get the header from the fits file
data37 = pf.getdata(file37) #Get the CCD data from the fits file

sliced37 = data37[660:703,:]
od37 = np.sum(sliced37, axis = 0)
mean37 = od37/43

sliced_av37 = np.sum(mean37/(len(mean37)))
Nsum37 = mean37 - sliced_av37

asum37 = np.append(z1, Nsum37)
NewSum37 = np.append(asum37, z2)

ham37 = np.hamming(len(NewSum37))
win37 = ham37 * NewSum37

al = []
for i in np.arange(-15, 15):
    rwin37 = np.roll(win37, i)
    totalsum = np.sum(win * rwin37)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents37 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents37)
################################################################################################################
file38 = './Nov-17-2017/solar_scan-038-2.fit'
header38 = pf.getheader(file38) #Get the header from the fits file
data38 = pf.getdata(file38) #Get the CCD data from the fits file

sliced38 = data38[660:703,:]
od38 = np.sum(sliced38, axis = 0)
mean38 = od38/43

sliced_av38 = np.sum(mean38/(len(mean38)))
Nsum38 = mean38 - sliced_av38

asum38 = np.append(z1, Nsum38)
NewSum38 = np.append(asum38, z2)

ham38 = np.hamming(len(NewSum38))
win38 = ham38 * NewSum38

al = []
for i in np.arange(-15, 15):
    rwin38 = np.roll(win38, i)
    totalsum = np.sum(win * rwin38)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents38 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents38)
################################################################################################################
file39 = './Nov-17-2017/solar_scan-039-2.fit'
header39 = pf.getheader(file39) #Get the header from the fits file
data39 = pf.getdata(file39) #Get the CCD data from the fits file

sliced39 = data39[660:703,:]
od39 = np.sum(sliced39, axis = 0)
mean39 = od39/43

sliced_av39 = np.sum(mean39/(len(mean39)))
Nsum39 = mean39 - sliced_av39

asum39 = np.append(z1, Nsum39)
NewSum39 = np.append(asum39, z2)

ham39 = np.hamming(len(NewSum39))
win39 = ham39 * NewSum39

al = []
for i in np.arange(-15, 15):
    rwin39 = np.roll(win39, i)
    totalsum = np.sum(win * rwin39)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents39 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents39)
################################################################################################################
file40 = './Nov-17-2017/solar_scan-040-2.fit'
header40 = pf.getheader(file40) #Get the header from the fits file
data40 = pf.getdata(file40) #Get the CCD data from the fits file

sliced40 = data40[660:703,:]
od40 = np.sum(sliced40, axis = 0)
mean40 = od40/43

sliced_av40 = np.sum(mean40/(len(mean40)))
Nsum40 = mean40 - sliced_av40

asum40 = np.append(z1, Nsum40)
NewSum40 = np.append(asum40, z2)

ham40 = np.hamming(len(NewSum40))
win40 = ham40 * NewSum40

al = []
for i in np.arange(-15, 15):
    rwin40 = np.roll(win40, i)
    totalsum = np.sum(win * rwin40)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents40 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents40)
################################################################################################################
file41 = './Nov-17-2017/solar_scan-041-2.fit'
header41 = pf.getheader(file41) #Get the header from the fits file
data41 = pf.getdata(file41) #Get the CCD data from the fits file

sliced41 = data41[660:703,:]
od41 = np.sum(sliced41, axis = 0)
mean41 = od41/43

sliced_av41 = np.sum(mean41/(len(mean41)))
Nsum41 = mean41 - sliced_av41

asum41 = np.append(z1, Nsum41)
NewSum41 = np.append(asum41, z2)

ham41 = np.hamming(len(NewSum41))
win41 = ham41 * NewSum41

al = []
for i in np.arange(-15, 15):
    rwin41 = np.roll(win41, i)
    totalsum = np.sum(win * rwin41)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents41 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents41)
################################################################################################################
file42 = './Nov-17-2017/solar_scan-042-2.fit'
header42 = pf.getheader(file42) #Get the header from the fits file
data42 = pf.getdata(file42) #Get the CCD data from the fits file

sliced42 = data42[660:703,:]
od42 = np.sum(sliced42, axis = 0)
mean42 = od42/43

sliced_av42 = np.sum(mean42/(len(mean42)))
Nsum42 = mean42 - sliced_av42

asum42 = np.append(z1, Nsum42)
NewSum42 = np.append(asum42, z2)

ham42 = np.hamming(len(NewSum42))
win42 = ham42 * NewSum42

al = []
for i in np.arange(-15, 15):
    rwin42 = np.roll(win42, i)
    totalsum = np.sum(win * rwin42)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents42 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents42)
################################################################################################################
file43 = './Nov-17-2017/solar_scan-043-2.fit'
header43 = pf.getheader(file43) #Get the header from the fits file
data43 = pf.getdata(file43) #Get the CCD data from the fits file

sliced43 = data43[660:703,:]
od43 = np.sum(sliced43, axis = 0)
mean43 = od43/43

sliced_av43 = np.sum(mean43/(len(mean43)))
Nsum43 = mean43 - sliced_av43

asum43 = np.append(z1, Nsum43)
NewSum43 = np.append(asum43, z2)

ham43 = np.hamming(len(NewSum43))
win43 = ham43 * NewSum43

al = []
for i in np.arange(-15, 15):
    rwin43 = np.roll(win43, i)
    totalsum = np.sum(win * rwin43)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents43 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents43)
################################################################################################################
file44 = './Nov-17-2017/solar_scan-044-2.fit'
header44 = pf.getheader(file44) #Get the header from the fits file
data44 = pf.getdata(file44) #Get the CCD data from the fits file

sliced44 = data44[660:703,:]
od44 = np.sum(sliced44, axis = 0)
mean44 = od44/43

sliced_av44 = np.sum(mean44/(len(mean44)))
Nsum44 = mean44 - sliced_av44

asum44 = np.append(z1, Nsum44)
NewSum44 = np.append(asum44, z2)

ham44 = np.hamming(len(NewSum44))
win44 = ham44 * NewSum44

al = []
for i in np.arange(-15, 15):
    rwin44 = np.roll(win44, i)
    totalsum = np.sum(win * rwin44)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents44 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents44)
################################################################################################################
file45 = './Nov-17-2017/solar_scan-045-2.fit'
header45 = pf.getheader(file45) #Get the header from the fits file
data45 = pf.getdata(file45) #Get the CCD data from the fits file

sliced45 = data45[660:703,:]
od45 = np.sum(sliced45, axis = 0)
mean45 = od45/43

sliced_av45 = np.sum(mean45/(len(mean45)))
Nsum45 = mean45 - sliced_av45

asum45 = np.append(z1, Nsum45)
NewSum45 = np.append(asum45, z2)

ham45 = np.hamming(len(NewSum45))
win45 = ham45 * NewSum45

al = []
for i in np.arange(-15, 15):
    rwin45 = np.roll(win45, i)
    totalsum = np.sum(win * rwin45)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents45 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents45)
################################################################################################################
file46 = './Nov-17-2017/solar_scan-046-2.fit'
header46 = pf.getheader(file46) #Get the header from the fits file
data46 = pf.getdata(file46) #Get the CCD data from the fits file

sliced46 = data46[660:703,:]
od46 = np.sum(sliced46, axis = 0)
mean46 = od46/43

sliced_av46 = np.sum(mean46/(len(mean46)))
Nsum46 = mean46 - sliced_av46

asum46 = np.append(z1, Nsum46)
NewSum46 = np.append(asum46, z2)

ham46 = np.hamming(len(NewSum46))
win46 = ham46 * NewSum46

al = []
for i in np.arange(-15, 15):
    rwin46 = np.roll(win46, i)
    totalsum = np.sum(win * rwin46)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents46 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents46)
################################################################################################################
file47 = './Nov-17-2017/solar_scan-047-2.fit'
header47 = pf.getheader(file47) #Get the header from the fits file
data47 = pf.getdata(file47) #Get the CCD data from the fits file

sliced47 = data47[660:703,:]
od47 = np.sum(sliced47, axis = 0)
mean47 = od47/43

sliced_av47 = np.sum(mean47/(len(mean47)))
Nsum47 = mean47 - sliced_av47

asum47 = np.append(z1, Nsum47)
NewSum47 = np.append(asum47, z2)

ham47 = np.hamming(len(NewSum47))
win47 = ham47 * NewSum47

al = []
for i in np.arange(-15, 15):
    rwin47 = np.roll(win47, i)
    totalsum = np.sum(win * rwin47)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents47 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents47)
################################################################################################################
file48 = './Nov-17-2017/solar_scan-048-2.fit'
header48 = pf.getheader(file48) #Get the header from the fits file
data48 = pf.getdata(file48) #Get the CCD data from the fits file

sliced48 = data48[660:703,:]
od48 = np.sum(sliced48, axis = 0)
mean48 = od48/43

sliced_av48 = np.sum(mean48/(len(mean48)))
Nsum48 = mean48 - sliced_av48

asum48 = np.append(z1, Nsum48)
NewSum48 = np.append(asum48, z2)

ham48 = np.hamming(len(NewSum48))
win48 = ham48 * NewSum48

al = []
for i in np.arange(-15, 15):
    rwin48 = np.roll(win48, i)
    totalsum = np.sum(win * rwin48)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents48 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents48)
################################################################################################################
file49 = './Nov-17-2017/solar_scan-049-2.fit'
header49 = pf.getheader(file49) #Get the header from the fits file
data49 = pf.getdata(file49) #Get the CCD data from the fits file

sliced49 = data49[660:703,:]
od49 = np.sum(sliced49, axis = 0)
mean49 = od49/43

sliced_av49 = np.sum(mean49/(len(mean49)))
Nsum49 = mean49 - sliced_av49

asum49 = np.append(z1, Nsum49)
NewSum49 = np.append(asum49, z2)

ham49 = np.hamming(len(NewSum49))
win49 = ham49 * NewSum49

al = []
for i in np.arange(-15, 15):
    rwin49 = np.roll(win49, i)
    totalsum = np.sum(win * rwin49)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents49 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents49)
################################################################################################################
file50 = './Nov-17-2017/solar_scan-050-2.fit'
header50 = pf.getheader(file50) #Get the header from the fits file
data50 = pf.getdata(file50) #Get the CCD data from the fits file

sliced50 = data50[660:703,:]
od50 = np.sum(sliced50, axis = 0)
mean50 = od50/43

sliced_av50 = np.sum(mean50/(len(mean50)))
Nsum50 = mean50 - sliced_av50

asum50 = np.append(z1, Nsum50)
NewSum50 = np.append(asum50, z2)

ham50 = np.hamming(len(NewSum50))
win50 = ham50 * NewSum50

al = []
for i in np.arange(-15, 15):
    rwin50 = np.roll(win50, i)
    totalsum = np.sum(win * rwin50)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents50 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents50)
################################################################################################################
file51 = './Nov-17-2017/solar_scan-051-2.fit'
header51 = pf.getheader(file51) #Get the header from the fits file
data51 = pf.getdata(file51) #Get the CCD data from the fits file

sliced51 = data51[660:703,:]
od51 = np.sum(sliced51, axis = 0)
mean51 = od51/43

sliced_av51 = np.sum(mean51/(len(mean51)))
Nsum51 = mean51 - sliced_av51

asum51 = np.append(z1, Nsum51)
NewSum51 = np.append(asum51, z2)

ham51 = np.hamming(len(NewSum51))
win51 = ham51 * NewSum51

al = []
for i in np.arange(-15, 15):
    rwin51 = np.roll(win51, i)
    totalsum = np.sum(win * rwin51)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents51 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents51)
################################################################################################################
file52 = './Nov-17-2017/solar_scan-052-2.fit'
header52 = pf.getheader(file52) #Get the header from the fits file
data52 = pf.getdata(file52) #Get the CCD data from the fits file

sliced52 = data52[660:703,:]
od52 = np.sum(sliced52, axis = 0)
mean52 = od52/43

sliced_av52 = np.sum(mean52/(len(mean52)))
Nsum52 = mean52 - sliced_av52

asum52 = np.append(z1, Nsum52)
NewSum52 = np.append(asum52, z2)

ham52 = np.hamming(len(NewSum52))
win52 = ham52 * NewSum52

al = []
for i in np.arange(-15, 15):
    rwin52 = np.roll(win52, i)
    totalsum = np.sum(win * rwin52)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents52 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents52)
################################################################################################################
file53 = './Nov-17-2017/solar_scan-053-2.fit'
header53 = pf.getheader(file53) #Get the header from the fits file
data53 = pf.getdata(file53) #Get the CCD data from the fits file

sliced53 = data53[660:703,:]
od53 = np.sum(sliced53, axis = 0)
mean53 = od53/43

sliced_av53 = np.sum(mean53/(len(mean53)))
Nsum53 = mean53 - sliced_av53

asum53 = np.append(z1, Nsum53)
NewSum53 = np.append(asum53, z2)

ham53 = np.hamming(len(NewSum53))
win53 = ham53 * NewSum53

al = []
for i in np.arange(-15, 15):
    rwin53 = np.roll(win53, i)
    totalsum = np.sum(win * rwin53)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents53 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents53)
################################################################################################################
file54 = './Nov-17-2017/solar_scan-054-2.fit'
header54 = pf.getheader(file54) #Get the header from the fits file
data54 = pf.getdata(file54) #Get the CCD data from the fits file

sliced54 = data54[660:703,:]
od54 = np.sum(sliced54, axis = 0)
mean54 = od54/43

sliced_av54 = np.sum(mean54/(len(mean54)))
Nsum54 = mean54 - sliced_av54

asum54 = np.append(z1, Nsum54)
NewSum54 = np.append(asum54, z2)

ham54 = np.hamming(len(NewSum54))
win54 = ham54 * NewSum54

al = []
for i in np.arange(-15, 15):
    rwin54 = np.roll(win54, i)
    totalsum = np.sum(win * rwin54)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents54 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents54)
################################################################################################################
file55 = './Nov-17-2017/solar_scan-055-2.fit'
header55 = pf.getheader(file55) #Get the header from the fits file
data55 = pf.getdata(file55) #Get the CCD data from the fits file

sliced55 = data55[660:703,:]
od55 = np.sum(sliced55, axis = 0)
mean55 = od55/43

sliced_av55 = np.sum(mean55/(len(mean55)))
Nsum55 = mean55 - sliced_av55

asum55 = np.append(z1, Nsum55)
NewSum55 = np.append(asum55, z2)

ham55 = np.hamming(len(NewSum55))
win55 = ham55 * NewSum55

al = []
for i in np.arange(-15, 15):
    rwin55 = np.roll(win55, i)
    totalsum = np.sum(win * rwin55)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents55 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents55)
################################################################################################################
file56 = './Nov-17-2017/solar_scan-056-2.fit'
header56 = pf.getheader(file56) #Get the header from the fits file
data56 = pf.getdata(file56) #Get the CCD data from the fits file

sliced56 = data56[660:703,:]
od56 = np.sum(sliced56, axis = 0)
mean56 = od56/43

sliced_av56 = np.sum(mean56/(len(mean56)))
Nsum56 = mean56 - sliced_av56

asum56 = np.append(z1, Nsum56)
NewSum56 = np.append(asum56, z2)

ham56 = np.hamming(len(NewSum56))
win56 = ham56 * NewSum56

al = []
for i in np.arange(-15, 15):
    rwin56 = np.roll(win56, i)
    totalsum = np.sum(win * rwin56)
    al.append(totalsum)

zal = np.append(z1, al)
nal = np.append(zal, z2)

cents56 = find_all_centroids(np.arange(0, len(nal)), nal)
alcents.append(cents56)

xcents = np.arange(0, len(alcents))

plt.figure()
plt.plot(xcents, alcents, 'b.')
plt.show()

Nfit, covN = np.polyfit(xcents, alcents, 1, full=False, cov = True)
y = Nfit[0]*xcents + Nfit[1]
plt.plot(xcents, alcents, 'o')
plt.plot(xcents, y, 'g-')
plt.title('Correlation Centroid Fit', fontsize=17)
plt.ylabel('Pixels', fontsize=17)
plt.xlabel('Time', fontsize=17)
plt.show()

print('Coeff = ', Nfit)
print('y-intercept = ', Nfit[0])
print('Slope = ', Nfit[1])

print('alcents = ', alcents)