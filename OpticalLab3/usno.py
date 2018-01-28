#######################################################################################################################
def usno(radeg,decdeg,fovam,epoch):
	import string as str
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
		kw =k.split('\t')
		ded0 = float(kw[2])
		pmrad = float(kw[6])/3600e3/np.cos(np.deg2rad(ded0)) #convert from mas/yr to deg/year
		pmded = float(kw[7])/3600e3

		name = np.append(name,kw[0])
		rad = np.append(rad, float(kw[1]) + pmrad*(epoch-2000.0))
		ded = np.append(ded,float(kw[2]) + pmded*(epoch-2000.0))

		if kw[12] != '  ':
			rmag = np.append(rmag,float(kw[12]))
		else:
			rmag = np.append(rmag,np.nan)
	return name,rad,ded,rmag




#######################################################################################################################
ras = '20:59:54.03'
des = '+07:16:50.6'

radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
dsgn = np.sign(float(des[0:3]))
dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.

fovam = 22.0
name, rad, ded, rmag = usno(radeg, dedeg, fovam, 2000)

w = np.where(rmag < 12.)[0]

plt.plot(rad[w],ded[w],'g.')
plt.locator_params(axis='x',nbins=4)
plt.locator_params(axis='y',nbins=4)
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]')
plt.ticklabel_format(useOffset=false)
plt.axis('scald')

ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])