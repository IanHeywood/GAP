import Pyxis
import pylab
import Tigger
import numpy
import os
import sys
import glob
import pyfits
from astLib import astWCS
from matplotlib.patches import Ellipse


myms = sys.argv[1]


def square_6x6():
        # stolen from
        # https://svn.atnf.csiro.au/askap/ACES/pythonlib/aces/footprint_class.py
        offsets = numpy.array([[-0.5,+0.5],[+0.5,+0.5],[-0.5,-0.5],[+0.5,-0.5],[-1.5,+1.5],[-0.5,+1.5],
        [+0.5,+1.5],[+1.5,+1.5],[+1.5,+0.5],[+1.5,-0.5],[+1.5,-1.5],[+0.5,-1.5],
        [-0.5,-1.5],[-1.5,-1.5],[-1.5,-0.5],[-1.5,+0.5],[-2.5,+2.5],[-1.5,+2.5],
        [-0.5,+2.5],[+0.5,+2.5],[+1.5,+2.5],[+2.5,+2.5],[+2.5,+1.5],[+2.5,+0.5],
        [+2.5,-0.5],[+2.5,-1.5],[+2.5,-2.5],[+1.5,-2.5],[+0.5,-2.5],[-0.5,-2.5],
        [-1.5,-2.5],[-2.5,-2.5],[-2.5,-1.5],[-2.5,-0.5],[-2.5,+0.5],[-2.5,+1.5]])
        return offsets


def dEtags(inputLsm,namelist):
	print '     Applying tags'
	cluster_leads = []
	skymodel = Tigger.load(inputLsm,verbose=False)
	for src in skymodel.sources:
		name = src.name
		cluster_size = src.getTag('cluster_size')
		if name in namelist:
			print '    ',name,'tagged'
			src.setAttribute('dE',True)
			if cluster_size > 1:
				print '    ',name,'is part of a cluster of size',cluster_size
				cluster_leads.append(src.getTag('cluster'))
	for src in skymodel.sources:
		cluster = src.getTag('cluster')
		if cluster in cluster_leads:
			print '          Tagging',src.name,'as part of cluster',cluster
			src.setAttribute('dE',True)
	Tigger.save(skymodel,inputLsm,verbose=False)


def getBeam(infits):
        part = infits.split('.')[1].split('_')[0]
        if part[:4] == 'beam':
                beam = part.replace('beam','')
        return int(beam)


def srcthumbs(inputFits,inputLsm):
	padding = 30
	srclist = []

	input_hdu = pyfits.open(inputFits)[0]
	hdr = input_hdu.header
	WCS = astWCS.WCS(hdr, mode = 'pyfits')
	if len(input_hdu.data.shape) == 2:
		image = numpy.array(input_hdu.data[:,:])
	elif len(input_hdu.data.shape) == 3:
		image = numpy.array(input_hdu.data[0,:,:])
	else:
		image = numpy.array(input_hdu.data[0,0,:,:])

	pitch = 0.9
	offsets = pitch*square_6x6()
	beam = getBeam(inputFits)
	offset = offsets[beam]

	ra_ptg = hdr.get('CRVAL1')
	dec_ptg = hdr.get('CRVAL2')

	ra0 = ra_ptg + offset[0]
	dec0 = dec_ptg + offset[1]

	print ra_ptg,dec_ptg,ra0,dec0,beam

	tmp_x = []
	tmp_y = []
	for src in Tigger.load(inputLsm, verbose = False).sources:
		tmp_x.append(src.pos.ra*180.0/numpy.pi)
		tmp_y.append(src.pos.dec*180.0/numpy.pi)
	mx = numpy.mean(numpy.array(tmp_x)) 
	my = numpy.mean(numpy.array(tmp_y)) 
	print 'mx,my =',mx,my

	for src in Tigger.load(inputLsm, verbose = False).sources:
		name = src.name
		flux = src.flux.I
		ra_d = src.pos.ra*180.0/numpy.pi
		dec_d = src.pos.dec*180.0/numpy.pi

		ra_pix,dec_pix = WCS.wcs2pix(ra_d,dec_d)
		x0,x1 = int(ra_pix-padding),int(1+ra_pix+padding)
		y0,y1 = int(dec_pix-padding),int(1+dec_pix+padding)
		if x0 < 1: 
			x0 = 0
		if y0 < 1:
			y0 = 0
		if x1 > image.shape[1]-1:
			x1 = image.shape[1]
		if y1 > image.shape[0]-1:
			y1 = image.shape[0]
		thumbnail = image[y0:y1,x0:x1]
		rms = numpy.std(thumbnail)

		# This probably works for anything other than ASKAP
		# r = src.r

		dx = ra_d - ra0
		dy = dec_d - dec0
		r = ((dx**2.0)+(dy**2.0))**0.5
		print r

		srclist.append((name,ra_d,dec_d,flux,rms,r))
	return srclist


def thresh(x0,y0,xt,yt,x,y):
	locus = (((x-x0)**2.0)/(xt**2.0))+(((y-y0)**2.0)/(yt**2.0))
	if locus > 0.9:
		return True
	else:
		return False


pitch = 0.9
offsets = pitch*square_6x6()
lsmlist = glob.glob('img_*'+myms+'*.html')

#fig = pylab.figure(figsize=(32,24))

count = 1

for mylsm in lsmlist:
	beam = mylsm.split('_')[3]
	print '--------------------------------------'
	print 'Processing ',beam

#	ax = fig.add_subplot(6,6,count)
	xx = []
	yy = []
	names = []
	namelist = []
	prefix = mylsm.split('.ms')[0]
	resid_fits = glob.glob(prefix+'*MFS-residual.fits')[0]
	srclist = srcthumbs(resid_fits,mylsm)

	# srclist.append((name,ra_d,dec_d,flux,rms,r))

	for src in srclist:
		rad_wt_flux = src[3]*(src[5]**0.5)
		xx.append(rad_wt_flux)
		yy.append(src[4])
		names.append(src[0])

	xx = numpy.array(xx)
	yy = numpy.array(yy)

	xx_std = numpy.std(xx)
	yy_std = numpy.std(yy)
	xx_m = numpy.median(xx)
	yy_m = numpy.median(yy)

	for n in [2,4,6,8]:
		ell = Ellipse(xy=(xx_m,yy_m),width=n*xx_std,height=n*yy_std)
		ax.add_artist(ell)
		ell.set_facecolor('none')
		ell.set_edgecolor('black')
		ell.set_alpha(0.8)
		ell.set_linestyle('dashed')

	for i in range(0,len(xx)):
#		if xx[i] > 3.5*xx_std and yy[i] > 3.0*yy_std:
		if thresh(xx_m,yy_m,3.0*xx_std,2.5*yy_std,xx[i],yy[i]):
			ax.plot(xx[i],yy[i],'.',markersize=10,alpha=0.5,color='red')
			print '     Identified '+names[i]+' for dE'
			namelist.append(names[i])
		else:
			ax.plot(xx[i],yy[i],'.',markersize=10,alpha=0.5,color='blue')

	dEtags(mylsm,namelist)

	# beam_median = numpy.median(xx)
	# frac_xx = xx / beam_median
	# frac_std = numpy.std(frac_xx)
	# flux_std = numpy.std(yy)
	# for i in range(0,len(xx)):
	# 	if frac_xx[i] > 4*frac_std and yy[i]/xx[i] > 5:
	# 		print '      Tag dE: '+names[i]
	# 		ax.plot(frac_xx[i],yy[i],'.',markersize=10,alpha=0.5,color='red')
	# 	else:
	# 		ax.plot(frac_xx[i],yy[i],'.',markersize=10,alpha=0.5,color='blue')
	count += 1
#	ax.set_title(beam)
	# ax.set_xscale('log')
	# ax.set_yscale('log')
	# ax.set_ylim((1e-4,10))
	# ax.set_xlim((0.1,10))
#	ax.set_xlabel('Rad. weighted total intensity')
#	ax.set_ylabel('Local RMS')

#fig.tight_layout()
#fig.savefig('localrms.png')
