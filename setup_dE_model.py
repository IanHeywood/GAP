import pickle
import numpy
import os
import Pyxis
import Tigger
import sys
import glob
import pyfits
import string
import random
from astLib import astWCS
from astLib import astCoords as ac
import bdsf

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Function definitions
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


def gi(message):
        print '\033[92m'+message+'\033[0m'


def ri(message):
        print '\033[91m'+message+'\033[0m'


def rad2deg(x):
	return 180.0*x/numpy.pi


def tempname(size=12,chars=string.ascii_uppercase+string.digits+string.ascii_lowercase):
	return ''.join(random.choice(chars) for _ in range(size))


def getfreq(f0,f1,chan,nchan):
	chan = int(chan)
	df = (float(f1) - float(f0))/float(nchan)
	return f0+(chan*df)+(df/2.0)


def tiggerConvert(gaul):
	args = []
	tigger_convert  = Pyxis.x("tigger-convert")
	#Dictionary for establishing correspondence between parameter names in gaul files produced by pybdsm, and pyxis parameter names
	dict_gaul2lsm = {'Gaus_id':'name', 'Isl_id':'Isl_id', 'Source_id':'Source_id', 'Wave_id':'Wave_id', 'RA':'ra_d', 'E_RA':'E_RA', 'DEC':'dec_d', 'E_DEC':'E_DEC', 'Total_flux':'i', 'E_Total_flux':'E_Total_flux', 'Peak_flux':'Peak_flux', 'E_Peak_flux':'E_Peak_flux', 'Xposn':'Xposn', 'E_Xposn':'E_Xposn', 'Yposn':'Yposn', 'E_Yposn':'E_Yposn', 'Maj':'Maj', 'E_Maj':'E_Maj', 'Min':'Min', 'E_Min':'E_Min', 'PA':'PA', 'E_PA':'E_PA', 'Maj_img_plane':'Maj_img_plane', 'E_Maj_img_plane':'E_Maj_img_plane', 'Min_img_plane':'Min_img_plane', 'E_Min_img_plane':'E_Min_img_plane', 'PA_img_plane':'PA_img_plane', 'E_PA_img_plane':'E_PA_img_plane', 'DC_Maj':'emaj_d', 'E_DC_Maj':'E_DC_Maj', 'DC_Min':'emin_d', 'E_DC_Min':'E_DC_Min', 'DC_PA':'pa_d', 'E_DC_PA':'E_DC_PA', 'DC_Maj_img_plane':'DC_Maj_img_plane', 'E_DC_Maj_img_plane':'E_DC_Maj_img_plane', 'DC_Min_img_plane':'DC_Min_img_plane', 'E_DC_Min_img_plane':'E_DC_Min_img_plane', 'DC_PA_img_plane':'DC_PA_img_plane', 'E_DC_PA_img_plane':'E_DC_PA_img_plane', 'Isl_Total_flux':'Isl_Total_flux', 'E_Isl_Total_flux':'E_Isl_Total_flux', 'Isl_rms':'Isl_rms', 'Isl_mean':'Isl_mean', 'Resid_Isl_rms':'Resid_Isl_rms', 'Resid_Isl_mean':'Resid_Isl_mean', 'S_Code':'S_Code', 'Total_Q':'q', 'E_Total_Q':'E_Total_Q', 'Total_U':'u', 'E_Total_U':'E_Total_U', 'Total_V':'v', 'E_Total_V':'E_Total_V', 'Linear_Pol_frac':'Linear_Pol_frac', 'Elow_Linear_Pol_frac':'Elow_Linear_Pol_frac', 'Ehigh_Linear_Pol_frac':'Ehigh_Linear_Pol_frac', 'Circ_Pol_Frac':'Circ_Pol_Frac', 'Elow_Circ_Pol_Frac':'Elow_Circ_Pol_Frac', 'Ehigh_Circ_Pol_Frac':'Ehigh_Circ_Pol_Frac', 'Total_Pol_Frac':'Total_Pol_Frac', 'Elow_Total_Pol_Frac':'Elow_Total_Pol_Frac', 'Ehigh_Total_Pol_Frac':'Ehigh_Total_Pol_Frac', 'Linear_Pol_Ang':'Linear_Pol_Ang', 'E_Linear_Pol_Ang':'E_Linear_Pol_Ang'}

	#Dictionary for classifying a parameter as a general parameter or a polarization-specific parameter
	dict_pol_flag = {'Gaus_id':0, 'Isl_id':0, 'Source_id':0, 'Wave_id':0, 'RA':0, 'E_RA':0, 'DEC':0, 'E_DEC':0, 'Total_flux':0, 'E_Total_flux':0, 'Peak_flux':0, 'E_Peak_flux':0, 'Xposn':0, 'E_Xposn':0, 'Yposn':0, 'E_Yposn':0, 'Maj':0, 'E_Maj':0, 'Min':0, 'E_Min':0, 'PA':0, 'E_PA':0, 'Maj_img_plane':0, 'E_Maj_img_plane':0, 'Min_img_plane':0, 'E_Min_img_plane':0, 'PA_img_plane':0, 'E_PA_img_plane':0, 'DC_Maj':0, 'E_DC_Maj':0, 'DC_Min':0, 'E_DC_Min':0, 'DC_PA':0, 'E_DC_PA':0, 'DC_Maj_img_plane':0, 'E_DC_Maj_img_plane':0, 'DC_Min_img_plane':0, 'E_DC_Min_img_plane':0, 'DC_PA_img_plane':0, 'E_DC_PA_img_plane':0, 'Isl_Total_flux':0, 'E_Isl_Total_flux':0, 'Isl_rms':0, 'Isl_mean':0, 'Resid_Isl_rms':0, 'Resid_Isl_mean':0, 'S_Code':0, 'Total_Q':1, 'E_Total_Q':1, 'Total_U':1, 'E_Total_U':1, 'Total_V':1, 'E_Total_V':1, 'Linear_Pol_frac':1, 'Elow_Linear_Pol_frac':1, 'Ehigh_Linear_Pol_frac':1, 'Circ_Pol_Frac':1, 'Elow_Circ_Pol_Frac':1, 'Ehigh_Circ_Pol_Frac':1, 'Total_Pol_Frac':1, 'Elow_Total_Pol_Frac':1, 'Ehigh_Total_Pol_Frac':1, 'Linear_Pol_Ang':1, 'E_Linear_Pol_Ang':1}

	lines = [line.strip() for line in open(gaul)]

	for line in range(len(lines)):
		if lines[line]:
			if lines[line].split()[0] is not '#': 
				gaul_params = lines[line-1].split()[1:] #Parameter list is last line in gaul file that begins with a '#'
				break

	# Initialize lists for general and polarization parameters 
	lsm_params_general = []
	lsm_params_polarization = []

	for param in gaul_params:
		if dict_pol_flag[param] is 0:
			lsm_params_general.append(dict_gaul2lsm[param])
		if dict_pol_flag[param] is 1:
			lsm_params_polarization.append(dict_gaul2lsm[param])

	general_params_string = ' '.join(lsm_params_general)
	pol_params_string = ' '.join(lsm_params_polarization)

	output = gaul.replace('.gaul','.lsm.html')

	cluster = 120.0
	MIN_EXTENT = 10.0

	tigger_convert(gaul,output,"-t","ASCII","--format", general_params_string,
		"-f","--rename",
		"--cluster-dist",cluster,
		"--min-extent",MIN_EXTENT,
		split_args=False,
		*args);
	return output


def makesubim(infits,ra,dec,size,outfits):
	# ra,dec,size all in degrees
	syscall = 'mSubimage '+infits+' '+outfits+' '
	syscall+= str(ra)+' '+str(dec)+' '+str(size)+' '
	gi('Writing subim: '+outfits)
	os.system(syscall)
	return outfits


def mosaic_directions(fitslist,outfits):
	gi('Making mosaic of subimages')
	tempdir = tempname()
	gi('Creating '+tempdir+'/repro')
	os.mkdir(tempdir)
	os.mkdir(tempdir+'/repro')
	os.chdir(tempdir)
	for item in fitslist:
		os.symlink('../'+item,item)
	os.system('mImgtbl . images.tbl')
	os.system('mMakeHdr images.tbl template.hdr')
	for item in fitslist:
		os.system('mProject '+item+' repro/'+item+' template.hdr')
	os.chdir('repro')
	os.system('mImgtbl . images.tbl')
	os.system('mAdd images.tbl ../template.hdr '+outfits)
	os.rename(outfits,'../../'+outfits)
	os.chdir('../../')
	if cleanup:
		gi('Removing '+tempdir)
		os.system('rm -rf '+tempdir)
	gi('Written '+outfits)
	return outfits


def makeLSM(infits):#,mybeam,myfreq):
	img = bdsf.process_image(infits,thresh_pix=9.0,thresh_isl=8.0)#,beam=mybeam,frequency=myfreq)
	foundsrcs = img.write_catalog(format='ascii',catalog_type='gaul',clobber=True,incl_empty=True)
	gaul = infits.replace('.fits','.pybdsm.gaul')
	oplsm = gaul.replace('gaul','lsm.html')
	if foundsrcs:
		tiggerConvert(gaul)
		if cleanup:
			os.system('rm '+gaul)
			os.system('rm '+gaul.replace('.pybdsm.gaul','.fits.pybdsm.log'))
		gi('Wrote LSM:     '+oplsm)
		return (oplsm,True)
	else:
		ri('No sources found in '+str(infits))
		ri('Could be corrupt image, could be trouble source entering null')
		return (oplsm,False)

	
def fixMontageHeaders(infile,outfile,axes):
	# Images produced by Montage do not have FREQ or STOKES axes
	# or information about the restoring beam. This confuses things like PyBDSM
	# infile provides the keywords to be written to outfile
	inphdu = pyfits.open(infile)
	inphdr = inphdu[0].header
	outhdu = pyfits.open(outfile,mode='update')
	outhdr = outhdu[0].header
	keywords = ['CTYPE','CRVAL','CDELT','CRPIX']
	for axis in axes: 
		for key in keywords:
			inkey = key+str(axis)
			outkey = key+str(axis)
			afterkey = key+str(axis-1)
			xx = inphdr[inkey]
			outhdr.set(outkey,xx,after=afterkey)
	outhdr.set('BUNIT',inphdr['BUNIT'],after=outkey)
	outhdr.set('BMAJ',inphdr['BMAJ'],after='BUNIT')
	outhdr.set('BMIN',inphdr['BMIN'],after='BMAJ')
	outhdr.set('BPA',inphdr['BPA'],after='BMIN')
	outhdu.flush()	
	

def tagsources(inlsm):
	model = Tigger.load(inlsm)
	gi('Reading LSM:   '+inlsm)
	srcs = model.sources
	counter = 0
	for src in srcs:
		src.setTag('dE',True)
		counter += 1
	gi('Tagged:        '+str(counter)+' source(s)')
	model.save(inlsm)


def add_dummy(inlsm):
	dummylsm = 'dummy.lsm.html'
	gi('Adding dummy:  '+inlsm) 
	dsrc = Tigger.load(dummylsm).sources
	model = Tigger.load(inlsm)
	model.sources.append(dsrc[0])
	model.save(inlsm)


def getImage(fitsfile):
	input_hdu = pyfits.open(fitsfile)[0]
	if len(input_hdu.data.shape) == 2:
	        image = numpy.array(input_hdu.data[:,:])
	elif len(input_hdu.data.shape) == 3:
	        image = numpy.array(input_hdu.data[0,:,:])
	else:
	        image = numpy.array(input_hdu.data[0,0,:,:])
	return image


def flushFits(newimage,fitsfile):
	gi('Flushing:      '+fitsfile)
	f = pyfits.open(fitsfile,mode='update')
	input_hdu = f[0]
	if len(input_hdu.data.shape) == 2:
	        input_hdu.data[:,:] = newimage
	elif len(input_hdu.data.shape) == 3:
	        input_hdu.data[0,:,:] = newimage
	else:
	        input_hdu.data[0,0,:,:] = newimage
	f.flush()


def maskimage(infits,directions,size):
	backup = True
	if backup:
		backupfits = infits+'.backup'
		if not os.path.isfile(backupfits):
			gi('Making backup: '+backupfits)
			os.system('cp '+infits+' '+backupfits)
		else:
			ri(backupfits+' already exists, will not overwrite')
	gi('Reading:       '+infits)
	input_hdu = pyfits.open(infits)[0]
	hdr = input_hdu.header
	WCS = astWCS.WCS(hdr,mode='pyfits')
	deg2pix = 1.0/hdr.get('CDELT2') # declination increment
	if len(input_hdu.data.shape) == 2:
		image = numpy.array(input_hdu.data[:,:])
	elif len(input_hdu.data.shape) == 3:
		image = numpy.array(input_hdu.data[0,:,:])
	else:
		image = numpy.array(input_hdu.data[0,0,:,:])
	dx = dy = extent*deg2pix/2.0
	for ra,dec in directions:
		gi('Direction:     '+str(ra)+' '+str(dec))
		xpix,ypix = WCS.wcs2pix(ra,dec)
		x0 = xpix - dx
		x1 = xpix + dx
		y0 = ypix - dy
		y1 = ypix + dy
		x0 = int(x0)
		x1 = int(x1)
		y0 = int(y0)
		y1 = int(y1)
		gi('Masking:       '+str(x0)+' '+str(x1)+' '+str(y0)+' '+str(y1))
		image[y0:y1,x0:x1] = 0.0
	flushFits(image,infits)


def predict(msname,imgbase):
	syscall = 'wsclean -predict -channelsout 4 -size 8192 8192 '
	syscall+= '-scale 3.5asec -name '+imgbase+' -mem 90 '
	syscall+= '-predict-channels 64 '+msname
	os.system(syscall)



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Switches
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# img_scienceData_SB4191_G23_T0-0B.beam03_averaged_fg_wtspec.ms_apcal-MFS-image.pybdsm.lsm.html_collapsed.txt

myms = sys.argv[1]
inp = glob.glob('*'+myms+'*.txt')[0]

mslist = [inp.split('.ms')[0].lstrip('img_')+'.ms']
imgbase = inp.split('-MFS')[0]

directions = []

f = open(inp)
line = f.readline()
while line:
	cols = line.split()
	dir = (float(cols[1]),float(cols[2]))
	directions.append(dir)
	line = f.readline()
f.close()

print mslist
print imgbase
print directions

extent = 0.06 # size of subimages and mask

cleanup = True
makelsms = True
maskmodels = True
runwspredict = True

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Guts
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

chan_images = sorted(glob.glob(imgbase+'-000*-image.fits'))
chans = []
for chan_image in chan_images:
	chan = chan_image.split('-')[2]
	print chan
	if chan not in chans:
		chans.append(chan)


#img_laduma_01_2048_wtspec_J033230-280757.ms_corr_pcal_fitsmask-0000-image.fits
if makelsms:
	full_merge_list = []
	for chan_image in chan_images:
		dir_merge_list = []
		for ra,dec in directions:
			# Make thumbnail images around troublemaker
			subim = chan_image.replace('.fits','')+'_'+str(ra)+'_'+str(dec)+'.fits'
			subfits = makesubim(chan_image,ra,dec,extent,subim)
			dir_merge_list.append(subfits)
		full_merge_list.append(dir_merge_list)

	# Merge images
	for dir_merge_list in full_merge_list:
		chan = chans[full_merge_list.index(dir_merge_list)]
		chan_dE_image = imgbase+'_'+chan+'_dE-sources.fits'
		if len(directions) > 1:
			mosaic_directions(dir_merge_list,chan_dE_image)
		else:
			os.rename(dir_merge_list[0],chan_dE_image)
		if cleanup:
			for tempfile in dir_merge_list:
				gi('Removing '+tempfile)
				os.system('rm '+tempfile)

	# Fix Montage headers
	for chan_image in chan_images:
		dE_img = chan_image.replace('-00','_00').replace('-image','_dE-sources')
		fixMontageHeaders(chan_image,dE_img,[3,4])
		
		
	# Find sources with a high threshold, tag LSMs with 'dE'
	for chan in chans:
		chan_dE_image = imgbase+'_'+chan+'_dE-sources.fits'
		# getfreq not needed anymore here because beam / freq info carried over from input
		# images with the fixMontageHeaders step above
		mylsm = makeLSM(chan_dE_image) #,(0.00125,0.00125,0.0),getfreq(1.0e9,2.0e9,chan,4))
		if mylsm[1]:
			tagsources(mylsm[0])
		#	add_dummy(mylsm[0])
		else:
			ri('Writing dummy lsm to '+mylsm[0])
			os.system('cp dummy.lsm.html '+mylsm[0])

			
# Mask the model images around the troublemaker
if maskmodels:
	for chan_image in chan_images:
		mod_image = chan_image.replace('image','model')
		maskimage(mod_image,directions,extent)

# Run wsclean in PREDICT mode on new model images
if runwspredict:
	for myms in mslist:
		predict(myms,imgbase)
