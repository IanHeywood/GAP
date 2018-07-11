import numpy
import glob
import os
import sys
from astropy.io import fits as pyfits

def gi(message):
        print '\033[92m'+message+'\033[0m'

def ri(message):
        print '\033[91m'+message+'\033[0m'

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
        f = pyfits.open(fitsfile,mode='update')
        input_hdu = f[0]
        if len(input_hdu.data.shape) == 2:
                input_hdu.data[:,:] = newimage
        elif len(input_hdu.data.shape) == 3:
                input_hdu.data[0,:,:] = newimage
        else:
                input_hdu.data[0,0,:,:] = newimage
        f.flush()

# image prefixes
resid = sys.argv[2]
model = sys.argv[1]

modelresid = model+'-MFS-residual.fits'
modelfull = model+'-MFS-image.fits'
resid = resid+'-MFS-image.fits'
modelonly = modelfull.replace('.fits','_modconv.fits')
restored = resid.replace('.fits','_restored.fits')

gi('Subtracting:')
gi('      '+modelresid)
gi('from:')
gi('      '+modelfull)
gi('to give:')
gi('      '+modelonly)

mr = getImage(modelresid)
mf = getImage(modelfull)
mo = mf - mr
os.system('cp '+modelresid+' '+modelonly)
flushFits(mo,modelonly)

gi('Restoring:')
gi('      '+modelonly)
gi('into:')
gi('      '+resid)
gi('to give:')
gi('      '+restored)

res = getImage(resid)
os.system('cp '+resid+' '+restored)
im = mo+res
flushFits(im,restored)

