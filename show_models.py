import numpy
import Tigger
import os
import sys

prefix = sys.argv[1]

fitsfile = prefix+'.fits'
lsmfile = prefix+'.pybdsm.lsm.html'
#pngname = '/nfs/wwwpeople/hey036/ASKAP/pcal_models/'+prefix+'.png'
pngname = prefix+'.png'

srcs = Tigger.load(lsmfile).sources

syscall = 'mViewer -ct 0 -color yellow -grid equatorial -color cyan '

def r2d(x):
	return 180.0*x/numpy.pi

for src in srcs:
	ra_d = str(r2d(src.pos.ra))
	dec_d = str(r2d(src.pos.dec))
	dE = src.getTag('dE')
	if dE:
		syscall += '-color red '
	else:
		syscall += '-color yellow '
	syscall += '-symbol 0.50 circle -mark '+ra_d+' '+dec_d+' '

syscall += '-gray '+fitsfile+' -2s max gaussian-log '
syscall += '-png '+pngname

os.system(syscall)
