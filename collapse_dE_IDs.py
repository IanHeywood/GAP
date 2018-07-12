import Tigger
import glob
import numpy
import sys

myms = sys.argv[1]

xx = glob.glob('*'+myms+'*html')


def rad2deg(xx):
	return 180.0*xx/numpy.pi

rf = open('setup_dE_runfile.sh','w')

count = 0
limit = 4

for lsm in xx:
	opdirs = lsm+'_collapsed.txt'
	print >>rf,'python setup_dE_model.py '+opdirs
	f = open(opdirs,'w')
	srcs = Tigger.load(lsm,verbose=False).sources
	for src in srcs:
		dE = src.getTag('dE')
		lead = src.getTag('cluster_lead')
		if dE and lead:
			name = src.name
			ra = rad2deg(src.pos.ra)
			dec = rad2deg(src.pos.dec)
			if ra < 0.0:
				ra += 360.0
			if count < limit:
				print >>f,name,ra,dec
				count += 1
	f.close()

rf.close()

