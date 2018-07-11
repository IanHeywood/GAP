import numpy
import glob
import subprocess
from pyrap.tables import table 

mslist = glob.glob('*wtspec.ms')

for myms in mslist:
	dryrun = False

	tt = table(myms,ack=False)
	fields = numpy.unique(tt.getcol('FIELD_ID')).tolist()
	field = fields[0]
	tt.done()

	print myms+' using field '+str(field)

	# wsclean myms suffix column field automask updatemodel dryrun
	cc = 'python run_wsclean.py '+myms+' data DATA '+str(field)+' True True '+str(dryrun)
	subprocess.call(cc,shell=True)

	# cubical myms parset prefix field dryrun
	cc = 'python run_cubical.py '+myms+' askap-phasecal.parset pcal '+str(field)+' '+str(dryrun)
	subprocess.call(cc,shell=True)

	# wsclean myms suffix column field automask updatemodel dryrun
	cc = 'python run_wsclean.py '+myms+' pcal CORRECTED_DATA '+str(field)+' True True '+str(dryrun)
	subprocess.call(cc,shell=True)

	# cubical myms parset prefix field dryrun
        cc = 'python run_cubical.py '+myms+' askap-apcal.parset apcal '+str(field)+' '+str(dryrun)
        subprocess.call(cc,shell=True)

        # wsclean myms suffix column field automask updatemodel dryrun
        cc = 'python run_wsclean.py '+myms+' apcal CORRECTED_DATA '+str(field)+' True True '+str(dryrun)
        subprocess.call(cc,shell=True)

	# run pybdsm
	cc = 'python find_sources.py img_'+myms+'_apcal'
	subprocess.call(cc,shell=True)

	# identify dE sources
	cc = 'python identify_dE_sources.py '+myms
	subprocess.call(cc,shell=True)

	# collapse models
	cc = 'python collapse_dE_IDs.py '+myms
	subprocess.call(cc,shell=True)

	collapse = glob.glob('*'+myms+'*.txt')
	if len(collapse) > 0:

		# setup dE models
		cc = 'python setup_dE_models.py '+myms
		subprocess.call(cc,shell=True)

		# dE calibration
	        cc = 'python run_cubical.py '+myms+' askap-dEcal.parset dEcal '+str(field)+' '+str(dryrun)
	        subprocess.call(cc,shell=True)	

		# image residuals
	        cc = 'python run_wsclean.py '+myms+' dEcal CORRECTED_DATA '+str(field)+' False False '+str(dryrun)
	        subprocess.call(cc,shell=True)

		# restore models
		model = glob.glob('img_'+myms+'*pcal*-MFS-image.fits')[0].split('-MFS')[0]
		resid = glob.glob('img_'+myms+'*dEcal*-MFS-image.fits')[0].split('-MFS')[0]
		cc = 'python restore_model.py '+model+' '+resid
		subprocess.call(cc,shell=True)
