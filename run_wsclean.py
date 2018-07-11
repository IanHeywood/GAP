import subprocess
import glob
import sys

myms = sys.argv[1]
suffix = sys.argv[2]
column = sys.argv[3]
field = sys.argv[4]
automask = sys.argv[5]
updatemodel = sys.argv[6]
dryrun = sys.argv[7]

def str2bool(str):
	if str == 'True':
		ret = True
	elif str == 'False':
		ret = False
	return ret

dryrun = str2bool(dryrun)
automask = str2bool(automask)
updatemodel = str2bool(updatemodel)

cleanup = True

imgname = 'img_'+myms+'_'+suffix

syscall = 'wsclean '
syscall+= '-log-time '
syscall+= '-field '+str(field)+' '
syscall+= '-size 8192 8192 '
syscall+= '-scale 3.5asec '
syscall+= '-niter 30000 '
syscall+= '-gain 0.1 '
syscall+= '-mgain 0.85 '
syscall+= '-weight briggs 0.0 '
#syscall+= '-minuv-l 500 '
#syscall+= '-taper-inner-tukey 50 '
syscall+= '-datacolumn '+column+' '
if automask:
	syscall+= '-local-rms '
	syscall+= '-auto-threshold 0.3 '
	syscall+= '-auto-mask 5.0 '
syscall+= '-name '+imgname+' '
syscall+= '-channelsout 4 '
syscall+= '-fit-spectral-pol 4 '
syscall+= '-joinchannels '
if not updatemodel:
	syscall+= '-no-update-model-required '
syscall+= '-mem 50 '
syscall+= myms

print syscall

if not dryrun:
	subprocess.call(syscall,shell=True)

	rem = ['first-residual','psf','dirty']
	chan_rem = ['residual']

	if cleanup:
		for product in rem:
			xx = sorted(glob.glob('img_'+myms+'_'+suffix+'*'+product+'*.fits'))
			for item in xx:
				subprocess.call('rm '+item,shell=True)
		for product in chan_rem:
			xx = sorted(glob.glob('img_'+myms+'_'+suffix+'*00*'+product+'*.fits'))
			for item in xx:
				subprocess.call('rm '+item,shell=True)
