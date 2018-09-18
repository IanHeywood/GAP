import glob
import subprocess
import numpy
from pyrap.tables import table
import datetime
import time
import sys


def getchans(msname,nfreqs):
	tt = table(msname.rstrip('/')+'/SPECTRAL_WINDOW')
	freqs = tt.getcol('CHAN_FREQ')
	nchans = len(freqs[0])
	dchan = nchans/nfreqs
	chansels = []
	for i in range(0,nfreqs):
		ch0 = i*dchan
		if i == nfreqs-1:
			ch1 = nchans-1
		else:
			ch1 = ((i+1)*dchan)-1
		chansels.append((ch0,ch1))
	tt.done()
	return chansels


myms = sys.argv[1]
parset = sys.argv[2]
suffix = sys.argv[3]
field = sys.argv[4]
dryrun = sys.argv[5]


nfreqs = 4


if dryrun == 'True':
	dryrun = True
elif dryrun == 'False':
	dryrun = False


lsmlist = sorted(glob.glob('*000*.html'))
chansels = getchans(myms,nfreqs)


for i in range(0,nfreqs):
	now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
	ch0 = chansels[i][0]
	ch1 = chansels[i][1]
	mylsm = lsmlist[i]
	syscall = 'gocubical '+parset+' '
	syscall += '--data-ms='+myms+' '
	syscall += '--sel-field='+str(field)+' '
	syscall += '--sel-chan='+str(ch0)+'~'+str(ch1)+' '
	syscall += '--model-list=MODEL_DATA:'+lsmlist[i]+'@dE '
	syscall += '--g-save-to={data[ms]}/G-field:{sel[field]}:fchunk'+str(i)+' '
	syscall += '--de-save-to={data[ms]}/dE-field:{sel[field]}:fchunk'+str(i)+' '
	syscall += '--de-fix-dir 0 '
	syscall += '--out-name cube_'+suffux+'_'+myms+'_'+now+'_fchunk'+str(i)
	if not dryrun:
		subprocess.call(syscall,shell=True)
