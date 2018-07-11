import numpy
import glob
import subprocess
import datetime
import time
import sys
from pyrap.tables import table

myms = sys.argv[1]
parset = sys.argv[2]
suffix = sys.argv[3]
field = sys.argv[4]
dryrun = sys.argv[5]

if dryrun == 'True':
	dryrun = True
elif dryrun == 'False':
	dryrun = False

now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
syscall = 'gocubical '+parset+' '
syscall += '--data-ms='+myms+' '
syscall += '--sel-field='+str(field)+' '
syscall += '--out-name cube_'+suffix+'_'+myms+'_'+now
print syscall
if not dryrun:
	subprocess.call(syscall,shell=True)
