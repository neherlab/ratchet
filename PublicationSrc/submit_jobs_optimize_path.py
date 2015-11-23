#!/usr/bin/python                                                                                                           
import os
import glob
from pylab import *

p=50

outdir = 'DataOptimalPath'
try:
    os.mkdir(outdir)
except:
    print "Cannot make directory"

for kbar in arange(1,31):
    qsub_command = 'qsub -cwd -p '+str(p)+' submit_script_optimize.py --lam '+str(kbar)+' --outdir '+outdir 
    print qsub_command
    os.system(qsub_command)

