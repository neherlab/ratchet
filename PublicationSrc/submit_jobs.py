#!/usr/bin/python                                                                                                           
import os
import glob
from pylab import *

slist = [0.001,  0.01]
p=50

outdir = 'DataPoissonMultinomial'
try:
    os.mkdir(outdir)
except:
    print "Cannot make directory"

for kbar in [2,5,8]:
    for s in slist:
        for N in linspace(exp(kbar)/s, 30*exp(kbar)/s, 30):
            if (N*exp(-kbar)*s>1) and (N*exp(-kbar)*s<30):
                qsub_command = 'qsub -cwd -p '+str(p)+' submit_script.py --pop '+str(int(N))+' --sel '+str(s)+' --lam '+str(kbar)+' --Trun 1e7 --Teq 1e4 --poisson 0 --speed --outdir '+outdir 
                print qsub_command
                os.system(qsub_command)

