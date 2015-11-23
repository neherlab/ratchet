#!/usr/bin/python                                                                                                           
import os
import glob
from pylab import *

slist = [0.001,  0.01]
p=50

outdir = 'DataPoissonLong'
try:
    os.mkdir(outdir)
except:
    print "Cannot make directory"

for kbar in [2,3,4,5,7,10,15,20,25,30]:
    for s in slist:
        for N in linspace(exp(kbar)/s, 29*exp(kbar)/s, 15):
            if (N*exp(-kbar)*s>0.5) and (N*exp(-kbar)*s<30):
                qsub_command = 'qsub -cwd -p '+str(p)+' submit_script.py --pop '+str(int(N))+' --sel '+str(s)+' --lam '+str(kbar)+' --Trun 1e8 --Teq 1e4 --poisson 1 --speed --outdir '+outdir 
                print qsub_command
                os.system(qsub_command)

