#!/usr/bin/python                                                                                                           
import os
import glob
from pylab import *

slist = [0.001,  0.01]
n0s_list = [100,1000]
p=50

outdir = 'DataTrajectories'
try:
    os.mkdir(outdir)
except:
    print "Cannot make directory"

for kbar in [1,2,3,4,5,7,10,15]:
    for s in slist:
        for n0s in n0s_list:
            N=n0s*exp(kbar)/s

            qsub_command = 'qsub -cwd -p '+str(p)+' submit_script.py --pop '+str(int(N))+' --sel '+str(s)+' --lam '+str(kbar)+' --Ttraj '+str(1e4/s)+' --Teq 1e4 --poisson 1 --traj --dt '+str(int(1.0/(10*s)))+' --outdir '+outdir 
            print qsub_command
            os.system(qsub_command)

