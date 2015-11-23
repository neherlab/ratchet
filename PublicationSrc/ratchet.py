'''
Created on May 18, 2011

@author: richard
'''
from __future__ import division
from pylab import *
from evolve_methods import *
from scipy import optimize
import sys
import argparse

#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a deleterious mutation selection balance and determine the rate of Muller's ratchet")
parser.add_argument('--pop', default=10000, type=float, help='Population size (N)')
parser.add_argument('--lam', default=None,type=float, help='Lambda (u/s)')
parser.add_argument('--sel',  default=None,type=float, help='Selection strength (s)')
parser.add_argument('--sigma', default=None,type=float, help='Sigma (sqrt(u*s))')
parser.add_argument('--mu', default=None,type=float, help='Mutation rate u')
parser.add_argument('--Trun', default=1000000,type=float, help='Time over which the ratchet rate is measured in generations')
parser.add_argument('--Ttraj', default=None,type=float, help='Length of trajectory in generations')
parser.add_argument('--Teq', default=1000,type=float, help='Equilibration time in generations')
parser.add_argument('--poisson',  default=1,type=int, help='Poisson [1] sampling vs multinomial [0]')
parser.add_argument('--speed',  default=0,action='count', help='Measure speed')
parser.add_argument('--traj',  default=0,action='count', help='Store trajectory')
parser.add_argument('--dt',  default=1,type =int, help='time increments of trajectory')
parser.add_argument('--outdir',  default='./', help='Directory into which output is directed')
params=parser.parse_args()

#there are multiple redundant was in which the parameters can be specified
#the following three are supported
if params.lam is not None and params.sigma is not None:
    params.sel = params.sigma/sqrt(params.lam)
    params.mu = params.sigma*sqrt(params.lam)
elif params.lam is not None and params.sel is not None:
    params.sigma = params.sel*sqrt(params.lam)
    params.mu = params.sigma*sqrt(params.lam)
elif params.sel is not None and params.mu is not None:
    params.sigma = sqrt(params.sel*params.mu)
    params.lam = params.mu/params.sel
else:
    print "Specify either lambda and sigma, lambda and s, or s and u"
    exit()
if params.Ttraj is None: params.Ttraj=params.Trun
if params.outdir[-1] is not '/': params.outdir=params.outdir+'/'

params.pop=int(params.pop)
params.Teq=int(params.Teq)
params.Trun=int(params.Trun)
params.Ttraj=int(params.Ttraj)

print 'Parameters:',params

filestr='N_'+str(params.pop)+'_sel_'+str(params.sel)+'_lam_'+str(params.lam)+'_Poisson_'+str(int(params.poisson))+'.dat'
if params.speed:
    print 'Calculate ratchet rate'
    v =speed(params.pop,params.lam, -log(1-params.sel), params.mu  ,params.Trun, params.Teq, params.poisson)
    
    #output the result to file, the file is labelled by s and lambda
    output_file = open(params.outdir+'ratchet_rate_'+filestr,'w')
    output_file.write(str(params.pop)+'\t'+str(params.sel)+'\t'+str(params.mu)+'\t'+str(params.lam)+'\t'+str(params.Trun)+'\t'+str(-v/params.sel))
    output_file.write('\n')
    output_file.close()

if params.traj:
    print 'Save trajectories of nose and mean fitness'
    meanfitness_traj, nose_traj =trajectories(params.pop,params.lam, -log(1-params.sel), params.mu  ,params.Ttraj, params.Teq,params.dt, params.poisson)
    savetxt(params.outdir+'trajectory_dt_'+str(params.dt)+'_'+filestr, array([meanfitness_traj, nose_traj]).T)


