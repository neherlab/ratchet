'''
Created on Oct 2, 2011

@author: richard
'''
from pylab import *
from math import factorial
from scipy import interpolate as ip
from scipy import optimize as op
import argparse

def dn(n, lam):
    '''
    Time derivative of the bin sizes
    '''
    ndot = zeros_like(n)
    meanfit = sum(arange(len(n))*n)/sum(n)
    ndot+=(meanfit-arange(len(n))-lam)*n
    ndot[1:]+=lam*n[:-1]
    return ndot

def solve_n(n0,dt, lam):
    '''solve Haigh's deterministic model using a simple forward Euler scheme
    while fixing the size of the top bin
    Arguments are the trajectory of the top bin n0, the time step and lambda = u/s
    '''
    kmax=int(ceil(10+4*lam))
    n = array([exp(-lam)*lam**k/factorial(k) for k in range(kmax)])
    n[0]=n0[0]
    mk = []
    
    for n0p in n0:
        n[0]=n0p
        n/=sum(n)
        mk.append(sum(arange(kmax)*n))
        n+=dt*dn(n, lam)
    
    #return distribution and the mean fitness trajectory    
    return n, array(mk)


def Action(n0, dt, lam):
    '''given a n_0 trajectory, calculate the associated action
    parameters: the trajecory n0, the time steps of in between values of the trajectory, and lambda'''
    
    #calculate the trajectory of the mean fitness
    n,mk=solve_n(n0, dt, lam)
    #integrate the action
    return 0.5*dt*sum(((n0[1:]-n0[:-1])/dt-(mk[:-1]-lam)*n0[:-1])**2/n0[:-1]/(1-n0[:-1]))


def make_n0(n0_pivots, t_pivots, dt):
    '''from the set of n_0, time pairs, produce a dense trajectory of n_0 by linear interpolation'''
    n0func=ip.interp1d(t_pivots, n0_pivots, kind='linear')
    tarray=linspace(min(t_pivots), max(t_pivots), (max(t_pivots)-min(t_pivots))/dt+1)
    return n0func(tarray)


def function_to_minimize(path_pivots, ta,na,tb,nb, lam,dt):
    '''
    wrapper that interpolates the path given the pivots and return the action
    '''
    if min(path_pivots)<0:
        return 1e10
    if max(path_pivots)>1:
        return 1e10
    n0_pivots=[na]
    n0_pivots.extend(path_pivots)
    n0_pivots.append(nb)
    t_pivots=[ta]
    t_pivots.extend(linspace(ta,tb,len(n0_pivots))[1:-1])
    t_pivots.append(tb)
    
    n0=make_n0(n0_pivots, t_pivots, dt)
    return Action(n0,dt,lam)

def return_path(path_pivots, ta,na,tb,nb, lam,dt):
    '''
    As above, this calculates the path given the pivots and returns the path along with the action
    '''
    if min(path_pivots)<0:
        return 1e10
    if max(path_pivots)>1:
        return 1e10
    n0_pivots=[na]
    n0_pivots.extend(path_pivots)
    n0_pivots.append(nb)
    t_pivots=[ta]
    t_pivots.extend(linspace(ta,tb,len(n0_pivots))[1:-1])
    t_pivots.append(tb)
    
    n0=make_n0(n0_pivots, t_pivots, dt)
    n,mk=solve_n(n0,dt, lam)
    S=Action(n0,dt,lam)
    #print S, path_pivots
    return S,n0,mk


#if __main__:
if True:
    parser=argparse.ArgumentParser(description="Find the optimal path to extinction that minimizes the path integral action")
    parser.add_argument('--lam',default=5, type=float, help="parameter lambda, ratio of mutation rate to mutation effect")
    parser.add_argument('--outdir',default='./',  help="directory into which output files are saved")
    params=parser.parse_args()
    if params.outdir[-1] is not '/': params.outdir=params.outdir+'/'
    lam=params.lam
    
    
    '''
    All calculations here are done in natural units with s=1 and N=1, 
    the result has to be multiplied by Ns to convert to actual units
    '''
    
    #set the time resolution for the numerical integration
    dt=0.005
    #set the endpoints and the time interval to be covered
    na=exp(-lam)
    ta=0
    nb=0 
    tb=20   #needs to large compared to log lamda, hence 20 should be plenty
    
    #determine the path, initially with few intermediate points
    npivots=4
    n0_init_pivots = (na*linspace(1,0, npivots))[1:-1]
    for npivots in range(4,25,2):    
        n0_opt, A= op.fmin(function_to_minimize, n0_init_pivots, args=(ta,na,tb,nb,lam,dt), full_output=1, ftol=1e-4*exp(-lam))[:2]
        S, n0, mk = return_path(n0_opt, ta,na,tb,nb, lam,dt)
        print "action:",S
        #generate the starting pivots for the next iteration
        n0_init_pivots = n0[::int(tb/dt/(npivots+1))][1:-1]
    
    file = open(params.outdir+'escape_path_lambda_'+str(lam)+'.dat', 'w')
    file.write('# lambda '+str(lam)+'\tS='+str(S)+'\t')
    for n,f in zip(n0,mk): file.write(str(n/na)+'\t'+str(f)+'\n')
    file.close()
    
    file = open(params.outdir+'Smin_lambda_'+str(lam)+'.dat', 'w')
    file.write(str(lam)+'\t'+str(S)+'\t')
    file.close()
    
