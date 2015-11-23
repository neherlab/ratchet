'''
Created on May 18, 2011

@author: richard
'''
from __future__ import division
from pylab import *
from scipy import stats
from math import factorial
import sys


def speed(popsize, kbar, dfit, mu, trun=1e6, teq=10000, poisson=True):
    '''
    Calculate the speed of Muller's ratchet for a given set of population parameters
    allows to choose between Poisson and multinomial resampling
    '''
    #initialize and equilibrate
    pop, fitness, meanfitness = equilibrate(popsize, kbar, dfit, mu, teq, poisson)
    #memorize the current fitness and evolve for trun generations
    initial_meanfitness=meanfitness
    pop, fitness, meanfitness = evolve_t_generations(trun, pop, fitness, popsize, mu, dfit, poisson)
    
    #return the average speed
    return (meanfitness-initial_meanfitness)/trun


def trajectories(popsize, kbar, dfit, mu, trun=1e6, teq=10000,dt=1, poisson=True):
    '''
    Provide trajectories of the mean fitness and the size of the top bin
    '''
    #initialize and equilibrate
    pop, fitness, meanfitness = equilibrate(popsize, kbar, dfit, mu, teq, poisson)

    #get trajectories
    pop, fitness, meanfitness_traj, nose_traj = get_trajectories(trun, pop, fitness, popsize, mu, dfit,dt, poisson)

    #return trajectories
    return meanfitness_traj, nose_traj


def equilibrate(popsize, kbar, dfit, mu, teq=10000, poisson=True):
    '''
    set up and equilibrate the population
    '''
    #set grid_length
    sum_pmf = 1
    grid=0
    while (sum_pmf> 0.0001):
        sum_pmf-=stats.poisson.pmf(grid,kbar)
        #print grid, sum_pmf
        grid+=1
    
    grid*=2
    
    #initialize the fitness distribution and the population
    fitgrid = dfit*arange(-grid,1)
    pop = popsize*stats.poisson.pmf(arange(grid+1),kbar)[::-1]
    
    #equilibrate (we start in the steady state, so this shouldn;t matter much)
    pop, fitness, meanfitness = evolve_t_generations(teq, pop, fitgrid, popsize, mu, dfit, poisson)   
    return pop, fitness, meanfitness

def evolve_t_generations(iterations, initialpop, initialfitness, targetN, mu, dfit, poisson=True):
    '''
    move the population forward for a specified number of generations
    '''
    fitness=array(initialfitness)
    pop=array(initialpop)
    meanfitness = sum(fitness*pop)/sum(pop)
    
    #decide on evolution function to use
    if poisson:
        evofunc=evolve_ratchet
    else:
        evofunc= evolve_ratchet_multinomial
    
    for t in xrange(int(iterations)):
        pop = evofunc(pop, fitness, targetN, mu)
        meanfitness = sum(fitness*pop)/sum(pop)

        #shift grid to center on average
        k=1
        while (pop[-k]==0):
            k+=1
        if(k>1):
            pop[(k-1):]=pop[:-(k-1)]
            fitness-=dfit*(k-1)
        
    return pop, fitness, meanfitness

def get_trajectories(iterations, initialpop, initialfitness, targetN, mu, dfit,dt=1, poisson=True):
    '''
    function that returns trajectories for future analysis. It need to be provided with an initial distribution that
    is evolved over iterations generations. Note that this might crash due to memory constraints in iterations is large than 
    1e7 or 1e8
    '''
    #allocate arrays to store mean fitness and nose size in 
    meanfitness_traj = zeros(iterations/dt)
    nose_traj = zeros(iterations/dt)
    fitness=array(initialfitness)
    pop=array(initialpop)
    
    #decide on evolution function to use
    if poisson:
        evofunc=evolve_ratchet
    else:
        evofunc= evolve_ratchet_multinomial
    
    
    #iterate and save trajectories
    count=0
    for t in xrange(int(iterations)):
        pop = evofunc(pop, fitness, targetN, mu)
        meanfitness = sum(fitness*pop)/sum(pop)
        #shift grid to center on average
        k=1
        while (pop[-k]==0):
            k+=1
        if(k>1):
            pop[(k-1):]=pop[:-(k-1)]
            fitness-=dfit*(k-1)

        if (t%dt)==0:
            meanfitness_traj[count]=meanfitness
            nose_traj[count]=pop[-1]
            count+=1
        
    return pop, fitness, meanfitness_traj,nose_traj


def evolve_ratchet(population, fit, targetN,mu):
    '''
    Function that takes a population stratified according to the fitness grid fit, selects on it, mutates, and
    resamples it according to a Poisson distribution. The population size is kept constant on average by a carrying capacity targetN
    '''
    popsize=sum(population)
    fit_average = sum(population*fit)/popsize
    exp_fit_average = sum(population*exp(fit-fit_average))/popsize
    population_regularizer = 0.5*(1.0-1.0*popsize/targetN)
    #Selection step
    newpopulation = population*exp(fit-fit_average+population_regularizer)/exp_fit_average

    #mutation step (determine the maximally plausible number of mutations and the associated Poisson weights)
    max_mut_number = min((6+int(mu+5*sqrt(mu))), len(population))
    mut_dis=stats.poisson.pmf(arange(max_mut_number),mu)

   #redistribute the population according to these Poisson weights
    mutants = zeros(newpopulation.shape)
    mutants += mut_dis[0]*newpopulation
    for k in range(1,max_mut_number):
        mutants[:-k] += mut_dis[k]*newpopulation[k:]
        
    #poisson resample the population (approximate the Poisson by a Gaussian whenever the number is greater than 10000)
    gauss_indicies = where(mutants>=10000)
    newpopulation[where(mutants<10000)] =poisson(mutants[where(mutants<10000)])
    newpopulation[gauss_indicies] =mutants[gauss_indicies] + sqrt(mutants[gauss_indicies]) * randn(len(gauss_indicies[0]))
    return newpopulation


def evolve_ratchet_multinomial(population, fit, targetN,mu):
    '''
    Function that takes a population stratified according to the fitness grid fit, selects on it, mutates, and
    produces a multinomially distributed sample
    '''
    popsize=sum(population)
    fit_average = sum(population*fit)/popsize
    exp_fit_average = sum(population*exp(fit-fit_average))/popsize
    #Selection step
    newpopulation = population*exp(fit-fit_average)/exp_fit_average

    #mutation step (determine the maximally plausible number of mutations and the associated Poisson weights)
    max_mut_number = min((6+int(mu+5*sqrt(mu))), len(population))
    mut_dis=stats.poisson.pmf(arange(max_mut_number),mu)

    #redistribute the population according to these Poisson weights
    mutants = zeros(newpopulation.shape)
    mutants += mut_dis[0]*newpopulation
    for k in range(1,max_mut_number):
        mutants[:-k] += mut_dis[k]*newpopulation[k:]
    
    #renormalize
    mutants/=sum(mutants)
    #generate a sample multinomial distribution by successively sampling binomially "bin k or the rest"
    #readjust the remaining sample size and the success probability in each round
    newpopulation=zeros_like(population)
    tempN=targetN;
    tempNorm=1.0
    for k in range(1,newpopulation.shape[0]):
        #calculate the renormalized success probability and check for rounding errors
        p = min(1,mutants[-k]/tempNorm)
        p = max(0,p)
        #sample
        newpopulation[-k]=stats.binom.rvs(tempN,p)
        #calculate the remaining individuals to sample and the remaining probability mass
        tempN-=newpopulation[-k]
        tempNorm-=mutants[-k]
        if (tempN==0): break
    
    return newpopulation

