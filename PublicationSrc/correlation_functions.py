import numpy as np
from scipy import integrate,optimize


def n0_kbar_cov_integrand(x,lam, dt):
    if dt>0:
        return lam*np.exp(-dt)*((x-1)*np.exp(exp(-dt)*lam*x**2-lam*(1+np.exp(-dt))*x)+np.exp(-np.exp(-dt)*lam*x))
    else:
        return lam*((x*np.exp(dt)-1)*np.exp(np.exp(dt)*lam*x**2-lam*(1+np.exp(dt))*x)+np.exp(-lam*x))
        
def n0_acorr_integrand(x,lam, dt):
    if dt>0:
        return (np.exp(np.exp(-dt)*lam*x**2-lam*(1+np.exp(-dt))*x)-np.exp(-np.exp(-dt)*lam*x)-np.exp(-lam*x)+1)/x
    else:
        return (np.exp(np.exp(dt)*lam*x**2-lam*(1+np.exp(dt))*x)-np.exp(-np.exp(dt)*lam*x)-np.exp(-lam*x)+1)/x

def kbar_acorr_integrand(x,lam, dt):
    if dt>0:
        return np.exp(np.exp(-dt)*lam*x**2-lam*(1+np.exp(-dt))*x)*np.exp(-dt)*(lam*x+(np.exp(-dt)*lam*x**2-lam*x)*(lam*x-lam))
    else:
        return np.exp(np.exp(dt)*lam*x**2-lam*(1+np.exp(dt))*x)*np.exp(dt)*(lam*x+(np.exp(dt)*lam*x**2-lam*x)*(lam*x-lam))



def n0_kbar_cov(lam, dt):
    return integrate.quad(n0_kbar_cov_integrand, 0, 1, args=(lam, dt))[0]

def n0_acorr(lam, dt):
    return integrate.quad(n0_acorr_integrand, 0, 1, args=(lam, dt))[0]

def kbar_acorr(lam, dt):
    return integrate.quad(kbar_acorr_integrand, 0, 1, args=(lam, dt))[0]

def n0_kbar_cov_min(dt, lam):
    return -integrate.quad(n0_kbar_cov_integrand, 0, 1, args=(lam, dt))[0]

def max_n0_kbar_cov(lam):
    return optimize.fmin(n0_kbar_cov_min, 0, args=(lam,))[0]


