from pylab import *
from scipy import stats, special, integrate, optimize

import glob
execfile('correlation_functions.py')

params = {'backend': 'ps',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 20,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
rcParams.update(params)

def parse_file_name(fname):
    params = {}
    entries=fname.split('_')
    print entries
    for i in range(len(entries)):
        try:
            params[entries[-i-2]]=float(entries[-i-1])
        except:
            pass
    return params


#load the ratchet rates determined by simulations
rates={}
clicks={}
file_list = glob.glob('DataPoissonLong/ratchet_rate_N*')
for filename in file_list:   
    params = parse_file_name(filename[:-4])
    s = float(params['sel'])
    lam = float(params['lam']) 
    N=int(params['N'])
    sigma=s*sqrt(lam)
    data=loadtxt(filename)
    label=(N,lam,s)
    rates[label] = data[-1]
    clicks[label] = data[-1]*data[-2]

 
#load the numerically determined minimal action and the associated optimal escape path
action={}
path={}
file_list = glob.glob('DataOptimalPath/escape_path_lambda*')
for filename in file_list:
    params = parse_file_name(filename[:-4])
    lam = float(params['lambda']) 
    path[lam]=loadtxt(filename)

file_list = glob.glob('DataOptimalPath/Smin_lambda*')
for filename in file_list:
    params = parse_file_name(filename[:-4])
    lam = float(params['lambda']) 
    data=loadtxt(filename)
    action[lam] = data[-1]

#summarize all data in a big list
rates_array=[]
for label in rates:
    N,lam,s=label
    n0s = N*exp(-lam)*s
    n0var=n0_acorr(lam, 0)

    if lam>0.5 and clicks[label]>50:
        rates_array.append([N, lam, s, n0s, rates[label],clicks[label], n0var, N*s*action[lam]])

Ni=0
lami=1
si=2
n0si=3
ratei=4
clicki=5
n0vari=6
Smini=7
ra=array(rates_array)
ms=45   #markersize

# produce figure comparing the simulated rate with the traditional Stephan/Gordo-Charlesworth/Jain formula
figure()
haigh = 0.5
prefactor = haigh*ra[:,si]*sqrt(ra[:,n0si]*haigh)/sqrt(pi)
scatter(ra[:,n0si]*haigh, (ra[:,ratei]/prefactor), c=(ra[:,lami]), s=ms)
title(r"Haigh's factor $\alpha="+str(haigh)+"$")
plot([0,10], exp(array([0,-10])), 'k', linewidth=2)
ax=gca()
ax.set_yscale('log')
xlim([0,32*haigh])
ylim(exp(array([-10,0])))
xlabel(r'$ \alpha Ns\bar{x}_0 $')
ylabel(r'$\gamma\times \sqrt{\pi/(\bar{x}_0 Ns \alpha^3)}$')
colorbar()
savefig('Figures/rates_GC.pdf')


# produce figure comparing the simulated rate to the approximation put forward in this manuscript
figure()
prefactor = ra[:,Smini]*ra[:,si]/sqrt(8*pi*ra[:,n0si]*ra[:,n0vari])
scatter(ra[:,Smini], (ra[:,ratei]/prefactor), c=(ra[:,lami]), s=ms)
plot([0,10], exp(array([0,-10])), 'k', linewidth=2)
ax=gca()
ax.set_yscale('log')
ylim(exp(array([-11,2])))
xlim([0,10])
xlabel(r'$Ns\; S^*_{\lambda}(0)$')
ylabel(r'$\gamma\times \sqrt{8\pi \sigma^2}/S^*_{\lambda}$')
colorbar()
savefig('Figures/rates_Smin.pdf')

# produce figure comparing the simulated rate to the approximation using the inverse lambda dependent variance as Haigh's factor
figure()
prefactor = ra[:,si]/ra[:,n0vari]*sqrt(ra[:,n0si]/ra[:,n0vari])/sqrt(pi)
scatter(ra[:,n0si]/ra[:,n0vari], (ra[:,ratei]/prefactor), c=(ra[:,lami]), s=ms)
plot([0,10], exp(array([0,-10])), 'k', linewidth=2)
ax=gca()
ax.set_yscale('log')
xlim([0,10])
ylim(exp(array([-10,1])))
xlabel(r'$Ns\bar{x}_0/\sigma^2$')
ylabel(r'$\frac{\gamma \times \sqrt{2\pi}}{s\sigma^{-2}\sqrt{\bar{x}_0 Ns\sigma^{-2}}}$')
colorbar()
savefig('Figures/rates_n0var.pdf')



# plot optimal escape path
figure()
dt=0.005
for lam in [2,10,25]:
    t50_n0 = (argmax(path[lam][:,0]<0.5)*dt)
    plot((linspace(0,20,path[lam].shape[0])-t50_n0), path[lam][:,0],color = cm.jet(int(255.0/13*lam)), linestyle='-', linewidth=2, label=r'$\lambda='+str(lam)+'$')
    plot((linspace(0,20,path[lam].shape[0])-t50_n0), 1-(path[lam][:,1]-lam),color = cm.jet(int(255.0/13*lam)) , linestyle='--', linewidth=2)

xlabel(r'$\tau\; [s^{-1}]$')
ylabel(r'$x_0(\tau)/\bar{x}_0$ (solid), $1-\delta \bar{k}(\tau)$ (dashed)')
legend(loc=1)
xlim([-8,8])
ylim([0, 1.05])
ax=axes([0.19,0.19, 0.3, 0.35])
for lam in range(2,30):
    plot(path[lam][:,0], (path[lam][:,1]-lam), color = cm.jet(int(255.0/30*lam)))

ylim([0,1])
xlim([0,1])
xticks([0,1])
yticks([0,1])
#colorbar()
xlabel(r'$x_0(\tau)/\bar{x}_0$')
ylabel(r'mean fitness $\delta k(\tau)$')
savefig('Figures/n0_and_dk_trajectories.pdf')

# compare Smin to the effective Haigh parameter
figure()
haigh = []
for lam in range(2,30):
    haigh.append([lam, action[lam]*exp(lam)])

haigh=array(haigh)
plot(haigh[:,0], haigh[:,1], label=r'$S^*_{\lambda}e^{\lambda}$', linewidth=2, color='k')
#plot([2,10], [0.6, 0.6], label='Stephan', linestyle='--', linewidth=2)
#plot([2,10], [0.5, 0.5], label='Gordo-Charlesworth', linestyle='-.', linewidth=2)
xlim([0,30])
ylim([0,0.8])
xlabel(r'$\lambda$')
ylabel(r"Haigh's factor $\alpha(\lambda)$")
legend()
savefig('Figures/approx_Smin.pdf')

#compare explicit formula for rate
# compare Smin to the effective Haigh parameter

figure()
scatter(1/ra[:,ratei], 2.5*sqrt(ra[:,n0vari]*ra[:,n0si])/ra[:,si]/ra[:,Smini]*exp(ra[:,Smini]) ,c=(ra[:,Smini]), s=ms)
xlabel(r'$T_{click}$')
ylabel(r"Theory")
text(2000,2e6,r'$T_{click}\approx \frac{2.5\zeta(\lambda)}{\alpha(\lambda)s \sqrt{Nse^{-\lambda}}}e^{Ns \alpha(\lambda)e^{-\lambda}}$', fontsize=22)
ax=gca()
ax.set_yscale('log')
ax.set_xscale('log')
xlim([1e3,1e7])
ylim([1e3,1e7])
plot([1e3,1e7], [1e3,1e7])
colorbar()
legend()
savefig('Figures/explicit_theory_simulation_comparison.pdf')





