from pylab import *
from scipy import stats, special
import glob


params = {'backend': 'ps',  
          'axes.labelsize': 18, 
          'text.fontsize': 18,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 14,
'xtick.labelsize': 14,
'ytick.labelsize': 14,
'text.usetex': True}
rcParams.update(params)
execfile('correlation_functions.py')

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


kbarset=set()
sset=set()
n0sset=set()
correlation_functions={}

file_list = glob.glob('DataTrajectories/trajectory_dt*')
for filename in file_list:   
    params = parse_file_name(filename[:-4])
    N=float(params['N'])
    s = float(params['sel'])
    kbar = float(params['lam'])
    n0s = round(N*s*exp(-kbar))    
    label = (n0s,kbar, s)

    kbarset.add(kbar)
    n0sset.add(kbar)
    sset.add(s)
    
    temp=loadtxt(filename)
    mean_k = mean(temp[:,0])
    mean_n0 = mean(temp[:,1])

    dt=int(params['dt'])
    #calculate autp correlations    
    data_n0_acorr = [mean((temp[:,1]-mean_n0)*temp[:,1])]
    data_k_acorr = [mean((temp[:,0]-mean_k)*temp[:,0])]
    data_n0_acorr.extend([mean((temp[i:,1]-mean_n0)*(temp[:-i,1]-mean_n0)) for i in range(1,60)])
    data_k_acorr.extend([mean((temp[i:,0]-mean_k)*(temp[:-i,0]-mean_k)) for i in range(1,60)])
    
    #calculate cross correlation
    data_cross_corr = [mean((temp[i:,1]-mean_n0)*(temp[:-i,0]-mean_k)) for i in range(100,0,-1)]
    data_cross_corr.append(mean((temp[:,1]-mean_n0)*(temp[:,0]-mean_k)))
    data_cross_corr.extend([mean((temp[:-i,1]-mean_n0)*(temp[i:,0]-mean_k)) for i in range(1,100)])
    
    #save suitably normalized correlation functions.
    #note that the mean fitness is in units of s, hence we need to divide data_k_acorr by s^2 
    correlation_functions[label]=[array(data_n0_acorr)/n0s*s**2, array(data_k_acorr)*n0s/s**2, array(data_cross_corr), dt]



kbarlist = sorted(list(kbarset))
slist = sorted(list(sset))


#plot top bin correlation function
figure()
n0s=100
linestyles = ['--', ':']
for si,s in enumerate(slist[1:]):
    for kbar in kbarlist:
        label=(n0s, kbar, s)
        tarray=linspace(0,6,61)
        dt = correlation_functions[label][-1]
        try:
            c=cm.jet(int(255.0/log(max(kbarlist))*log(kbar)))
            plot(tarray[:-1], correlation_functions[label][0], color=c, label=r'$\lambda='+str(kbar)+'$', linewidth=2, linestyle=linestyles[si])
            theory=[]
            for t in tarray: theory.append(n0_acorr(kbar,t))
            plot(tarray, theory, linestyle='-', color=c, linewidth=2)
        except:
            print label

xlabel(r'$\tau\; [s^{-1}]$')
ylabel(r'$Nse^{\lambda}\langle \delta x_0(0)\delta x_0(\tau)\rangle$')
text(-0.7, 3, 'A', fontsize=22)
ylim([0,3])
legend()
savefig('Figures/n0_correlation.pdf')


#plot top bin meanfitness function
figure()
n0s=100
s=0.01
for kbar in kbarlist:
    label=(n0s, kbar, s)
    tarray=linspace(0,6,61)
    dt = correlation_functions[label][-1]
    try:
        c=cm.jet(int(255.0/log(max(kbarlist))*log(kbar)))
        plot(tarray[:-1], correlation_functions[label][1], color=c, label=r'$\lambda='+str(kbar)+'$', linewidth=2, linestyle=linestyles[si])
        theory=[]
        for t in tarray: theory.append(kbar_acorr(kbar,t))
        plot(tarray, theory, linestyle='-', color=c, linewidth=2)
    except:
        print label
xlabel(r'$\tau \; [s^{-1}]$')
ylabel(r'$Nse^{-\lambda}\langle \delta \bar{k}(0)\delta \bar{k}(\tau)\rangle$')
legend()
savefig('Figures/mean_fitness_correlation.pdf')


#plot cross correlation function
figure()
n0s=100
s=0.01
for kbar in kbarlist:
    label=(n0s, kbar, s)
    tarray=linspace(-10,10,201)
    dt = correlation_functions[label][-1]
    try:
        c=cm.jet(int(255.0/log(max(kbarlist))*log(kbar)))
        plot(tarray[:-1], correlation_functions[label][2], color=c, label=r'$\lambda='+str(kbar)+'$', linewidth=2, linestyle=linestyles[si])
        theory=[]
        for t in tarray: theory.append(n0_kbar_cov(kbar,t))
        plot(tarray, theory, linestyle='-', color=c, linewidth=2)
    except:
        print label
xlabel(r'$\tau \; [s^{-1}]$')
ylim([0,1])
ylabel(r'$Ns\langle  \delta x_0(0)\delta \bar{k}(\tau)\rangle$')
text(-12, 1, 'B', fontsize=22)
legend()
savefig('Figures/cross_correlation.pdf')



