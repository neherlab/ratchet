from pylab import *
import pickle
params = {'backend': 'ps',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 16,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
rcParams.update(params)

file=open('escape.data', 'r')
nose, mean_fitness= pickle.load(file)

N=5e7
s=0.01
kbar=10





fig=figure(figsize=(12,5))
fig.subplots_adjust(bottom=0.15)

subplot(132)
plot(arange(3760)*10*s, nose[:3760]/N*exp(kbar), color='k')
plot([0,40000*s], [1,1], c='r', linewidth=2)
xticks([0,10000*s,20000*s, 30000*s])
ylim([0,2.2])
text(850*s, 1.97, 'B', fontsize=22)
text(22000*s, 1.97, r'$x_0(\tau)/\bar{x}_0$', fontsize=22)
xlabel(r'time $\tau\; [s^{-1}]$')
#ylabel(r'$n/\bar{n}$')


subplot(133)
plot(arange(3630, 3760)*10*s, nose[3630:3760]/N*exp(kbar), color='k')
plot([36100*s, 37600*s], [1,1], c='r', linewidth=2)
ylim([0,2.2])
xlim([36300*s,37600*s])
text(36350*s, 1.97, 'C', fontsize=22)
text(37050*s, 1.97, r'$x_0(\tau)/\bar{x}_0$', fontsize=22)
xticks([36500*s,37000*s, 37500*s])
xlabel(r'time $\tau\; [s^{-1}]$')
yticks([])

subplot(131)
hist(nose/N*exp(kbar), 100, normed=True, color='k', facecolor='k')
text(0.05, 1.25, 'A', fontsize=22)
text(1.8, 1.25, r'$p(x_0)$', fontsize=22)
xlabel(r'size of fittest class $x_0/\bar{x}_0$')

savefig('Figures/click_illustration.pdf')
