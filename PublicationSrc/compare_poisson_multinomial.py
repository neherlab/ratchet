from pylab import *
import glob

params = {'backend': 'ps',  
          'axes.labelsize': 16, 
          'text.fontsize': 16,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 16,
'xtick.labelsize': 14,
'ytick.labelsize': 14,
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




#read data
rates=[{},{}]
clicks=[{},{}]
file_list = glob.glob('DataPoissonMultinomial/ratchet_rate_N*')
for filename in file_list:
    params = parse_file_name(filename[:-4])
    s = params['sel']
    lam = params['lam']
    N=int(params['N'])
    poi = int(params['Poisson'])
    data = loadtxt(filename)
    
    label = (N,lam, s)
    rates[poi][label] = data[-1]
    clicks[poi][label] = data[-2]*data[-1]

poisson_over_multinomial=[]
for label in rates[1]:
    try:
        poisson_over_multinomial.append([ label[0], label[1], label[2], rates[1][label], rates[0][label],rates[1][label]/rates[0][label], clicks[1][label]])
    except KeyError:
        print "Multinomial: Data point",label,"is missing"
poisson_over_multinomial=array(poisson_over_multinomial)


well_sampled_rates = where(poisson_over_multinomial[:,-1]>10)
figure()
plot([1e-7, 1e-2],[1e-7, 1e-2], color='k', linewidth=2)
scatter(poisson_over_multinomial[well_sampled_rates,3], poisson_over_multinomial[well_sampled_rates,4],c=log(poisson_over_multinomial[well_sampled_rates,0]), s=55)
colorbar()
xlabel("Muller's ratchet rate, Poisson Sampling")
ylabel("Muller's ratchet rate, Multinomial Sampling")
ax=gca()
ax.set_yscale('log')
ax.set_xscale('log')
xlim([1e-7,1e-2])
ylim([1e-7,1e-2])
savefig('Figures/comparePoissonMultinomial.pdf')

