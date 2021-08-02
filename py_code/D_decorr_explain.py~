import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit


def func(x,m):
  return m*x
  

plt.figure(figsize=(2.5,2))
ax = plt.axes([0.2,0.2,0.7,0.7])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
  
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

K_array = [50,100,200]
nu_array = [0.2,0.5,1,2,5,10,20]
J_array = [-1,-2,-5]

steps = len(nu_array)

D_decorr = np.zeros(len(K_array)*len(nu_array)*len(J_array))
D_decorr_err = np.zeros(len(K_array)*len(nu_array)*len(J_array))
phase_var = np.zeros(len(K_array)*len(nu_array)*len(J_array))


marker_array = ['o','v','s']

for j in range(len(J_array)):
  
  for i in range(len(K_array)):
    col = i/float(len(K_array)-1)
    c = (col,col,col)
    dat_D_decorr = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[1000],'K':[K_array[i]],'rateWnt':nu_array,'J0':[J_array[j]],'NeuronType':[2],'special':[2]},mode_return=1)

    dat_stats = scriptFrame('r','statistics',1,alpha=[[0],[0],[0],[0]],netConst={'N':[1000],'K':[K_array[i]],'rateWnt':nu_array,'J0':[J_array[j]],'NeuronType':[2],'special':[0]},mode_return=1)
    
    D_decorr[steps*(i+j*len(K_array)):steps*(i+j*len(K_array)+1)] = dat_D_decorr['D_decorr'][:,0,0]
    D_decorr_err[steps*(i+j*len(K_array)):steps*(i+j*len(K_array)+1)] = np.sqrt(dat_D_decorr['D_decorr'][:,0,1])	
    phase_var[steps*(i+j*len(K_array)):steps*(i+j*len(K_array)+1)] = np.sqrt(dat_stats['phase_moments'][:,0,1])
    
    ax.errorbar(np.sqrt(dat_stats['phase_moments'][:,0,1]),dat_D_decorr['D_decorr'][:,0,0],yerr=dat_D_decorr['D_decorr'][:,0,1],color=c,fmt=marker_array[j],markersize=5)
    
    if not j:
      ax.plot([0],[0],color=c,marker='o',label=r'$\displaystyle K=%d$'%K_array[i])
  
  ax.plot([0],[0],color=c,marker=marker_array[j],label=r'$\displaystyle J_0=%d$'%J_array[j])
  

#print D_decorr
#print phase_var

popt,pcov = curve_fit(func,phase_var,D_decorr,sigma=D_decorr_err)
perr = np.sqrt(np.diag(pcov))

print popt
print pcov
print "error: ", perr

#m, b, r_value, p_value, std_err = stats.linregress(phase_var,D_decorr)

#m = np.sqrt(2000)
#b=0

#m = ret[0][0]
#b = ret[0][1]
print r'slope: $\displaystyle %5.3g\pm %5.3g$' % (popt,perr)

pltArray = np.linspace(0,1,1000)
ax.plot(pltArray,popt*pltArray,'--',color='grey')


ax.set_xlabel(r'$\displaystyle \sigma_{\phi}$',fontsize=12)
ax.set_ylabel(r'$\displaystyle D_{\phi}(decorr)$',fontsize=12)
#ax.set_xscale('log')
ax.set_xlim(0.1,0.3)
ax.set_ylim(5,20)
ax.legend(loc=(-0.15,0.66),prop={'size':10},ncol=2,numpoints=1)
plt.savefig('./pics/ERn/distance_phasevar_corr.pdf')
plt.show(block=False)