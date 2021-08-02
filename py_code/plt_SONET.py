import matplotlib.lines as mlines
import matplotlib.pyplot as plt

N = int(raw_input("N: "))
K= int(raw_input("K: "))
rateWnt= float(raw_input("rate: "))

cluster = int(raw_input("cluster: "))
net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[0]}
cluster_dict = {'scratch':cluster,'par':0,'q':0}

ana = {}
ana['recip'] = scriptFrame('r','eft',9,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict,mode_return=True)
ana['conv'] = scriptFrame('r','eft',11,alpha=[[0],[0,1],[0],[0]],netConst=net_dict,cluster=cluster_dict,mode_return=True)
ana['div'] = scriptFrame('r','eft',13,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict,mode_return=True)
ana['chain'] = scriptFrame('r','eft',11,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict,mode_return=True)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

mpl.rcParams['font.size'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12

plt.figure(figsize=(3,2.5))
ax = plt.axes([0.2,0.2,0.75,0.75])
#print ana['recip'].keys()

#print ana['recip']
mask = np.invert(np.isnan(ana['recip']['eft'][0,:,0]))
plt_bootstrap(np.linspace(-1,7,9),ana['recip']['eft'][0,:],ax,'k','-',ms=None,label=r'$\displaystyle \alpha_{recip}$',mask=mask)

mask = np.invert(np.isnan(ana['conv']['eft'][0,:,0]))
plt_bootstrap(np.linspace(0,1,11),ana['conv']['eft'][0,:],ax,'r','-',ms=None,label=r'$\displaystyle \alpha_{conv}$',mask=mask)

mask = np.invert(np.isnan(ana['div']['eft'][0,:,0]))
plt_bootstrap(np.linspace(0,6,13),ana['div']['eft'][0,:],ax,'b','-',ms=None,label=r'$\displaystyle \alpha_{div}$',mask=mask)

mask = np.invert(np.isnan(ana['chain']['eft'][0,:,0]))
plt_bootstrap(np.linspace(-0.5,0.5,11),ana['chain']['eft'][0,:],ax,'y','-',ms=None,label=r'$\displaystyle \alpha_{chain}$',mask=mask)

ax.plot([-2,10],[ana['conv']['eft'][0,0,0],ana['conv']['eft'][0,0,0]],'--',color='lightgrey',linewidth=0.5)
ax.set_yscale('log')
#ax.set_ylim([max(10**(-3),10**(np.floor(np.log10(np.nanmin(ana['conv']['eft'][0,:,0]))))),10**(np.ceil(np.log10(np.nanmax(ana['div']['eft'][0,:,0]))))])
ax.set_ylim([10**(-2),1])

ax.set_xlim([-1,7])
#plt.ylim([0,0.1])

ax.set_xticks(np.arange(-1,8))
#plt.yticks(np.linspace(0,0.1,6))
ax.set_xlabel(r'$\displaystyle \alpha$',fontsize=14)
ax.set_ylabel(r'$\displaystyle \varepsilon_{FT}$',fontsize=14)
ax.legend(loc=2,prop={'size':12})

save_str = './pics/SONET/total_fluxtube_N=%d_K=%d_rate=%d.pdf'%(N,K,rateWnt)
plt.savefig(save_str)
print "Figure saved to %s." % save_str
plt.show()
