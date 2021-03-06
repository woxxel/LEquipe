

J_array = [-1,-2,-5]

plt.figure(figsize=(4,3))
  
ax = plt.axes([0.15,0.55,0.8,0.4])
ax1 = plt.axes([0.15,0.1,0.8,0.4])

    
for j in range(len(J_array)):
  dat_tmp = scriptFrame(mode,'samples',1,hPconst={'cplg':[J_array[j]],'rate':10**np.linspace(2,3,17)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},cluster=cluster_dict,mode_return=1)
  
  dat_T[j] = dat_tmp
  
    
  for i in range(steps):
    col = i/float(steps)
    c = (col,col,col)
    ax.plot(pltRange,analysis['tau_RC'][0][i::steps],'-',color=c)
    ax1.plot(pltRange,analysis['T_min'][0][i::steps])
    
ax.set_xscale('log')
ax1.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([10**(-3),10**(2)])
ax.set_xlim([10**1.9,10**3.1])
ax1.set_xlim([10**1.9,10**3.1])

ax.set_xlabel(r'$\displaystyle \nu_p/Hz$',fontsize=14)
ax.set_ylabel(r'$\displaystyle \tau_{RC}/s$',fontsize=14)

  
if save:
  save_str = './pics/%s/J_nu_phase_N=%d.pdf' % (pltDat['topo'],meta_array[0]['N'])
  print 'figure saved as: %s ' % save_str
  plt.savefig(save_str)

plt.show(block=False)