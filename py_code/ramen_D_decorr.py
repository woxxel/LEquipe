
def ramen_cplg_influence(steps,scratch=2,save=0):
  
  N = 1000
  K = 100
  rateWnt = 1
  
  #cplg = -10
  
  cluster_dict = {'scratch':scratch,'par':0,'q':0}
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  plt_range = 10**np.linspace(0,3,steps)
  
  cplg_array = -np.array([1,2,5])
  
  plt.figure(figsize=(6,5))
  ax = plt.axes([0.12,0.12,0.77,0.77])
  ax2 = ax.twinx()
  
  #dat_decorr = scriptFrame('r','samples',steps,hPconst=['cplg',cplg],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict,mode_return=1)
  
  for i in range(len(cplg_array)):
    col = i/float(len(cplg_array))
    c = (col,col,col)
    
    J_p = cplg_array[i]
    
    dat_stat = scriptFrame('r','statistics',1,hPconst={'cplg':[J_p],'rate':plt_range},netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[0]},cluster={'scratch':2,'par':0,'q':1},mode_return=1)
    
    print dat_stat.keys()
  
    dat_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[J_p],'rate':plt_range},netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[2]},cluster={'scratch':2,'par':0,'q':1},mode_return=1)
    
    dat_T_conv = scriptFrame('r','samples',1,hPconst={'cplg':[J_p],'rate':plt_range},netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[30]},cluster={'scratch':2,'par':0,'q':1},mode_return=1)
    
    dat_eft = scriptFrame('r','eft',1,hPconst={'cplg':[J_p],'rate':plt_range},netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[0]},cluster={'scratch':2,'par':0,'q':1},mode_return=1)
    
    #print dat_decorr['D_decorr'][0][:,0]
    #print dat_eft['eft']
    #print dat_eft.keys()
    #print dat_T_conv
    mask_decorr = np.invert(np.isnan(dat_decorr['D_decorr_bs'][0,:,0]))
    mask_eft = np.invert(np.isnan(dat_eft['eft'][0,:,0]))
    #mask_T = np.invert(np.isnan(dat_T_conv['tau_RC'][0,:]))
    
    
    #mask_sync = np.where(dat_decorr['D_decorr'][0,:,0]<dat_eft['eft'][0,:,0])
    #print dat_eft['eft'][0,:,0]
    #factor = np.sqrt(-J_p)
    factor = 1
    #plt.errorbar(plt_range[mask]*factor,dat_decorr['D_decorr'][0][mask],dat_decorr['D_decorr_var'][0][mask],color=c,fmt='--',label='$\displaystyle J_p = %g$' % J_p)
    
    mask_reconv = np.isnan(dat_T_conv['tau_RC'][0])
    
    mask_inf = np.invert(dat_T_conv['tau_RC'][0] > 10**2)
    
    #for n in range(steps):	# don't plot solutions with 1 trajectory only
      #if np.isfinite(dat_T_conv['tau_RC'][0,n]):
      ##if dat_eft['eft'][0,n,0] > dat_decorr['D_decorr_bs'][0,n,0]:
	##dat_eft['eft'][0,n:,0] = dat_decorr['D_decorr_bs'][0,n,0]
	##dat_eft['eft'][0,n+1:,0] = np.nan #dat_decorr['D_decorr'][0,n,0]
	#dat_eft['eft'][0,n:,:] = np.nan
	#dat_decorr['D_decorr_bs'][0,n:,:] = np.nan
	#break
    
    plt_bootstrap(plt_range,dat_eft['eft'][0],ax,col=c,ls='-',ms=None,label=r'$\displaystyle J_p = %g$'%J_p,mask=mask_eft&mask_reconv&mask_inf)
    
    plt_bootstrap(plt_range,dat_decorr['D_decorr_bs'][0],ax,col=c,ls='--',ms=None,mask=mask_decorr&mask_reconv&mask_inf)
    #ax.plot(plt_range[mask_decorr]*factor,dat_decorr['D_decorr_bs'][0,:,0][mask_decorr],'--',color=c)
    
    
    #plt.errorbar(plt_range[mask_eft]*factor,dat_eft['eft'][0,:,0][mask_eft],yerr=dat_eft['eft'][0,:,1][mask_eft],color=c,fmt='-',linewidth=3,label=r'$\displaystyle J_p = %g$'%J_p)
    
    mask_stat = np.invert(np.isnan(dat_stat['D_spike_fail'][0]))
    print dat_stat['D_spike_fail']
    ax.plot(plt_range[mask_stat&mask_reconv],dat_stat['D_spike_fail'][0][mask_stat&mask_reconv],'-.',color=c,linewidth=2)
    
    print dat_eft['eft'][0]
    ax.plot(plt_range,dat_T_conv['tau_RC'][0],':',color=c,linewidth=2)
    
  
  #plt.plot(plt_range,D_decorr_ER*np.ones(len(plt_range)),'--r',linewidth=2,label=r'$\displaystyle D_{\phi}(ER)$')
  #plt.plot(plt_range,D_spike_fail*np.ones(len(plt_range)),'--r',linewidth=2,label='spike_failure')
  #ax.plot(plt_range,np.ones(len(plt_range)),'--r',linewidth=1,label=r'$\displaystyle D_{\phi}(spike)$')
  ax.plot([0,0],[-1,-1],'-k',label=r'$\displaystyle \varepsilon_{FT}$')
  ax.plot([0,0],[-1,-1],'--k',label=r'$\displaystyle D_{\phi}^{dc}$')
  ax.plot([0,0],[-1,-1],'-.k',label=r'$\displaystyle D_{\phi}^{sf}$')
  ax.plot([0,0],[-1,-1],':k',label=r'$\displaystyle \tau_{RC}$')
  
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim([10**0,10**3.1])
  ax.set_ylim([10**(-2),10**2])
  
  ax.set_xlabel(r'$\displaystyle \nu_p$ in Hz',fontsize=14)
  ax.set_ylabel(r'$\displaystyle D_{\phi}$',fontsize=14)
  #ax.legend(loc='center left',bbox_to_anchor=(1,0.5),prop={'size':12})
  ax.legend(loc=2,ncol=3,prop={'size':12},handlelength=2)
  
  ax2.set_yscale('log')
  ax2.set_ylim([10**(-2),10**2])
  ax2.set_ylabel(r'$\displaystyle \tau_{RC}$',fontsize=14)
  
  if save:
    save_str = './pics/poisson/overall.pdf'
    print 'figure saved as: %s ' % save_str
    plt.savefig(save_str)
  plt.show(block=False)
    
  #if mode=='D_decorr':
 
    #for i in range(3):
      
      #cplg = cplg_array[i]
      
      #col = i/float(len(cplg_array))
      #c = (col,col,col)
      
      #dat_tmp = scriptFrame('r','samples',steps,hPconst=['cplg',cplg],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict,mode_return=1)

      #Data['mean'][i] = dat_tmp['D_decorr_num'][0]
      #Data['var'][i] = dat_tmp['D_decorr_var'][0]
      
      #plt.errorbar(plt_range,Data['mean'][i],Data['var'][i],color=c,linewidth=2,ecolor=c,fmt='-o',label=r'$\displaystyle J_p = %g$'%(-10**i))
      
  #plt.plot(plt_range,D_decorr_ER*np.ones(len(plt_range)),'--r',linewidth=2,label='ER')
  
    #plt.xscale('log')
    #plt.ylim([0,10])
    #plt.xlabel(r'$\displaystyle \bar{\nu}_p$ in Hz',fontsize=18)
    #plt.ylabel(r'$\displaystyle D_{decorr}$',fontsize=18)
    #plt.legend(prop={'size':14},loc=3)
    #plt.show()
  
  #if mode=='eft':
    #for i in range(0,3):
      
      #col = i/2.2
      #c = (col,col,col)
      
      #dat_tmp = scriptFrame('r','eft',steps,hPconst=['cplg',-10**i],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2]},cluster=cluster_dict,mode_return=1)
      
      ##print dat_tmp['eft'][0]

      #Data['mean'][i] = dat_tmp['eft'][0]
      ##print dat_tmp['p']
      #for j in range(len(Data['mean'][i])):
	#if (Data['mean'][i][j] >= D_decorr_ER):# or np.isnan(Data['mean'][i][j]):
	  #Data['mean'][i][j] = D_decorr_ER
      
      #Data['var'][i] = dat_tmp['var_eft'][0]
      
      #mask = np.invert(np.isnan(Data['mean'][i]))
      #ax.errorbar(plt_range[mask],Data['mean'][i][mask],Data['var'][i][mask],color=c,linewidth=1,ecolor=c,fmt='-o',label=r'$\displaystyle J_p = %g$'%(-10**i))
      #ax.set_ylabel(r'$\displaystyle \varepsilon_{FT}$',fontsize=16)
      #ax.set_xlim([1,10**4])
    #dat_tmp = scriptFrame('r','eft',1,alpha=[[0],[0],[0],[0]],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2]},cluster=cluster_dict,mode_return=1)
    
    #ER_eft = D_decorr_ER#dat_tmp['eft'][0]
    ##ER_eft_var = dat_tmp['var_eft'][0]
    
    #ax.plot(plt_range,D_decorr_ER*np.ones(len(plt_range)),'--r',linewidth=2)#,label=r'$\displaystyle D_{decorr}$ from ER')
    #ax.plot(plt_range,D_spike_fail*np.ones(len(plt_range)),'--r',linewidth=2,label=r'$\displaystyle D_{spike fail}$')
    
    #ax.plot(plt_range,ER_eft*np.ones(len(plt_range)),'--k',linewidth=1,label=r'$\displaystyle \varepsilon_{FT}$ from ER')
    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_ylim([0,10])
    #ax.set_xlabel(r'$\displaystyle \bar{\nu}_p$ in Hz',fontsize=14)
    
    #ax.legend(prop={'size':12},loc=2)
    #plt.show()
  