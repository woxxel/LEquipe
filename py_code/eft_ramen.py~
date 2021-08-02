
# mode can be: (r)ead or (s)imulate
# mode_sim can be: 'eft','moments','synchrony','pert_spread'

def normal_eft_ramen(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}

  #scriptFrame(mode,'eft',1,alpha=[[0],[0],[0],[0]],netConst={'N':[5000,10000,20000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  scriptFrame(mode,'eft',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[100,200,500,1000,2000,4000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster={'scratch':1,'par':0,'q':q,'name':cluster})
  #scriptFrame(mode,'eft',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[1,2,5,10,20],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  scriptFrame(mode,'eft',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[100,200,500,1000,2000,4000],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster={'scratch':2,'par':0,'q':q,'name':cluster})
  
  
def normal_samples_ramen(mode,special,q,cluster):
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}

  #scriptFrame(mode,'samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[500,1000,2000,5000,10000,20000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[special]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[2000,5000,10000,20000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[special]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[special]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[100,200,500,1000,2000],'rateWnt':[0.2,0.5,1,2,5,10,20],'J0':[-1],'NeuronType':[2],'special':[special]},cluster=cluster_dict)
  scriptFrame(mode,'samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[100,500],'rateWnt':[1],'J0':[-0.5,-1,-2,-5,-10],'NeuronType':[2],'special':[special]},cluster=cluster_dict)


def normal_LE_ramen(mode,q,cluster):
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}

  scriptFrame(mode,'LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[500,1000,5000,10000,20000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,hide=0)
  scriptFrame(mode,'LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[100,200,500,1000,2000],'rateWnt':[0.2,0.5,1,2,5,10,20],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,hide=0)
  scriptFrame(mode,'LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[1,2,5,10,20],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,hide=0)
  
def SONET_eft_ramen(mode,N,K,rateWnt,q,cluster):

  net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[0]}
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}
  
  scriptFrame(mode,'eft',6,alpha=[[0],[0,1],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'eft',5,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'eft',7,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'eft',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict)

def SONET_eft_tmp_ramen(mode,N,K,rateWnt,q,cluster):

  net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[0]}
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}
  
  scriptFrame(mode,'eft_tmp',1,alpha=[[0],[0.5],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  #scriptFrame(mode,'eft',5,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  #scriptFrame(mode,'eft',7,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict)
  #scriptFrame(mode,'eft',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict)
  

def SONET_samples(mode,N,K,rateWnt,special,q,cluster):

  net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[special]}
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}
  
  #if N > 1000:
    #conv_steps=6
    #conv_max = 0.5
  #else:
  conv_steps = 6
  conv_max = 1
  
  #scriptFrame(mode,'samples',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict)
  #scriptFrame(mode,'samples',6,alpha=[[0],[0,1],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  #scriptFrame(mode,'samples',5,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'samples',7,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict)
  

def SONET_LE(mode,N,K,rateWnt,special,q,cluster):

  net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[special]}
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}
  
  #scriptFrame(mode,'LE',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict,hide=0)
  scriptFrame(mode,'LE',6,alpha=[[0],[0,1],[0],[0]],netConst=net_dict,cluster=cluster_dict,hide=0)
  #scriptFrame(mode,'LE',5,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict,hide=0)
  #scriptFrame(mode,'LE',7,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict,hide=0)
  
  

def SONET_statistics(mode,N,K,rateWnt,special,q,cluster):

  net_dict = {'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[-1],'NeuronType':[2],'special':[special]}
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':cluster}
  
  scriptFrame(mode,'statistics',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'statistics',6,alpha=[[0],[0,1],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'statistics',5,alpha=[[-1,7],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
  scriptFrame(mode,'statistics',7,alpha=[[0],[0],[0,6],[0]],netConst=net_dict,cluster=cluster_dict)
  
  
  
def poisson_ramen(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}
  
  J_array = -np.linspace(1,20,20)#[-1,-2,-3,-5,-7,-10,-15,-20]
  nu_array = 1000*np.array([0.1,1,2])
  
  scriptFrame(mode,'eft',1,hPconst={'cplg':np.linspace(-2,-20,10),'rate':[nu_array[0]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  scriptFrame(mode,'eft',1,hPconst={'cplg':np.linspace(-1,-10,10),'rate':[nu_array[1]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  scriptFrame(mode,'eft',1,hPconst={'cplg':np.linspace(-1,-6,6),'rate':[nu_array[2]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[2000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[2],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'eft',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[5],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)

def poisson_eft_tmp_ramen(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}
  
  scriptFrame(mode,'eft_tmp',1,hPconst={'cplg':[-2],'rate':[10**2.5]},hPcorr=2,netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  
def poisson_LE(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}
  
  scriptFrame(mode,'LE',1,hPconst={'cplg':[-1,-5],'rate':10**np.linspace(0,3,7)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,hide=0)


def poisson_samples_ramen(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}
  
  J_array = -np.linspace(1,20,20)#[-1,-2,-3,-5,-7,-10,-15,-20]
  nu_array = 1000*np.array([0.1,1,2])
  #scriptFrame(mode,'samples',1,hPconst={'cplg':J_array,'rate':nu_array},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-20],'rate':[nu_array[0]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  
  scriptFrame(mode,'samples',1,hPconst={'cplg':np.linspace(-14,-10,3),'rate':[nu_array[1]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  scriptFrame(mode,'samples',1,hPconst={'cplg':np.linspace(-20,-16,3),'rate':[nu_array[1]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=1,cluster=cluster_dict)
  
  scriptFrame(mode,'samples',1,hPconst={'cplg':np.linspace(-8,-4,3),'rate':[nu_array[2]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  scriptFrame(mode,'samples',1,hPconst={'cplg':np.linspace(-20,-10,6),'rate':[nu_array[2]]},hPcorr=10,netConst={'N':[2000],'K':[1000],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=1,cluster=cluster_dict)
  
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[2000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)

  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[400],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[100],'rateWnt':[2],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[100],'rateWnt':[5],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,13)},netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-5],'rate':10**np.linspace(3,2.5,9)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=0.5,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-5],'rate':10**np.linspace(2.5,2.25,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=1,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-5],'rate':10**np.linspace(2.25,2,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-2],'rate':10**np.linspace(3,2.75,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=0.5,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-2],'rate':10**np.linspace(2.75,2.5,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=1,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-2],'rate':10**np.linspace(2.5,2,9)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(3,2.75,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=1,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(2.75,2.5,5)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=2,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(2.5,2,9)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[30]},TC=10,cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},TC=0.2,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[500],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},TC=0.2,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[500],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},TC=0.2,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},TC=0.2,cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)

  ## special == 0, convergence rate
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,25)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[2],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[5],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1],'rate':10**np.linspace(0,3,13)},netConst={'N':[2000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)

  
  ### all done!!!
  #different K
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.2,-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(3,0,25)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[400],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  ##different rate
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[2],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[5],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  ##different N
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-0.5,-1,-2,-5,-10,-20],'rate':10**np.linspace(0,3,13)},netConst={'N':[2000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[2]},cluster=cluster_dict)
  
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(2,3,9)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(2,3,9)},netConst={'N':[500],'K':[50],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(2,3,9)},netConst={'N':[500],'K':[200],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict)
  
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(3,2,9)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict)
  #scriptFrame(mode,'samples',1,hPconst={'cplg':[-1,-2,-5],'rate':10**np.linspace(3,2,9)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[3]},cluster=cluster_dict)
  

def poisson_statistics_ramen(mode,q,cluster):
  cluster_dict = {'scratch':2,'par':0,'q':q,'name':cluster}
  
  #scriptFrame(mode,'statistics',1,hPconst={'cplg':[-1,-5],'rate':10**np.linspace(3,0,7)},netConst={'N':[500],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  scriptFrame(mode,'statistics',1,hPconst={'cplg':[-1,-5],'rate':10**np.linspace(3,0,7)},netConst={'N':[1000],'K':[100],'rateWnt':[1],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  
  
execfile('Code/py_code/scriptRahmen.py')


cluster = sys.argv[1]
q = int(sys.argv[2])
mode = sys.argv[3]
if len(sys.argv)>4:
  comment = sys.argv[4]

poisson_samples_ramen(mode,q=q,cluster=cluster)
#poisson_ramen(mode,q=q,cluster=cluster)
#poisson_eft_tmp_ramen(mode,q=q,cluster=cluster)
#normal_samples_ramen(mode,special=0,q=q,cluster=cluster)
#normal_eft_ramen(mode,q=1,cluster=cluster)

#normal_LE_ramen(mode,q=q,cluster=cluster)

#SONET_LE('s',N=1000,K=100,rateWnt=1,special=2,q=q,cluster=cluster)
#SONET_samples('s',N=1000,K=100,rateWnt=1,special=8,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=1000,rateWnt=1,special=8,q=q,cluster=cluster)
#SONET_samples('s',N=1000,K=100,rateWnt=1,special=0,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=1000,rateWnt=1,special=2,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=1000,rateWnt=1,special=0,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=1000,rateWnt=10,special=2,q=q,cluster=cluster)
##SONET_samples('s',N=10000,K=1000,rateWnt=10,special=0,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=500,rateWnt=1,special=2,q=q,cluster=cluster)
##SONET_samples('s',N=10000,K=500,rateWnt=1,special=0,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=2000,rateWnt=1,special=2,q=q,cluster=cluster)
#SONET_samples('s',N=10000,K=2000,rateWnt=1,special=0,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=1000,K=100,rateWnt=1,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=1000,K=100,rateWnt=10,q=q,cluster=cluster)
#SONET_eft_ramen(mode,N=10000,K=1000,rateWnt=1,q=q,cluster=cluster)
#SONET_eft_tmp_ramen(mode,N=10000,K=1000,rateWnt=1,q=q,cluster=cluster)
#SONET_eft_ramen(mode,N=1000,K=100,rateWnt=10,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=10000,K=1000,rateWnt=10,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=10000,K=100,rateWnt=1,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=10000,K=500,rateWnt=1,q=q,cluster=cluster)
#SONET_eft_ramen('s',N=10000,K=2000,rateWnt=1,q=q,cluster=cluster)
#SONET_statistics('s',N=1000,K=100,rateWnt=10,special=0,q=q,cluster=cluster)
#SONET_statistics('s',N=10000,K=1000,rateWnt=10,special=0,q=q,cluster=cluster)

#poisson_ramen(mode,q=q,cluster=cluster)
#poisson_LE(mode,q=q,cluster=cluster)


#poisson_statistics_ramen(mode,q=q,cluster=cluster)

