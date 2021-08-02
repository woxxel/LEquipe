#execfile('Code/py_code/scriptRahmen.py')
def ramen_samples_ER(q,clustername='skadi',mode=None):
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':clustername}
  
  if mode == 'r':
    plot = 0
    save = 1
  else:
    plot = 0
    save = 0
  N_array = np.array([1000,2000,5000,10000,20000])#,20000])
  
  for i in range(len(N_array)):
    if (mode == 's') or not mode:
      scriptFrame('s','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[N_array[i]],'K':[1000],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  if (mode == 'r') or not mode:
    scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':N_array,'K':[1000],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  if (mode == 's'):
    scriptFrame('s','LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':N_array,'K':[1000],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  
  K_array = np.array([20,50,100,200,500,1000,2000,4000])
  
  for i in range(len(K_array)):
    if (mode == 's') or not mode:
      scriptFrame('s','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[K_array[i]],'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  if (mode == 'r') or not mode:
    scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':K_array,'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  #scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':K_array,'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  if (mode == 's'):
    scriptFrame('s','LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':K_array,'rateWnt':[10],'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  
  nu_array = np.array([0.5,1,2,3,4,5,10,20,50])
  
  for i in range(len(nu_array)):
    if (mode == 's') or not mode:
      scriptFrame('s','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':nu_array,'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  if (mode == 'r') or not mode:
    scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':nu_array,'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  if (mode == 's'):
    scriptFrame('s','LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':nu_array,'J0':[-1],'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  
  J_array = np.array([-0.1,-0.2,-0.5,-1,-2,-5])
  
  for i in range(len(J_array)):
    if (mode == 's') or not mode:
      scriptFrame('s','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[10],'J0':J_array,'NeuronType':[2],'special':[0]},cluster=cluster_dict)
  if (mode == 'r') or not mode:
    scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[10],'J0':J_array,'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)
  if (mode == 's'):
    scriptFrame('s','LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[10000],'K':[1000],'rateWnt':[10],'J0':J_array,'NeuronType':[2],'special':[0]},cluster=cluster_dict,plot=plot,save=save)



def ramen_samples_SONET(q):
  net_dict = {'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]}
  cluster_dict = {'scratch':1,'par':0,'q':q}
  
  recip = np.linspace(-1,7,5)
  conv = np.linspace(0,2,6)
  div = np.linspace(0,6,7)
  chain = np.linspace(-0.5,0.5,6)
  
  for i in conv:
    try:
      scriptFrame('s','samples',1,alpha=[[0],[i],[0],[0]],netConst=net_dict,cluster=cluster_dict)
      scriptFrame('r','samples',1,alpha=[[0],[i],[0],[0]],netConst=net_dict,cluster=cluster_dict)
    except:
      print "something went wrong"
  
  for i in div:
    try:
      scriptFrame('s','samples',1,alpha=[[0],[0],[i],[0]],netConst=net_dict,cluster=cluster_dict)
      scriptFrame('r','samples',1,alpha=[[0],[0],[i],[0]],netConst=net_dict,cluster=cluster_dict)
    except:
      print "something went wrong"
  
  for i in chain:
    try:
      scriptFrame('s','samples',1,alpha=[[0],[0.5],[1],[i]],netConst=net_dict,cluster=cluster_dict)
      scriptFrame('r','samples',1,alpha=[[0],[0.5],[1],[i]],netConst=net_dict,cluster=cluster_dict)  
    except:
      print "something went wrong"
  
  for i in recip:
    try:
      scriptFrame('s','samples',1,alpha=[[i],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
      scriptFrame('r','samples',1,alpha=[[i],[0],[0],[0]],netConst=net_dict,cluster=cluster_dict)
    except:
      print "something went wrong"
  

def ramen_samples_drive(q=1,name='frigg'):
  reread = 0
  save=0
  cluster_dict = {'scratch':1,'par':0,'q':q,'name':name}
  
  #scriptFrame('s','samples',9,hPconst=['cplg',-10],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_min=10,bound_max=5000)
  #scriptFrame('r','samples',9,hPconst=['cplg',-10],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_min=10,bound_max=5000)
  
  scriptFrame('s','samples',9,hPconst=['cplg',-100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_min=10,bound_max=1000)
  scriptFrame('r','samples',9,hPconst=['cplg',-100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_min=10,bound_max=1000)
  
  scriptFrame('s','samples',9,hPconst=['cplg',-1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=1000)
  scriptFrame('r','samples',9,hPconst=['cplg',-1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=1000)
  
  scriptFrame('s','samples',4,hPconst=['rate',10],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=10000)
  scriptFrame('r','samples',4,hPconst=['rate',10],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=10000)
  
  scriptFrame('s','samples',4,hPconst=['rate',100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=500)
  scriptFrame('r','samples',4,hPconst=['rate',100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=500)
  
  scriptFrame('s','samples',4,hPconst=['rate',1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=500)
  scriptFrame('r','samples',4,hPconst=['rate',1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'NeuronType':[2],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save,bound_max=500)
  
  #scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[500,1000,2000,5000],'K':[100],'rateWnt':[10],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save)
  #scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[20,50,100,200,500],'rateWnt':[10],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save)
  #scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[100],'rateWnt':[1,2,5,10,20],'J0':[-1],'special':[0]},cluster=cluster_dict,reread=reread,save=save)
  #scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-0.1,-1,-2,-5,-10],'special':[0]},cluster=cluster_dict,reread=reread,save=save)
  
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[10,20,50,100,200],'rateWnt':[2],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[10,20,50,100,200],'rateWnt':[10],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[10,20,50,100,200],'rateWnt':[20],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[20],'rateWnt':[1,2,5,10,20],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[50],'rateWnt':[1,2,5,10,20],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[100],'rateWnt':[1,2,5,10,20],'J0':[-1]},reread=reread,save=save)
  ##scriptFrame(mode,'samples',1,hPconst={'rate':[10],'cplg':[-10]},netConst={'N':[1000],'K':[200],'rateWnt':[1,2,5,10,20],'J0':[-1]},reread=reread,save=save)
  
  
def ramen_eft_drive(mode,q=1):
  
  steps = 13
  scriptFrame(mode,'eft',steps,hPconst=['cplg',-10],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['cplg',-100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['cplg',-1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  
  scriptFrame(mode,'eft',steps,hPconst=['current',100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['current',1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['current',10000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  
  scriptFrame(mode,'eft',steps,hPconst=['variance',100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['variance',1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['variance',10000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  
  scriptFrame(mode,'eft',steps,hPconst=['rate',100],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['rate',1000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})
  scriptFrame(mode,'eft',steps,hPconst=['rate',10000],netConst={'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]},cluster={'scratch':1,'par':0,'q':q})