def ramen_statistics(mode,q=1):
  reread = 0
  save = 1
  scratch = 0
  plot = 0
  network_dict = {'N':[1000],'K':[100],'rateWnt':[10],'J0':[-1],'NeuronType':[2]}
  cluster_dict = {'scratch':scratch,'par':0,'q':q}
  
  scriptFrame(mode,'statistics',5,alpha=[[-1,7],[0],[0],[0]],netConst=network_dict,cluster=cluster_dict,reread=reread,save=save,plot=plot)
  scriptFrame(mode,'statistics',6,alpha=[[0],[0,2],[0],[0]],netConst=network_dict,cluster=cluster_dict,reread=reread,save=save,plot=plot)
  scriptFrame(mode,'statistics',7,alpha=[[0],[0],[0,6],[0]],netConst=network_dict,cluster=cluster_dict,reread=reread,save=save,plot=plot)
  scriptFrame(mode,'statistics',6,alpha=[[0],[0.5],[1],[-0.5,0.5]],netConst=network_dict,cluster=cluster_dict,reread=reread,save=save,plot=plot)
  
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_rate',constVal=10,mode=mode)
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_rate',constVal=100,mode=mode)
  
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_cplg',constVal=10,mode=mode)
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_cplg',constVal=100,mode=mode)
  
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_variance',constVal=10,mode=mode)
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_variance',constVal=100,mode=mode)
  
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_current',constVal=10,mode=mode)
  #scriptPoisson(4,[0,3],['statistics'],N=200,hPconst='hP_current',constVal=100,mode=mode)