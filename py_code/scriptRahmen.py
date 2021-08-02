import numpy as np
import pylab as plt
#import matplotlib.pyplot as plt
import os, imp
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

imp.load_source('fts_simu', 'Code/py_code/tubesize.py')
imp.load_source('fit_func', 'Code/py_code/fit_functions.py')
from fts_simu import *
from fit_func import *


def scriptPertSpread(mode=None,mode_sim=['pert_spread'],alpha=None,hPvars=None,TC=1,N=1000):
  if not mode:
    mode = raw_input('(r)ead or (s)imulate data? ')
  
  str_alpha_label = ['recip','conv','div','chain']
  str_alpha = ['']*4
  
  for i in range(4):
    if alpha:
      str_alpha[i] = ('alpha_%s=%4.2f' % (str_alpha_label[i],alpha[i])).rstrip('0').rstrip('.') + '_'
    else:
      str_alpha[i] = ('alpha_%s=%4.2f' % (str_alpha_label[i],0)).rstrip('0').rstrip('.') + '_'
  
  str_hP_label = ['hPrate','hPcplg']
  str_hP = ['']*2
  
  for i in range(2):
    if alpha:
      str_hP[i] = ('%s=%4.2f' % (str_hP_label[i],0)).rstrip('0').rstrip('.') + '_'
    else:
      str_hP[i] = ('%s=%4.2f' % (str_hP_label[i],hPvars[i])).rstrip('0').rstrip('.') + '_'

    
  #do it with phases or at spiketimes? or both?
  if mode == 's':
    if 'plot' in mode_sim:
      steps = {'ps':1,'co':1,'in':1,'pd':1}
    
    # obtain average flux tube radius from simulation
    for dat in os.listdir('data/info/'):
      if ((str_alpha[0] in dat) and (str_alpha[1] in dat) and (str_alpha[2] in dat) and (str_alpha[3] in dat) and ('SONET_para' in dat)):
	eft,vareft,numsim = analyze(['eft'],fileName_in='data/info/'+dat,scriptName_in='SONET',call=1)
	if eft > 0.1:
	  eft = 0.1
	eft *= 5
    
    # then simulate diverging trajetory
    if 'plot' in mode_sim:
      diverged = 0
      while not diverged:
	if alpha:
	  net = network(mode='t',topo='S',alpha=alpha,N=N)
	if hPvars:
	  net = network(mode='t',topo='p',hPrate=hPvars[0],hPcplg=hPvars[1],N=N)
	sim = simulation_data(net,mode_sim=mode_sim,steps=steps,TC=TC)
	break_it = net.setup(sim,check=0)
	if break_it == 2:
	  diverged = 1
	elif break_it == 1:
	  erase(net.Path['script'],net.Path['info_file'],net.Hash)
	  diverged = 1
	elif break_it == 0:
	  net.simulation(sim,pert_spread_sz=eft)
	  break
	  for dat in os.listdir(net.Path['results']):	# Reference- and perturbed trajectory
	    if 'DataOut-' in dat:
	      Data = readDataOut(net.Path['results'] + dat)
	    elif 'RefTraj-' in dat:
	      DataRef = readDataOut(net.Path['results'] + dat)
	  if np.linalg.norm(DataRef['measurePhases'][-1]-Data['measurePhases'][-1]) > 1:
	    diverged = 1
	    print 'trajectory diverged'
	  else:
	    erase(net.Path['script'],net.Path['info_file'],net.Hash)
  
  else:
    sessionPath = 'data/'
    if alpha:
      scriptName = 'SONET_'
    else:
      scriptName = 'poisson_'
    for mode_iter in mode_sim:
      scriptName += mode_iter + '_'
    scriptName += 'para_'
    
    found = 0
    for dat in os.listdir(sessionPath + '/info'):
      #print dat
      if alpha:
	if ((str_alpha[0] in dat) and (str_alpha[1] in dat) and (str_alpha[2] in dat) and (str_alpha[3] in dat) and (scriptName in dat)):
	#print dat
	  pert_data = analyze(mode_sim,fileName_in=sessionPath + 'info/' + dat,scriptName_in=scriptName)
	  break
      else:
	if ((str_hP[0] in dat) and (str_hP[1] in dat) and (scriptName in dat)):
	  pert_data = analyze(mode_sim,fileName_in=sessionPath + 'info/' + dat,scriptName_in=scriptName)
	  break
    
    #print 'Numbers of unsuccessful switches before divergence: %5.3f' % np.mean(pert_data[0])
    #print 'Average time of divergence initiation: %5.3f' % pert_data[1]
    #print 'Average velocity of divergence: %5.3f' % pert_data[2]
    
    return pert_data
    
    
def scriptFrame(mode='r',mode_sim='eft',steps=1,alpha=None,hPconst=None,hPcorr=1,TC=0.5,netConst=None,cluster=None,perturbation=None,save=0,hide=1,plot=0,reread=0,mode_return=False,bound_min=None,bound_max=None):
  
  if cluster:
    if cluster['scratch']:
      clusterstr = '/scratch%02d/aschmidt/' % cluster['scratch']
    else:
      clusterstr = ''
  else:
    clusterstr = ''
    cluster = {'scratch':0,'par':0,'q':0}
    
  infoDir = clusterstr + 'data/info/'
  
  meta_steps = 1
  assert type(netConst) == dict, 'Please provide variables of hPconst in form of a dictionary!'
  assert len(netConst.keys()) >= 5, 'Please provide values for all 5 (6) network parameters! (N,K,J0,rateWnt,NeuronType(,special))'
  
  sim_range = None
  
  str_dict = {}
  str_dict['title'] = 'network constants:\t'
  str_dict['plot'] = ''
  # set up the meta iteration (over network parameters)
  
  meta_lists = []
  meta_steps_array = []
  for key in netConst.keys():
    len_tmp = len(netConst[key])
    if len_tmp > 1:
      meta_lists.append(key)
      meta_steps *= len_tmp
      meta_steps_array.append(len_tmp)
      sim_range = netConst[key]
      if key == 'rateWnt':
	str_dict['plot'] = r'$\displaystyle \bar{\nu}$'
      else:
	str_dict['plot'] = key
      
      str_dict['save'] = key
    else:
      str_dict['title'] += ' %s = %g,' % (key,netConst[key][0])
  
  str_dict['title'] = str_dict['title'].rstrip(',')
  
  if len(meta_lists) == 2:
    #print meta_steps_array
    #print meta_steps
    meta_values = np.zeros((meta_steps,2))
    meta_values[:,0] = np.kron(netConst[meta_lists[0]],np.ones(meta_steps_array[1]))
    meta_values[:,1] = np.kron(np.ones(meta_steps_array[0]),netConst[meta_lists[1]])
  
  meta_array = []
  for meta_step in range(meta_steps):
    meta_array.append({})
    for key in netConst.keys():
      
      if len(meta_lists) == 1:
	if key in meta_lists:#len(netConst[key])>1:
	  meta_array[meta_step][key] = netConst[key][meta_step]
	else:
	  meta_array[meta_step][key] = netConst[key][0]
      
      elif len(meta_lists) == 2:
	if key in meta_lists:#len(netConst[key])>1:
	  meta_array[meta_step][key] = meta_values[meta_step,np.where([key==i for i in meta_lists])[0][0]]	#netConst[key][meta_step%len(netConst[key])]
	else:
	  meta_array[meta_step][key] = netConst[key][0]
      else:
	meta_array[meta_step][key] = netConst[key][0]
  
  
  # initiate the arrays for iteration of topo-params and drive-params
  if hPconst:
    steps = len(hPconst['rate'])*len(hPconst['cplg'])
  
  alpha_sim = np.zeros((steps,4))
  alpha_sim[:] = None
  
  hP_sim = np.zeros((steps,2))
  hP_sim[:] = None
    
  if alpha:
    simVars = 4
    assert type(alpha)==list, 'Please specify the second order parameters in a list structure'
    assert len(alpha)==4, 'Please specify the Second order motif parameters properly! (recip,conv,div,chain)'
    
    scriptName = 'SONET'
    str_para_list = ['recip','conv','div','chain']
    
    # check, whether there should be an iteration over an alpha
    try:
      str_para = str_para_list[np.where(np.array([len(item) for item in alpha])>1)[0][0]]
      alpha_idx = np.where(np.array([len(item) for item in alpha])>1)[0][0]
    except:
      str_para = ''
      alpha_idx = 0
    
    alpha_sim = np.zeros((simVars,steps))
    
    for i in range(simVars):	# construct simulation ranges
      alpha_sim[i] = np.linspace(alpha[i][0],alpha[i][-1],steps)
    
    if steps > 1:
      sim_range = np.linspace(alpha_sim[alpha_idx][0],alpha_sim[alpha_idx][-1],steps)
      print sim_range
    if len(meta_array)==1:
      str_dict['title'] = r'$\displaystyle \alpha_{%s}$' % str_para
      str_dict['plot'] = r'$\displaystyle \alpha_{%s}$' % str_para
      str_dict['save'] = str_para
      str_dict['range'] = sim_range
    
    alpha_sim = alpha_sim.transpose()
    simulation_paras = alpha_sim
    
  elif hPconst:		# declare variance as "relative to balanced state recurrent variance" (=J^2*nu = 10)
    
    assert type(hPconst) == dict, 'Please provide variables of hPconst in form of a dictionary!'
    assert len(hPconst.keys()) == 2, 'Please specify the poisson parameters properly! (hPrate,hPcplg)' 
    
    scriptName = 'poisson'  
    simVars = 2
    str_para_list = ['rate','cplg']
          
    str_dict['title'] += '\n drive parameters:\t'
    
    simulation_paras = np.zeros((steps,simVars))
    for i in range(len(hPconst['rate'])):
      for j in range(len(hPconst['cplg'])):
	idx = i*len(hPconst['cplg'])+j
	simulation_paras[idx][:] = [hPconst['rate'][i],hPconst['cplg'][j]]
    
    hP_sim = np.array(simulation_paras)
    
    print str_dict['plot']
    str_dict['save'] = str_dict['plot']			#'nu_p=%3.1f' % hPconst[1]
    if str_dict['save'] == 'rate':
      str_dict['save'] += '_cplg=%g' % hPconst['cplg'][0]
    elif str_dict['save'] == 'rateWnt':
      str_dict['save'] += '_K=%g' % netConst['K'][0]
    elif str_dict['save'] == 'K':
      str_dict['save'] += '_rateWnt=%g' % netConst['rateWnt'][0]
    #str_dict['plot'] = 'plotting beta'			#r'$\displaystyle \log_{10}(-J_{p})$'
    
    str_dict['rateWnt'] = 'rateWnt=1 (naja)'
    str_dict['range'] = sim_range
    str_para = 'bla'
    
      
    #else:
      #str_para = ''
      #hP_sim, sim_range, str_dict = poisson_set_paras(hPconst,steps)	# construct simulation ranges
      #str_dict['const'] = hPconst
      #simulation_paras = hP_sim
    
  if not mode:
    mode = raw_input('(r)ead or (s)imulate data? ')
    
  eft = None
  
  str_mode = {'s':'\nNow simulating network with:','r':'\nNow reading simulation with:'}
  
  #do_it = {}
  #for i in range(len(hPconst['cplg'])):
    #do_it[hPconst['cplg'][i]] = 0
    #print hPconst['cplg'][i]
  #print "he"
  if mode == 'r':
    #print "steps: ", steps
    ##print meta_steps
    analysis = {}
    pltDat = initialize_plt(sim_range,steps,scriptName,str_para,scriptName[0])
    #print pltDat
    #if 'special' in netConst.keys():
      #if netConst['special'][0] in [2,3]:
    if hPconst:
      steps_array = [len(hPconst['rate']),len(hPconst['cplg'])]
      pltDat['steps'] = [len(hPconst['rate']),len(hPconst['cplg'])]
      pltDat['range'] = np.column_stack((np.kron(hPconst['rate'],np.ones(len(hPconst['cplg']))),np.kron(np.ones(len(hPconst['rate'])),hPconst['cplg'])))
  
  str_array_net = []
  
  for meta_step in range(meta_steps):
    
    str_array_net = []
    print "network parameters:\n"
    for key in meta_array[meta_step].keys():
      str_tmp = '%s=%g' % (key,meta_array[meta_step][key])
      print '\t%s' % str_tmp
      str_array_net.append(str_tmp + '_')
      
    #print str_array_net
    
    steps_iter = iter(np.arange(steps))
    
    for step in steps_iter:
      ### bounds if I want:
      if not bound_min:
	bound_min = -2
      if not bound_max:
	bound_max = 100000
      if simulation_paras[step][0] < bound_min or simulation_paras[step][0] > bound_max:
	print "Not now"
	continue
      
      simulated = 0
      
      str_array_para = []
      print "further parameters:\n"
      for i in range(simVars):
	str_array_para.append(('%s=%5.3f' % (str_para_list[i],simulation_paras[step][i])).rstrip('0').rstrip('.'))
	print '\t%s %s' % (scriptName,str_array_para[i])
	#print str_array_para[i]
	str_array_para[i] += '_'
      
	  
      if mode_sim in ['eft']:# or (mode_sim in ['mosaic','samples'] and not (mode == 'r')):
	str_script = 'SimInfo_' + scriptName + '_para'
	infoPath = infoDir + scriptName + '/'
      else:
	str_script = 'SimInfo_' + scriptName + '_' + mode_sim + '_para'
	infoPath = infoDir + scriptName + '_' + mode_sim + '/'
      
      if not os.path.exists(infoPath):
	os.mkdir(infoPath)
      
      #print str_array_para
      #print str_array_net
      #print infoPath
      
      for dat in os.listdir(infoPath):
	
	#print dat
	if (mode=='r') and (('special' in netConst.keys()) and not 'special' in dat):
	  continue
	if (mode=='r') and (not ('special' in netConst.keys()) and 'special' in dat):
	  continue
	if 'special' in netConst.keys():
	  compare = 6
	else:
	  compare = 5
	
	##print [str_test in dat for str_test in str_array_para]
	#print [str_test in dat for str_test in str_array_net]
	#print dat
	#if all([str_test in dat for str_test in str_array_para]):
	  #print dat
	if (str_script in dat) and all([str_test in dat for str_test in str_array_para]) and (np.sum([str_test in dat for str_test in str_array_net])==compare):
	  if mode == 's':
	    simulated = 1
	  else:
	    print "found: %s" % dat
	    net = network('r',infoPath + dat)
	    net.scriptName = scriptName
	    net.analysis = {}
	    sim = simulation_data('r',infoPath + dat,net,cluster=cluster,mode=mode_sim)
	    
	    break_it = net.setup(sim,hide=hide,call=1)
	    if not break_it:
	      
	      if sim.mode == 'samples':
		cluster_read = {'scratch':cluster['scratch'],'par':0,'q':0}
		if net.Const['special'] == 0:
		  if net.topo == 'p':
		    try:
		      ana_LE = scriptFrame('r','LE',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster={'scratch':2,'par':0,'q':0},mode_return=1)
		      net.analysis['1stLE'] = ana_LE['subLyapunovExponents'][0,0,:2]
		      print 'LE: ', net.analysis['1stLE']
		    except:
		      1
		  
		if net.Const['special'] == 1:
		  
		  print "\tnow reading D decorr..."
		  try:
		    ana_D_decorr = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
		    net.analysis['D_decorr_ER'] = ana_D_decorr['D_decorr_bs'][0,0]
		    print 'D: ', net.analysis['D_decorr_ER']
		  except:
		    1
		  print "\tnow reading tau conv..."
		  try:
		    ana_tau_conv = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
		    net.analysis['tau_conv_ER'] = ana_tau_conv['tau_conv_later'][0,0]
		    print 'tau: ', net.analysis['tau_conv_ER']
		  except:
		    1

		  print "\tnow reading LE..."
		  #ana_LE = scriptFrame('r','LE',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[1]},cluster=cluster_read,mode_return=1)
		  if net.topo == 'p':
		    try:
		      ana_LE = scriptFrame('r','LE',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster={'scratch':2,'par':0,'q':0},mode_return=1)
		      
		      net.analysis['1stLE'] = ana_LE['subLyapunovExponents'][0,0,:2]
		      print 'LE: ', net.analysis['1stLE']
		    except:
		      1
		  else:
		    try:
		      net.analysis['1stLE'] = ana_LE['LyapunovExponents'][0,0,1]
		      print 'LE: ', net.analysis['1stLE']
		    except:
		      1
		
		  print "\treading done!"
		    
		  
		  
		  
		
		if net.Const['special'] == 3:
		  if net.topo == 'p':
		    cluster_tmp = {'scratch':cluster['scratch'],'par':0,'q':0}
		    ana_D_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_tmp,mode_return=1)
		    #ana_tau_conv = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster=cluster_tmp,mode_return=1)
		  
		  else:
		    ana_D_decorr = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_read,mode_return=1)
		    ana_tau_conv = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster=cluster_read,mode_return=1)
		  
		  net.analysis['D_decorr'] = ana_D_decorr['D_decorr_bs'][0,0]
		  #net.analysis['tau_conv'] = ana_tau_conv['tau_conv'][0,0]
		  
		if net.Const['special'] == 4:
		  cluster_tmp = {'scratch':cluster['scratch'],'par':0,'q':0}
		  print "\tnow reading D decorr..."
		  try:
		    if not ('D_decorr_ER' in net.analysis.keys()):
		      ana_D_decorr = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
		      net.analysis['D_decorr_ER'] = ana_D_decorr['D_decorr_bs'][0,0]
		      print 'D: ', net.analysis['D_decorr_ER']
		  except:
		    1
		  
		  try:
		    if not ('D_decorr_poisson' in net.analysis.keys()):
		      ana2_D_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_tmp,mode_return=1)
		      net.analysis['D_decorr_poisson'] = ana2_D_decorr['D_decorr_bs'][0,0]
		      print 'Dp: ', net.analysis['D_decorr_poisson']
		  except:
		    1
		    
		  try:
		    if not ('tau_conv' in net.analysis.keys()):
		      ana_tau_conv = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster=cluster_tmp,mode_return=1)
		      
		      #print ana_tau_conv.keys()
		      #print ana_tau_conv
		      net.analysis['tau_conv'] = ana_tau_conv['tau_conv_first'][0,0]
		      #net.analysis['D_decorr_poisson'] = ana2_D_decorr['D_decorr_bs'][0,0]
		      print 'tau p: ', net.analysis['tau_conv']
		  except:
		    1
		  #try:
		    #print analysis.keys()
		  #if not ('phi_corr_ER' in net.analysis.keys()):
		    #ana_phi = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[4]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
		    #net.analysis['phi_corr_ER'] = ana_phi['phi_corr']
		    #print 'phi_corr:', net.analysis['phi_corr_ER']
		  #except:
		    #1
		    
		  #net.analysis = np.copy(analysis_tmp)
		  #analysis = analysis_tmp
		#if analysis:
		  ##for key in analysis_tmp.keys():
		    #net.analysis[key] = analysis_tmp[key]
	      #try:
	      analysis_tmp = analyze(mode_sim,net,sim,plot=plot,call=1,save=save,reread=reread)
	      #print 'D: ', net.analysis['D_decorr_ER']
	      #print analysis_tmp
	      #print analysis

	      analysis_handover(analysis_tmp,analysis,step,meta_step,steps,meta_steps)
	      #except:
		#1
	    break
      
      if mode == 's' and not simulated:
	

	#print 'he'
	##
	#if mode_sim == 'samples':
	  #if netConst['special'][0] == 3:
	    #print "reading..."
	    #print meta_array[meta_step]
	    #res_D_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster={'scratch':2,'par':0,'q':0},mode_return=1)
	    #print 'ho'
	    #res_eft = scriptFrame('r','eft',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster={'scratch':2,'par':0,'q':0},mode_return=1)
	    #print 'he'
	    #print res_D_decorr
	    #print res_eft['eft'][0,0,0]
	    #if 10*res_eft['eft'][0,0,0] > res_D_decorr['D_decorr'][0,0,0]:
	      #do_it[hP_sim[step][1]] = 1
	    #if not do_it[hP_sim[step][1]]:
	      #try:
		#if res_D_decorr['D_decorr'][0,0,0] > 10*res_eft['eft'][0,0,0]:
		  #print "no reconverging trajectories to be found here"
		  #continue
	      #except:
		#print "One of the simulations is not yet finished"
		#continue
	    #return
	net = network(alpha=alpha_sim[step],hPVars=hP_sim[step],hPcorr=hPcorr,netConst=meta_array[meta_step])
	sim = simulation_data(net=net,cluster=cluster,mode=mode_sim,TC=TC)
	
	if sim.mode == 'samples':
	  cluster_read = {'scratch':cluster['scratch'],'par':0,'q':0}
	  if net.Const['special'] == 3:
	    ana_D_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_read,mode_return=1)
	    
	    #if (ana_D_decorr['D_decorr_bs'][0,0,0] > 3) and (ana_D_decorr['D_decorr_bs'][0,0,0] < 10):
	      #continue
	    
	#if (mode_sim in ['mosaic','samples']) and (not perturbation):
	  #try:
	    #if data_tmp['eft'] > 1:
	      #perturbation = 0.1
	    #elif mode_sim == 'mosaic':
	      #perturbation = 4*data_tmp['eft']	# sample 4 times the size of the fluxtube
	    #elif mode_sim == 'samples':
	      #perturbation = data_tmp['eft']
	    #assert perturbation, "No perturbation could be set from data!"
	    #assert not np.isnan(perturbation), "No perturbation could be set from data!"
	    #print "Perturbation set to %g" % perturbation
	  #except:
	perturbation = 0.02
	    #print "No fitting simulation could be found, maximum perturbation size set to %g!" % perturbation
	if sim.mode == 'samples':
	  net.analysis = {}
	  if net.Const['special'] == 0:
	    cluster_read = {'scratch':cluster['scratch'],'par':0,'q':0}
	    anaTmp = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_read,mode_return=1)
	  
	    for key in anaTmp.keys():
	      net.analysis[key] = anaTmp[key][0,0]
	  
	if sim.mode == 'eft':
	  net.analysis = {}
	  
	  cluster_read = {'scratch':cluster['scratch'],'par':0,'q':0}
	  
	  if net.topo == 'p':
	    ana_D_decorr = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster=cluster_read,mode_return=1)
	    ana_tau_conv = scriptFrame('r','samples',1,hPconst={'cplg':[hP_sim[step][1]],'rate':[hP_sim[step][0]]},netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster=cluster_read,mode_return=1)
	  else:
	    ana_D_decorr = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[2]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
	    ana_tau_conv = scriptFrame('r','samples',1,alpha=[[0],[0],[0],[0]],netConst={'N':[meta_array[meta_step]['N']],'K':[meta_array[meta_step]['K']],'rateWnt':[meta_array[meta_step]['rateWnt']],'J0':[meta_array[meta_step]['J0']],'NeuronType':[meta_array[meta_step]['NeuronType']],'special':[0]},cluster={'scratch':1,'par':0,'q':0},mode_return=1)
	  
	  if 'D_decorr' in ana_D_decorr.keys():
	    net.analysis['D_decorr'] = ana_D_decorr['D_decorr'][0,0]
	  if 'tau_conv' in ana_D_decorr.keys():
	    net.analysis['tau_conv'] = ana_tau_conv['tau_conv'][0,0]

	# break_it == 1 for no topology found, == 2 for already existent
	break_it = net.setup(sim,perturbation,hide=hide)
	
	if break_it == 1:
	  erase(net.Path['script'],net.Path['info_file'],net.Hash)
	elif (break_it == 0) and (mode_sim in ['eft','eft_tmp','mosaic','samples']):
	  net.simulation(sim,hide=hide)
	  
	  #if ('special' in netConst.keys()) and not cluster['q']:
	    #if netConst['special'][0]==3:
	      #analysis_tmp = analyze(mode_sim,net,sim,plot=plot,call=1,save=save,reread=reread)
	      #print analysis_tmp['T_conv']
	      #print analysis_tmp['p']
	      #if not analysis_tmp['p']:
		## get number of steps to jump
		#steps_skip = len(hPconst['cplg']) - step%len(hPconst['cplg']) - 1
		#for i in range(steps_skip):
		  #next(steps_iter)
	perturbation=None
  
  print ' mode_return: %d' % mode_return
  #plt.close('all')
  if mode_return == 1:
    #print analysis.keys()
    return analysis
  
  #print analysis
  if mode == 'r' and (steps > 1 or meta_steps > 1):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['font.size'] = 12
    
    sessionPath = clusterstr + 'data/'
    
    #plot_all = {'eft':plt_eft(analysis,steps),'LE':plt_LE(analysis,steps),'statistics':plt_statistics(analysis,steps)}
    if mode_sim == 'eft':
      plt_eft(analysis,pltDat,str_dict,meta_array,net,save)
      
    elif mode_sim == 'LE':
      plt_LE(analysis,pltDat,str_dict,meta_array,net,save)

    elif mode_sim == 'statistics':
      plt_statistics(analysis,pltDat,str_dict,meta_array,meta_steps_array,net,save,mode_return)
      
    elif mode_sim == 'samples':
      dat_tmp = plt_samples(analysis,pltDat,str_dict,meta_array,meta_steps_array,net,save,mode_return)
      if mode_return == 2:
	return dat_tmp
      
    elif mode_sim == 'drive_cur':
      plt_drive_cur(analysis,pltDat,str_dict,meta_array,net,save,simulation_paras)



def analysis_handover(analysis_tmp,analysis,step,meta_step,steps,meta_steps):
  #print analysis_tmp
  for key in analysis_tmp.keys():
    shrink = 0
    try:
      if not (key in analysis.keys()):		# initiate item if not yet present
	try:
	  analysis[key] = np.zeros(sum(((meta_steps,),(steps,),analysis_tmp[key].shape),()))
	  assert analysis[key].shape[-1]>1
	except:
	  analysis[key] = np.zeros((meta_steps,steps))
	analysis[key][:] = np.nan
      elif not (analysis[key].shape[2:] == analysis_tmp[key].shape) and (not len(analysis_tmp[key])==1 and not (analysis_tmp[key].shape[0] == 1)):	
	new_shape = np.array(analysis[key].shape)
	
	for i in range(len(analysis_tmp[key].shape)):
	  new_shape[-(i+1)] = analysis_tmp[key].shape[-(i+1)]
	  if new_shape[-(i+1)] < analysis[key].shape[-(i+1)]:
	    shrink = 1
	if not shrink:
	  analysis[key].resize(new_shape)	
    except:
      1
    
    if not shrink:
      analysis[key][meta_step][step] = analysis_tmp[key]	# hand over item
    else:
      analysis[key][meta_step][step][:new_shape[2]] = analysis_tmp[key]
      
  return analysis



def initialize_plt(sim_range,steps,scriptName,str_para,topo):
  
  if sim_range == None:
    sim_range = np.array([1])
  
  pltDat = {}
  
  pltDat['topo'] = scriptName
  
  pltDat['str_para'] = str_para
  
  pltDat['steps'] = steps
  
  #if pltDat['topo'] == 'SONET':
  pltDat['min_dat'] = sim_range[0]
  pltDat['max_dat'] = sim_range[-1]
  
  pltDat['range'] = sim_range
  
  #else:
  #pltDat['min_dat'] = np.log10(sim_range.min())
  #pltDat['max_dat'] = sim_range.max())
  #print sim_range
  pltDat['range'] = sim_range
  
  if topo == 'p':
    pltDat['offset'] = (np.log10(pltDat['max_dat'])-np.log10(pltDat['min_dat']))/float(steps)
    pltDat['min_plot'] = np.exp(np.log(pltDat['min_dat'])-pltDat['offset']*2)
    pltDat['max_plot'] = np.exp(np.log(pltDat['max_dat'])+pltDat['offset']/2)
  else:
    pltDat['offset'] = (pltDat['max_dat']-pltDat['min_dat'])/float(steps)
    pltDat['min_plot'] = pltDat['min_dat']-pltDat['offset']
    pltDat['max_plot'] = pltDat['max_dat']+pltDat['offset']
  
  #print pltDat['min_plot']
  #print pltDat['max_plot']
  
  return pltDat



def poisson_set_paras(hPconst,steps):
  
  str_dict = {}
  hP_sim = np.zeros((len(hPconst),steps))

  rate_range = 10**np.linspace(0,4,steps)
  cplg_range = 10**np.linspace(0,4,steps)
  sim_range = rate_range
  
  # fixed poisson variance, varying poisson firing rate
  if hPconst[0] == 'variance':
    hP_sim = np.array([rate_range,-np.sqrt(float(hPconst[1])/rate_range)])
    
    str_dict['title'] = r'$\displaystyle \sigma_{p}^2 = %3.1f$' % hPconst[1]
    str_dict['save'] = 'sigma_p=%3.1f' % hPconst[1]
    str_dict['plot'] = r'$\displaystyle \nu_{p}$'

  # fixed poisson mean input current, varying poisson firing rate
  elif hPconst[0] == 'current':
    hP_sim = np.array([rate_range,-float(hPconst[1])/rate_range])
    
    str_dict['title'] = r'$\displaystyle I_{p} = %3.1f$' % hPconst[1]
    str_dict['save'] = 'I_p=%3.1f' % hPconst[1]
    str_dict['plot'] = r'$\displaystyle \nu_{p}$'
  
  # fixed poisson coupling, varying poisson firing rate 
  elif hPconst[0] == 'cplg':
    assert hPconst[1] < 0, "Please define inhibitory poisson coupling strength (< 0)"
    hP_sim = np.array([rate_range,np.ones(steps)*hPconst[1]])
    
    str_dict['title'] = r'$\displaystyle J_{p} = %3.1f$' % hPconst[1]
    str_dict['save'] = 'J_p=%3.1f' % hPconst[1]
    str_dict['plot'] = r'$\displaystyle \nu_{p}$'
    
  # fixed poisson firing rate, varying poisson coupling
  elif hPconst[0] == 'rate':
    hP_sim = np.array([np.ones(steps)*hPconst[1],-cplg_range])
    sim_range = cplg_range
    str_dict['title'] = r'$\displaystyle \nu_{p} = %3.1f$' % hPconst[1]
    str_dict['save'] = 'nu_p=%3.1f' % hPconst[1]
    str_dict['plot'] = r'$\displaystyle -J_{p}$'
    
  str_dict['rateWnt'] = 'rateWnt=10'
  str_dict['range'] = np.arange(steps)
  
  print str_dict['range']
    
  return hP_sim.transpose(), sim_range, str_dict



def plt_eft(analysis,pltDat,str_dict,meta_array,net,save):
  
  N = meta_array[0]['N']
  y_max = 0.1
  
  if net.topo == 'p':
    steps_cplg = pltDat['steps'][1]
    steps_rate = pltDat['steps'][0]
  #else:
  
  plt.figure(figsize=(3,2.5))
  
  ax = plt.axes([0.21,0.2,0.75,0.75])
  
  if len(meta_array) > 1:
    print pltDat['range']
    pltRange = np.array(pltDat['range'])
    print analysis['eft'][:,0,0]
    if not np.isnan(analysis['eft']).all():
      mask = np.invert(np.isnan(analysis['eft'][:,0,0]))
      
      plt_bootstrap(pltRange,analysis['eft'][:,0],ax,col='k',ls='-',fit_func=func_poly,mask=mask)
  else:
    for meta_step in range(len(meta_array)):
      steps = analysis['p'].shape[1]
      
      N=meta_array[0]['N']
      K=meta_array[0]['K']
      J=meta_array[0]['J0']
      rateWnt=meta_array[0]['rateWnt']
      NeuronType=meta_array[0]['NeuronType']
      
      fileName = 'get_D_decorr.nc'
      D_get_str = './Code/other/calc_I %s %d %d %g' % (fileName,N,K,rateWnt)
      
      runIt(None,D_get_str)
      
      try:
	ncid = netcdf.netcdf_file(fileName,'r')
	
	D_decorr = ncid.variables['D_decorr'].getValue()
	D_spike_fail = ncid.variables['D_spike_failure'].getValue()
	ncid.close()
      except:
	D_decorr = 100
	D_spike_fail = 10
      
      
      else:
	
	print pltDat['range']
	if net.topo == 'p':
	  for j in range(steps_cplg):
	    
	    col = j/float(steps_cplg)
	    c = (col,col,col)
	    
	    J_p = pltDat['range'][j][1]
	    
	    mask_val = np.where(pltDat['range'][:,1]==J_p)[0]
	    mask = np.invert(np.isnan(analysis['eft'][0,mask_val,0]))
	    mask_pos = analysis['eft'][0,mask_val,0] > 0
	    
	    plt_bootstrap(pltDat['range'][mask_val,0],analysis['eft'][0,mask_val],ax,col=c,ls='-',ms=None,mask=mask&mask_pos)
	
	
	#for stp in range(steps):
	  #if not np.isnan(analysis['eft'][0][stp][0]):
	    #if all(analysis['p'][0][stp]==1):
	      #fmt_plt = '-D'
	    #else:
	      #fmt_plt = '-o'
	    
	    #plt.plot(pltDat['range'][stp],analysis['eft'][0][stp],fmt_plt,color='k')
	
	#if net.topo == 'p':
	  #eft = np.zeros((steps_rate,steps_cplg))
	  #eft_var = np.zeros((steps_rate,steps_cplg))
	  #for stp_rate in range(steps_rate):
	    #for stp_cplg in range(steps_cplg):
	      #eft[stp_rate][stp_cplg] = analysis['eft'][0][stp_cplg + stp_rate*steps_cplg,0]
	      #eft_var[stp_rate][stp_cplg] = analysis['eft'][0][stp_cplg + stp_rate*steps_cplg,1]
	      
	  #np.place(eft,eft>D_decorr,D_decorr)
	  #np.place(eft_var,eft>D_decorr,np.nan)
	  #pltDat['range'] = pltDat['range'][::steps_cplg,0]
	  
	  #print analysis['eft']
	  #for stp_cplg in range(steps_cplg):
	    #col = stp_cplg/float(steps_cplg)
	    #c = (col,col,col)
	    #mask = np.invert(np.isnan(eft[:,stp_cplg]))
	    
	    #plt.plot(pltDat['range'][mask],eft[:,stp_cplg][mask],yerr=eft_var[:,stp_cplg][mask],color=c,ecolor='r',fmt='-',linewidth=2)
	    #pltDat['range']
	else:
	  mask = np.invert(np.isnan(analysis['eft'][0,:,0]))
	  
	  plt_bootstrap(pltDat['range'],analysis['eft'][0],ax,col='k',ls='-',mask=mask)
	
      
      #plt.show()
  if net.topo == 'p':
    ax.set_xscale('log')
    ax.set_yscale('log')
  else:
    if len(meta_array) > 1:
      ax.set_xscale('log')
      print pltDat['min_plot'],pltDat['max_plot']
    else:
      plt.ylim([10**(np.floor(np.log10(np.nanmin(analysis['eft'][0])))),10**(np.ceil(np.log10(np.nanmax(analysis['eft'][0]))))])
    
    
    ax.set_yscale('log')
    #plt.xlim([pltDat['min_plot'],pltDat['max_plot']])
    
  #if len(meta_array)>1:
    #plt.xscale('log')
    #eft_plot = analysis['eft'].transpose()[0]
    #plt_range = np.array(pltDat['range'])
    #mask = np.invert(np.isnan(eft_plot))
    
    #m, b, r_value, p_value, std_err = curve_fit(func_exp,plt_range,analysis['eft'][:,0,0])
    #print plt_range**m*np.exp(b)
    #plt.plot(plt_range,plt_range**m*np.exp(b),'--k')
  plt.ylim([10**(-2),10**(1)])
    
  #plt.yscale('log')
  
  
  
  plt.xlabel(r'$\displaystyle \nu_p$ in Hz')
  plt.ylabel(r'$\displaystyle\varepsilon_{FT}$')
  if save:
    save_str = './pics/%s/fluxtube/eft_N=%d.pdf' % (pltDat['topo'],N)
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
  plt.show(block=False)
  
  #mask = np.invert(np.isnan(analysis['pert_size']))
  #epsilon = np.exp(np.linspace(np.log(np.min(analysis['pert_size'][mask])),np.log(np.max(analysis['pert_size'][mask])),1000))
  
  ## construct colorbar
  #steps = analysis['p'].shape[1]
  #ccmap = mcolors.LinearSegmentedColormap.from_list(name='custom',colors=[(0,(0,0,0)),(1,(1-1/steps,1-1/steps,1-1/steps))],N=steps)
  #cbar_dummy = dummy_plot(ccmap,steps,steps+1)
  
  #fmt_plot = ['--o']
  
  #plt.figure()
  #for stp in range(steps):
    
    ##if np.isnan(analysis['eft'][0][stp]):
    #np.place(analysis['p'][0][stp],analysis['pert_size'][0][stp]==0,np.nan)
	   
    #col = stp/float(steps)
    #c = (col,col,col)
    
    #plt.plot(analysis['pert_size'][0][stp],analysis['p'][0][stp],'o',color=c)
    #plt.errorbar(analysis['pert_size'][0][stp],analysis['p'][0][stp],yerr=analysis['p_err'][0][stp],color=c,ecolor=c)
    
    #plt.plot(epsilon,np.exp(-epsilon/analysis['eft'][0][stp,0]),'--',color=c)
  
  #plt.plot([D_decorr,D_decorr],[0.01,1],'k',linewidth=2)
  #plt.plot([D_spike_fail,D_spike_fail],[0.01,1],'r',linewidth=2)
  
  #cbar0 = plt.colorbar(cbar_dummy)
  #cbar0.ax.set_ylabel(str_dict['plot'],fontsize=18)
  #cbar0.set_ticks(np.linspace(1./2,steps-1./2,int(steps/2)+1))
  #cbar0.set_ticklabels(np.log10(pltDat['range'])[::2])
  
  #plt.xlim([10**(-3),10**(1)])
  #if net.topo == 'p':
    #plt.xscale('log')
  
  #plt.ylim([10**(-2),1])
  
  #plt.xlabel(r'$\displaystyle \varepsilon \rightarrow$',fontsize=20)
  #plt.ylabel(r'$\displaystyle f_{RP} \rightarrow$',fontsize=20)
  
  #if save:
    #save_str = './pics/%s/fluxtubes/change_%s_N=%d.pdf' % (pltDat['topo'],str_dict['save'],N)
    #plt.savefig(save_str)
    #print "Figure saved to %s." % save_str
  #plt.show(block=False)
  
  #plt.figure()
  
  
  ##pltDat['max_plot'] = 200
  #plt.xlabel(str_dict['plot'],fontsize=20)
  #plt.ylabel(r'$\displaystyle\varepsilon_{FT}$',fontsize=20)
  ##plt.ylim([0,y_max])
  #plt.yticks(np.linspace(0,y_max,5))
  ##plt.title(r'Impact of $\displaystyle\alpha_{' + str_label[alpha_idx] + '}$')
  #plt.tick_params(labelsize=14)
  #if pltDat['topo'] == 'poisson':
    #pltDat['min_plot'] = max(1,pltDat['min_plot'])
    #plt.xscale('log')
    #plt.yscale('log')
  #plt.xlim([pltDat['min_plot'],pltDat['max_plot']])
  #plt.ylim([10**(-3),10**1])
  #plt.xlabel
  #plt.legend(loc=2)
  
  #plt.show(block=False)
  
  
def plt_LE(analysis,pltDat,str_dict,meta_array,net,save):
  
  
  plt.figure(figsize=(3,2.5))
  
  mp = plt.axes([0.25,0.2,0.6,0.6])
  sp = plt.axes([.65, .6, .3, .3])
  plt_str = r'$\displaystyle %s$' % str_dict['title']
  meta_steps = len(meta_array)
  
  #print analysis.keys()
  for meta_step in range(meta_steps):
    N=meta_array[meta_step]['N']
    K=meta_array[meta_step]['K']
    J=meta_array[meta_step]['J0']
    rateWnt=meta_array[meta_step]['rateWnt']
    NeuronType=meta_array[meta_step]['NeuronType']

    
    #print pltDat['range']
    #print len(meta_array)
    #print analysis['LyapunovExponents'].shape
    #print pltDat['steps']
    if net.topo == 'p':
      for step in range(pltDat['steps'][0]):
	#print step
	#print analysis['LyapunovExponents'][meta_step,step]
	plt_color = (step/float(pltDat['steps'][0]),step/float(pltDat['steps'][0]),step/float(pltDat['steps'][0]))
	
	mp.plot(range(N),analysis['subLyapunovExponents'][meta_step,step],color=plt_color)#,label='%s=%d' % (label_str,hP_range[step]))
	sp.plot(pltDat['range'][step,0],analysis['subLyapunovExponents'][meta_step,step,0],'s',color=plt_color,markersize=3)
	sp.plot(pltDat['range'][step,0],analysis['subLyapunovExponents'][meta_step,step,1],'o',color=plt_color,markersize=3)
    
    else:
      print pltDat['steps']
      for step in range(pltDat['steps']):
	plt_color = (step/float(pltDat['steps']),step/float(pltDat['steps']),step/float(pltDat['steps']))
        mp.plot(range(N),analysis['LyapunovExponents'][meta_step,step,:],color=plt_color)#,label='%s=%d' % (label_str,hP_range[step]))
    #print pltDat['range'][step].shape
    #print analysis['LyapunovExponents'][meta_step,0,1]
    #print pltDat['range'][step,0]
	sp.plot(pltDat['range'][step],analysis['LyapunovExponents'][meta_step,step,1],'.',color=plt_color,markersize=10)
  
  mp.plot([0,N],[-100,-100],'--',color='lightgrey')
  #mp.set_xscale('log')
  mp.set_xlabel('Neuron Index')
  mp.set_xlim([1,N])
  mp.set_xticks([1,10,N])
  mp.set_xticklabels(['1','10','N'],fontsize=12)
  mp.set_ylabel(r'$\displaystyle \lambda_i / s^{-1}$')
  max_y = 0
  min_y = -100
  mp.set_ylim([min_y-10,max_y+30])
  mp.set_yticks(np.linspace(min_y,max_y,3))
  #mp.set_yticklabels(np.linspace(min_y,max_y,3))
  
  # inplot-plot
  #sp.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
  sp.plot(-1,1,'sk',label=r'$\displaystyle \lambda_0$')  
  sp.plot(-1,1,'ok',label=r'$\displaystyle \lambda_1$')
  #sp.set_xticks(pltDat['range'])
  #sp.set_xticklabels(str_dict['range'],fontsize=12)
  sp.set_xlabel(r'$\displaystyle \nu_p$',fontsize=12)
  #sp.set_xlabel(str_dict['plot'],fontsize=12)
  if net.topo == 'p':
    sp.set_xscale('log')
  sp.set_ylim([min_y-10,max_y+10])
  print pltDat['range']
  #sp.set_xticks(pltDat['range'])
  #sp.set_xticks(np.linspace(0,1,3))
  sp.set_yticks(np.linspace(min_y,max_y,3))
  #sp.set_yticklabels(np.linspace(min_y,max_y,3))
  str_title = r'$\displaystyle 1^{st}$ LE'
  #sp.set_title(str_title)
  sp.legend(numpoints=1,bbox_to_anchor=(-0.2,1.2),loc='upper right',borderaxespad=0,prop={'size':10})
  
  mp.spines['top'].set_color('none')
  mp.spines['right'].set_color('none')
  mp.xaxis.set_ticks_position('bottom')
  mp.yaxis.set_ticks_position('left')
  
  #plt.suptitle(r'Lyapunov Spectra for %s' % str_dict['title'],fontsize=18) #, N=%d, K=%d, rate=%3.1f' % (str_dict['title'],N,K,nu),fontsize=20)
  #plt.suptitle(r'Lyapunov Spectra for ER-networks',fontsize=20)
  if save:
    save_str = './pics/%s/LE/spectra_Jp=%d.pdf' % (pltDat['topo'],pltDat['range'][0,1])
    #save_str = './pics/%s/LE/spectra_%s.pdf' % (pltDat['topo'],str(str_dict['save']))
    plt.savefig(save_str)
    print 'Figure saved to ' + save_str
  plt.show()

  
def plt_statistics(analysis,pltDat,str_dict,meta_array,meta_steps_array,net,save,mode_return):
  
  print analysis.keys()
  
  steps = len(pltDat['range'])
  
  N = meta_array[0]['N']
  ## moments, synchrony, input/output correlation
    
  try:
    analysis['phase_moments'] = analysis['phase_moments'].reshape((analysis['phase_moments'].shape[0]/steps,steps,analysis['phase_moments'].shape[1],analysis['phase_moments'].shape[2]))[:,:,0,:]
  except:
    1
    
  N_all = analysis['rateNeurons'].shape[-1]
  if net.topo == 'p':
    analysis['rateNetwork'] = np.nansum(analysis['rateNeurons'][0],axis=1)/N_all	# calculate mean value of each simulation
    analysis['cvNetwork'] = np.nanmean(analysis['cvNeurons'][0],axis=1)	# only take into account spiking neurons for CV
  else:
    analysis['rateNetwork'] = np.nansum(analysis['rateNeurons'][0],axis=1)/N_all	# calculate mean value of each simulation
    analysis['cvNetwork'] = np.nanmean(analysis['cvNeurons'][0],axis=1)	# only take into account spiking neurons for CV
  
  
  #box_plt_rate = [[y for y in row] for row in analysis['rateNeurons'][0]]

  if net.topo == 'p':
    box_plt_CV = [[y for y in row if not np.isnan(y)] for row in analysis['cvNeurons'][0]]
  else:
    box_plt_CV = [[y for y in row if not np.isnan(y)] for row in analysis['cvNeurons'][0]]
  #if net.topo == 'p':
    #if 'const' in str_dict.keys():
      #if str_dict['const'][0]=='rate':
	#contraction(analysis['finalCurrents'][0],-pltDat['range']/np.sqrt(meta_array[0]['K']),plt_range=pltDat['range'])
      #else:
	#contraction(analysis['finalCurrents'][0],np.ones(steps)*str_dict['const'][1]/np.sqrt(meta_array[0]['K']),plt_range=pltDat['range'])
  #else:
    #contraction(analysis['finalCurrents'][0],np.ones(steps)*meta_array[0]['J0']/np.sqrt(meta_array[0]['K']),plt_range=pltDat['range'])
  
  #print pltDat['range']
  if net.topo == 'p':
    pltRange = pltDat['range'][:,0]
  else:
    pltRange = pltDat['range']
  
  plt.figure(figsize=(2,1.5))
  ax = plt.axes([0.25,0.3,0.7,0.65])
  
  print analysis['phase_moments'][0]
  #for i in range(steps):
    #col = i/float(steps)
    #c = (col,col,col)
  ax.errorbar(pltRange,analysis['phase_moments'][0,:,0],yerr=np.sqrt(analysis['phase_moments'][0,:,1]),fmt='-ok',markersize=1)
  ax.plot([0,5000],[0,0],'--',color='grey')
  ax.plot([0,5000],[1,1],'--',color='grey')
  ax.set_xlim([10**(-0.1),10**3.1])
  ax.set_ylim([-1.5,1.2])
  ax.set_yticks(np.linspace(-1,1,3))
  ax.set_xscale('log')
  #ax.set_yscale('log')
  ax.set_xlabel(r'$\displaystyle \nu_p$ in Hz')
  ax.set_ylabel(r'$\displaystyle \langle \phi \rangle$')
  print pltDat['range']
  if save:
    save_str = './pics/%s/statistics/mean_phases_Jp=%d.pdf' % (pltDat['topo'],pltDat['range'][0,0])
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
  plt.show(block=False)
  
  print analysis['D_spike_fail'][0,:,0]
  print pltRange
  
  plt.figure(figsize=(3,1.5))
  ax = plt.axes([0.25,0.3,0.7,0.65])
  plt_bootstrap(pltRange,analysis['D_spike_fail'][0],ax,col='k')
  #ax.plot(pltRange,analysis['D_spike_fail'][0],'k')
  ax.set_xscale('log')
  ax.set_xlim([10**(-0.1),10**3.1])
  ax.set_xlabel(r'$\displaystyle \nu_p$ in Hz')
  ax.set_yscale('log')
  ax.set_ylim([10**-1,10**1.1])
  ax.set_ylabel(r'$\displaystyle D_{\phi}^{sf}$')
  if save:
    save_str = './pics/%s/statistics/D_spike_fail_Jp=%d.pdf' % (pltDat['topo'],pltDat['range'][0,1])
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
  plt.show(block=False)
  
  #print 'currents: ', analysis['finalCurrents'][0]
  #tauM = 0.01
  
  #T_free = -tauM*np.log(1.-1./(1+analysis['finalCurrents'][0]))
  #print 'T free: ', T_free
  
  #phi_mean = 0.5
  #print 'mean phase: ', analysis['phase_moments'][0,:,0]
  #PRC = - tauM/T_free * np.log(np.exp(-analysis['phase_moments'][0,:,0]*T_free/tauM) + net.Const['J0']/(net.Const['K']*(1+analysis['finalCurrents'][0]))) - analysis['phase_moments'][0,:,0]
  
  #print 'PRC: ', PRC
  #print 'average spike failure distance: ', PRC*net.Const['K']
  
  
  #print analysis['D_spike_fail']
  
  
  
  plt.figure(figsize=(6,5))
  
  ax0 = plt.axes([0.13,0.6,0.37,0.25])
  ax1 = plt.axes([0.13,0.1,0.37,0.45],sharex=ax0)
  ax2 = plt.axes([0.53,0.8,0.37,0.18],sharex=ax0)
  ax21 = plt.axes([0.53,0.6,0.37,0.18],sharex=ax0)
  ax3 = plt.axes([0.53,0.1,0.37,0.45],sharex=ax0)
  
  plt.setp(ax0, xticks=pltRange, xticklabels=pltRange)
  plt.setp(ax1.get_xticklabels(), visible=False)
  plt.setp(ax2.get_xticklabels(), visible=False)
  plt.setp(ax3.get_xticklabels(), visible=False)
  
  ax_labels = plt.axes([0.125,0.075,0.8,0.7],frameon=False)
  #ax_labels.set_xlabel(str_dict['plot'],fontsize=14)
  ax_labels.set_xlabel(r'$\displaystyle \nu_p$ in Hz',fontsize=14)
  #ax_labels.get_yaxis().set_visible(False)
  plt.setp(ax_labels, xticks=[], xticklabels=[])
  plt.setp(ax_labels, yticks=[], yticklabels=[])
      
  if net.topo == 'p':
    wid = pltRange/1.5
  else:
    wid = (pltDat['range'][-1]- pltDat['range'][0])/10.
  
  
  if net.topo == 'p' or not len(meta_steps_array):
    ax0.plot(pltRange,analysis['corr_KNu'][0],'--ok',label=r'$\displaystyle cor(K^{in},\nu)$')
    ax0.plot(pltRange,analysis['corr_NuCV'][0],'--oy',label=r'$\displaystyle cor(CV,\nu)$')
  else:
    #if meta_steps > 1:
    ax0.plot(pltRange,analysis['corr_KNu'],'--ok',label=r'$\displaystyle cor(K^{in},\nu)$')
    ax0.plot(pltRange,analysis['corr_NuCV'],'--oy',label=r'$\displaystyle cor(CV,\nu)$')
  
  ax0.legend(bbox_to_anchor=(0.5,1.05),loc='lower center',borderaxespad=0,prop={'size':12})
  #ax0.legend(loc=2,prop={'size':12})
  
  
  boxprops = dict(linestyle='-', linewidth=1, color='k')
  flierprops = dict(marker='+', markerfacecolor='k',color='k', markersize=2,linestyle='none')
  whiskerprops = dict(linestyle='--', linewidth=1, color='k')
  
  if net.topo == 'p':
    ax1.boxplot(analysis['rateNeurons'][0].transpose(),notch=1,positions=pltRange,widths=wid,boxprops=boxprops,flierprops=flierprops,whiskerprops=whiskerprops)
    ax1.plot(pltRange,analysis['rateNeurons_peaks'][0,:,0],'.',color='grey',markersize=10)
    rate_max = 5
  else:
    ax1.boxplot(analysis['rateNeurons'][0].transpose(),notch=1,positions=pltRange,widths=wid,boxprops=boxprops,flierprops=flierprops,whiskerprops=whiskerprops)
    ax1.plot(pltRange,analysis['rateNeurons_peaks'][0,:,0],'.',color='grey',markersize=10)
    rate_max = 4
  
  
  #if net.topo == 'p':
  ax2.plot(pltRange,analysis['finalCurrents'][0],'--or')
  #else:
    #ax2.plot(pltRange,analysis['finalCurrents'][0],'--or')
  ax2.yaxis.set_label_position('right')
  ax2.yaxis.set_ticks_position('right')
  
  #print analysis['chi']
  ax21.plot(pltRange,analysis['chi'][0],'--ok')
  ax21.set_ylim([0,1])
  ax21.yaxis.set_label_position('right')
  ax21.yaxis.set_ticks_position('right')
  ax21.set_yticks(np.linspace(0,1,3))
  
  
  boxprops = dict(linestyle='-', linewidth=1, color='k')
  flierprops = dict(marker='+', markerfacecolor='k',color='k', markersize=2,linestyle='none')
  whiskerprops = dict(linestyle='--', linewidth=1, color='k')
  
  #print analysis['cvNetwork']
  ax3.boxplot(box_plt_CV,notch=1,positions=pltRange,widths=wid,boxprops=boxprops,flierprops=flierprops,whiskerprops=whiskerprops)
  ax3.plot(pltRange,analysis['cvNetwork'],'.',color='k',markersize=10,label='mean')
  ax3.plot(-10,-10,'o',color='grey',label='peaks')
  ax3.plot(pltRange,analysis['cvNeurons_peaks'][0,:,0],'.',color='grey',markersize=10)
  ax3.yaxis.set_label_position('right')
  ax3.yaxis.set_ticks_position('right')
  ax3.legend(loc=1,prop={'size':12})
  
  if net.topo == 'p':
    plt.setp(ax0, xscale='log', xlim=[10**(-0.5),10**3.5], yticks=np.linspace(-1,1,5), ylim=[-1,1])
    plt.setp(ax1, xscale='log', yticks=np.linspace(0,rate_max,6),ylim=[0,rate_max])
    plt.setp(ax2, xscale='log', yticks=np.linspace(0,rate_max,6), yscale='log')
    plt.setp(ax21, xscale='log')
    plt.setp(ax3, xscale='log', yticks=np.linspace(0,2,5), ylim=[0,3])
  else:
    plt.setp(ax0, xlim=[pltDat['min_plot'],pltDat['max_plot']], yticks=np.linspace(-1,0,3), ylim=[-1,0.25])
    plt.setp(ax1, yticks=np.linspace(0,rate_max,5),ylim=[0,rate_max])
    currents_max = max(analysis['finalCurrents'][0])*1.2
    plt.setp(ax2, yscale='log')#,yticks=np.linspace(0,currents_max,6), ylim=[0,currents_max])
    plt.setp(ax3, yticks=np.linspace(0,2,5), ylim=[0,2])
  
  ax0.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '%g'%y))
  
  plt.setp(ax0, ylabel='Pearson correlation')
  plt.setp(ax1, ylabel=r'$\displaystyle \nu$ in Hz')
  plt.setp(ax2, ylabel=r'$\displaystyle I_{ext}$')
  plt.setp(ax21, ylabel=r'$\displaystyle \chi$')
  plt.setp(ax3, ylabel=r'$\displaystyle CV_{ISI}$')
  #[pltDat['min_plot'],pltDat['max_plot']]
  #plt.suptitle(str_dict['title'],fontsize=20)
  if save:
    save_str = './pics/%s/statistics/moments_%s_N=%d.pdf' % (pltDat['topo'],str_dict['save'],N)
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
  plt.show(block=False)

  if 'Xcor' in analysis.keys():
    print analysis['Xcor']
    #print analysis['high_Xcor']
    #print analysis['max_Xcor']
    plt.figure()
    gs = gridspec.GridSpec(4,2)
  
    ax0 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax1 = plt.subplot(gs[1:4,0:2])
    #ax2 = plt.subplot(gs[2,0]) 
    
    col = [(x,x,x) for x in np.linspace(0,0.8,steps)]
    for i in range(steps):
      norm = np.sum(analysis['yhist_Xcor'][0][i])
      
      ax0.plot(pltDat['range'][i],np.mean(analysis['Xcor'][0][i]),'o',color=col[i])
      ax1.plot(analysis['xhist_Xcor'][0][0],analysis['yhist_Xcor'][0][i]/float(norm),color=col[i])
      ax2.plot(pltDat['range'][i],analysis['chi'][0][i],'o',color=col[i],linewidth=4)

    #ax1.plot(pltDat['range'],np.mean(analysis['max_Xcor'][0],axis=1))
    ax0.plot([-10,10],[0,0],'--k',linewidth=0.2)
    ax0.xaxis.tick_top()
    ax0.xaxis.set_label_position('top')
    ax0.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
    ax0.set_yticks(np.linspace(-0.01,0.01,3))
    ax0.set_ylim([-0.01,0.01])
    ax0.set_xlabel(str_dict['plot'],fontsize=18)
    ax0.set_ylabel(r'$\displaystyle [cor(\{t\}_i,\{t\}_j)]_{i,j}$',fontsize=18)
    
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
    ax2.set_xlabel(str_dict['plot'],fontsize=18)
    ax2.set_ylabel(r'Synchrony $\displaystyle\chi$',fontsize=18)
    ax2.set_ylim([0,0.4])
    ax2.set_yticks(np.linspace(0,0.4,5))
    ax2.set_yticklabels(np.linspace(0,0.4,5),fontsize=14)
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_ticks_position('right')
    
    ax1.set_yscale('log')
    ax1.set_xlim([-0.7,0.7])
    ax1.set_xlabel(r'$\displaystyle cor(\{t\}_i,\{t\}_j)$',fontsize=18)
    ax1.set_ylim([10**(-6),1])
    ax1.set_ylabel('p(cor)',fontsize=18)
    #ax1.set_ylim([-1,1])
    #ax2.plot(pltDat['range'],np.mean(analysis['high_Xcor'][0],axis=1))
    if save:
      save_str = './pics/%s/statistics/spike_corr_%s.pdf' % (pltDat['topo'],str_dict['save'])
      plt.savefig(save_str)
      print "Figure saved to %s." % save_str
    plt.show(block=False)
    
  
  if not net.topo == 'p':
    print 'plotting degrees'
    ## in- and out-degree
    col = [(x,x,x) for x in np.linspace(0,0.8,steps)]
    
    ccmap = mcolors.LinearSegmentedColormap.from_list(name='custom',colors=[(0,(0,0,0)),(1,(0.8,0.8,0.8))],N=steps)
    cbar_dummy = dummy_plot(ccmap,steps,steps+1)
    
    plt.figure(figsize=(4,3))
    
    ax0 = plt.axes([0.15,0.75,0.66,0.15])
    ax1 = plt.axes([0.15,0.15,0.32,0.575])
    ax2 = plt.axes([0.49,0.15,0.32,0.575])
    ax3 = plt.axes([0.85,0.15,0.015,0.75])
    #gs = gridspec.GridSpec(4, 21)
    
    #ax0 = plt.subplot(gs[0,0:20])
    #ax1 = plt.subplot(gs[1:4,0:10])
    #ax2 = plt.subplot(gs[1:4,10:20])
    #ax3 = plt.subplot(gs[0:4,20])
    
    
    y_border = 0.1
    
    K = 0.1*N
    ax_border = [2*K,2*K]
    
    for i in range(steps):
      print col[i]
      ax0.plot(pltRange[i],analysis['corr_KK'][0,i],'.',markersize=8,color=col[i])
    ax1.hist(analysis['inDegree'][0].transpose(),bins=100,range=(0,ax_border[0]),linewidth=1,color=col,alpha=0.8,histtype='step')#,edgecolor="none")
    ax2.hist(analysis['outDegree'][0].transpose(),bins=100,range=(0,ax_border[1]),linewidth=1,color=col,alpha=0.8,histtype='step')#,edgecolor="none")
      
    ax1.set_xlabel(r'$\displaystyle K^{in}$',fontsize=14)
    ax1.set_xlim([0,ax_border[0]])
    ax1.set_xticks(np.linspace(0,ax_border[0],5))
    #ax1.set_xticklabels(np.linspace(0,ax_border[0],5).astype('int'))
    ax1.set_ylabel(r'$\displaystyle \rho(K)$',fontsize=14)
    ax1.set_ylim([0,y_border*N_all])
    ax1.set_yticks([])
    #ax1.set_yticks(np.linspace(0,y_border*N_all,3))
    #ax1.set_yticklabels(np.linspace(0,0.1,3))
    
    ax2.set_xlabel(r'$\displaystyle K^{out}$',fontsize=14)
    ax2.set_xlim([0,ax_border[1]])      
    ax2.set_xticks(np.linspace(0,ax_border[1],5))
    #ax2.set_xticklabels(np.linspace(0,ax_border[1],5).astype('int'))
    ax2.set_ylim([0,y_border*N_all])
    #ax2.set_yticks(np.linspace(0,y_border*N_all,3))
    ax2.set_yticklabels([])

    ax0.plot([pltDat['min_plot'],pltDat['max_plot']],np.zeros(2),'--',color='grey')
    ax0.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
    ax0.set_xticks(pltDat['range'])
    #ax0.set_xticklabels(pltDat['range'])
    ax0.set_ylim([-1,1])
    ax0.set_yticks(np.linspace(-1,1,3))
    #ax0.set_yticklabels(np.linspace(-1,1,3))
    ax0.set_ylabel(r'$\displaystyle cor(K^{in},K^{out})$')
    ax0.xaxis.set_label_position('top')
    ax0.xaxis.set_ticks_position('top')
    
    cbar0 = plt.colorbar(cbar_dummy,cax=ax3)
    cbar0.ax.set_ylabel(str_dict['plot'],fontsize=14)

    cbar0.set_ticks(np.linspace(1./2,steps-1./2,steps))
    cbar0.set_ticklabels(pltDat['range'])
    cbar0.ax.tick_params(labelsize=10)
    
    if save:
      save_str = './pics/%s/statistics/in_out_degree_%s_N=%d.pdf' % (pltDat['topo'],str_dict['save'],N)
      plt.savefig(save_str)
      print "Figure saved to %s." % save_str
    plt.show(block=False)
  
  
  ### spikecrossing statistics
  if not 'mean_ISIdecorr' in analysis.keys():
    return
  for key in ['mean_ISIdecorr','cv_ISIdecorr']:
    analysis[key] = analysis[key].transpose((1,2,0))	# rearrange, such that (data,simulation)  
    mask = np.isnan(analysis[key])
    analysis[key] = np.ma.array(analysis[key],mask=mask)
    analysis['net'+key] = np.mean(analysis[key],axis=0)	# only take into account spiking neurons for CV
  
  box_plt_meanISI0 = [[y for y in row if y] for row in analysis['mean_ISIdecorr'][0].T]
  box_plt_meanISI1 = [[y for y in row if y] for row in analysis['mean_ISIdecorr'][1].T]
  box_plt_cvISI0 = [[y for y in row if y] for row in analysis['cv_ISIdecorr'][0].T]
  box_plt_cvISI1 = [[y for y in row if y] for row in analysis['cv_ISIdecorr'][1].T]
  
  plt.figure()
  gs = gridspec.GridSpec(2, 1)

  ax0 = plt.subplot(gs[0,0])
  ax1 = plt.subplot(gs[1,0])
  x_max = 20
  
  ax0.plot(pltDat['range'],10*np.log(1+1/(np.sqrt(20)*analysis['finalCurrents'])),'Dr',markersize=10,label=r'$\displaystyle T_{free}$')
  
  
  wid=0.25
  ylim0 = 20
  ylim1 = 5
  
  bp00 = ax0.boxplot(box_plt_meanISI0,notch=1,positions=pltDat['range']-pltDat['offset']/5.,widths=wid,sym='')
  bp01 = ax0.boxplot(box_plt_meanISI1,notch=1,positions=pltDat['range']+pltDat['offset']/5.,widths=wid,sym='')
  
  ax0.plot(-1,-1,color='black',label='pre-post spike')
  ax0.plot(-1,-1,color='red',label='post-pre spike')
  
  bp10 = ax1.boxplot(box_plt_cvISI0,notch=1,positions=pltDat['range']-pltDat['offset']/5.,widths=wid,sym='')
  bp11 = ax1.boxplot(box_plt_meanISI1,notch=1,positions=pltDat['range']+pltDat['offset']/5.,widths=wid,sym='')
  
  plt.setp([bp00['boxes'],bp10['boxes']], color='black')
  #plt.setp([bp00['whiskers'],bp10['whiskers']], color='black')
  plt.setp([bp00['fliers'],bp10['fliers']], color='black')

  plt.setp([bp01['boxes'],bp11['boxes']], color='red')
  #plt.setp([bp01['whiskers'],bp11['whiskers']], color='red')
  plt.setp([bp01['fliers'],bp11['fliers']], color='red')

  ax0.set_ylim([0,ylim0])
  ax1.set_ylim([0.5,ylim1])
  
  ax0.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
  ax1.set_xlim([pltDat['min_plot'],pltDat['max_plot']])
  
  ax1.set_xlabel(str_dict['plot'],fontsize=18)
  
  ax0.legend(loc=2,prop={'size':12})
  ax0.set_ylabel('average ISI decorr',fontsize=18)
  ax1.set_ylabel('CV of ISI decorr',fontsize=18)
  plt.setp([ax0,ax1], xticks=pltDat['range'], xticklabels=pltDat['range'])
  plt.setp(ax0, yticks=np.linspace(0,ylim0,6), yticklabels=np.linspace(0,ylim0,6))
  plt.setp(ax1, yticks=np.linspace(0.5,ylim1,3), yticklabels=np.linspace(0.5,ylim1,3))
  
  #if pltDat['topo'] == 'poisson':
    #ax0.set_xscale('log')
    #ax1.set_xscale('log')
    
  if save:
    save_str = './pics/%s/statistics/ISIdecorr_%s.pdf' % (pltDat['topo'],str_dict['save'])
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
    
  plt.show(block=False)
    
    
def plt_samples(analysis,pltDat,str_dict,meta_array,meta_steps_array,net,save,mode_return):
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  ###plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
  
  mpl.rcParams['xtick.labelsize'] = 12
  mpl.rcParams['ytick.labelsize'] = 12
  mpl.rcParams['font.size'] = 12
  
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Arial')
  
  #mpl.rcParams['xtick.labelsize'] = 25
  #mpl.rcParams['ytick.labelsize'] = 25
  #mpl.rcParams['font.size'] = 30
    
  if net.topo == 'p':
    steps_cplg = pltDat['steps'][1]
    steps_rate = pltDat['steps'][0]
    
    steps = steps_cplg
  else:
    steps = 1
    #pltDat['steps']
  
  
  
  if net.Const['special'] == 0:
  #if 'convergeVelocity' in analysis.keys():
    print 'plot convergence statistics'
    print analysis.keys()
    print analysis['1stLE']
    
    plt.figure(figsize=(6,4))
      
    ax0 = plt.axes([0.22,0.2,0.73,0.52])
    
    if net.topo == 'p':
      
      for j in range(steps-1,-1,-1):
	J_p = pltDat['range'][j][1]
	
	mask_val = np.where(pltDat['range'][:,1]==J_p)[0]
	
	D_decorr = analysis['tau_conv_first'][0,mask_val]
	pltRange = pltDat['range'][mask_val,0]
      
	col = 1-(j+1)/float(steps)
	ck = (col,col,col)
	cr = ((j+1)/float(steps),0,0)
	
	
	mask = np.invert(np.isnan(analysis['tau_conv_first'][0,mask_val,0]))
	mask_sigma = np.abs(analysis['tau_conv_later'][0,mask_val,0]) > 0.02
	plt_bootstrap(np.abs(pltRange),np.abs(analysis['tau_conv_first'][0,mask_val]),ax0,ck,'--',mask=mask)#&mask_sigma)
	
	mask = np.invert(np.isnan(analysis['tau_conv_later'][0,mask_val,0]))
	plt_bootstrap(np.abs(pltRange),np.abs(analysis['tau_conv_later'][0,mask_val]),ax0,ck,'-',label=r'$\displaystyle J_p=%g$'%J_p,mask=mask)#,fit_func=func_poly,p0=[-1,10])#&mask_sigma)
	
	print analysis['1stLE'][0,mask_val,1]
	print analysis['1stLE'][0,mask_val,0]
	ax0.plot(np.abs(pltRange),1./np.abs(analysis['1stLE'][0,mask_val,0]),'s',color=ck,markersize=4)
	ax0.plot(np.abs(pltRange),1./np.abs(analysis['1stLE'][0,mask_val,1]),'o',color=ck,markersize=4)
	
      plt.xlabel(r'$\displaystyle \nu_p$ in Hz')
      plt.plot(0,0,'--k',label=r'$\displaystyle \tau_{C}^{(1)}$')
      plt.plot(0,0,'-k',label=r'$\displaystyle \tau_{C}^{(2)}$')
      
      plt.plot(0,0,'s',color=ck,markersize=4,label=r'$\displaystyle \lambda_0$')
      plt.plot(0,0,'o',color=ck,markersize=4,label=r'$\displaystyle \lambda_1$')
      
      #plt.legend(loc=1,prop={'size':10},ncol=2)
      
    else:
      #if len(meta_steps_array):
	##pltRange = np.array([0.2,0.5,1,2,5,10,20])
	#meta_steps = meta_steps_array[0]
	pltRange = abs(np.array(pltDat['range']))
	#print pltRange
	#print analysis['tau_conv_later'][:,0]
	#for meta_step in range(meta_steps):
	  #col = 1-(meta_step+1)/float(meta_steps_array[0])
	  #ck = (col,col,col)
	  ##print meta_array[meta_step]
	  #mask = np.invert(np.isnan(analysis['tau_conv_later'][meta_step::meta_steps_array[0],0,0]))
	  #plt_bootstrap(pltRange,np.abs(analysis['tau_conv_later'][meta_step::meta_steps_array[0],0]),ax0,ck,'-',r'$\displaystyle K=%d$'%meta_array[meta_step]['K'],mask=mask)
	  #ax0.set_xlim([0.1,50])
	#plt.xlabel(r'$\displaystyle \bar{\nu}$ in Hz',fontsize=12)
      #else:
      
	mask = np.invert(np.isnan(analysis['tau_conv_later'][:,0,0]))
	plt_bootstrap(pltRange,np.abs(analysis['tau_conv_later'][:,0]),ax0,'k','-',mask=mask)
	  
	  #plt.errorbar(np.abs(pltRange),np.abs(analysis['tau_conv_2nd'][i::meta_steps_array[1],0,0]),yerr=analysis['tau_conv_var'][i::meta_steps_array[1],0,0],fmt='o-',color='k',ecolor='k',label=r'initial $\displaystyle \tau_{conv}$')
	  #plt.errorbar(np.abs(pltRange),np.abs(analysis['tau_conv'][i::meta_steps_array[1],0,1]),yerr=analysis['tau_conv_var'][i::meta_steps_array[1],0,1],fmt='o-',color='r',ecolor='r',label=r'latter $\displaystyle \tau_{conv}$')
	plt.xlim([pltRange[0]*0.5,pltRange[-1]*1.5])
	#plt.xlabel(r'$\displaystyle \bar{\nu}$ in Hz',fontsize=12)
	plt.xlabel(r'$N$',fontsize=12)
	#plt.plot(np.abs(pltDat['range']),np.abs(pltDat['range'])**slope1*np.exp(intercept1),'--k',label=r'%s$\displaystyle^{%5.3f}$'%(str_dict['plot'],slope1))
    #plt.plot(np.abs(pltDat['range']),-conv1,'ko')
    plt.xlim([10**-0.1,10**3.1])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([10**(-2.5),10**(1)])
      
      #label_plot = 'later convergence'
    #else:
      #label_plot = r'$\displaystyle v_{conv}$'
    #plt.show()
    ###plt.plot(np.abs(pltDat['range']),np.abs(pltDat['range'])**slope2*np.exp(intercept2),'--r',label=r'%s$\displaystyle^{%5.3f}$'%(str_dict['plot'],slope2))
    
    plt.plot([0,100000],[0.01]*2,'--r',linewidth=2)
    plt.annotate(r'$\displaystyle -\tau_m$',xy=[0.12,0.01],xytext=[0.2,0.0035],arrowprops=dict(arrowstyle="->"),fontsize=10)
    #plt.errorbar(np.abs(pltDat['range']),-mean2,yerr=var2,color='r',ecolor='r')###,label=label_plot)
    
    plt.legend(ncol=3,bbox_to_anchor=(-0.2,1.07),loc='lower left',borderaxespad=0,prop={'size':10})
    ###plt.xlim([10**1,10**4])
    #plt.ylim([10**(-4),10**(0)])

    
    
    plt.ylabel(r'$\displaystyle \tau_{C}/s$',fontsize=12)
    #plt.xticks(pltDat['range'][::2])
    #plt.xticklabels(pltDat['range'],fontsize=14)
    
    
    #plt.setp(plt, xticks=, xticklabels=pltDat['range'])
      #mask = np.invert(np.isnan(analysis['convergeVelocity1'][0]))
      #print analysis['convergeVelocity1'].shape
      #print mask.shape
      #print np.mean(analysis['convergeVelocity1'][0],axis=1)
    
    #plt.xscale('log')
    #plt.yscale('log')
    
    #plt.title(str_dict['title'])
    if save:
      save_str = './pics/%s/samples/convergenceSpeed_N=%d_nu=%g.pdf' % (pltDat['topo'],meta_array[0]['N'],meta_array[0]['rateWnt'])
      plt.savefig(save_str)
      print 'figure saved as: %s ' % save_str
    plt.show(block=False)
    
    
  if net.Const['special'] == 2:
    
    if net.topo == 'p':
      slope_J = np.zeros(pltDat['steps'][1])
      slope_J_err = np.zeros(pltDat['steps'][1])
      
      J_p_array = np.zeros(pltDat['steps'][1])

    else:
      D_decorr = np.zeros((len(meta_array),steps))
      
      if len(meta_array) > 1:
	#print analysis['D_decorr']
	D_decorr = analysis['D_decorr_bs'][:,:,0]
	D_decorr_var = analysis['D_decorr_bs'][:,:,1]
      else:
	D_decorr = analysis['D_decorr_bs'][0,:,0]
	D_decorr_var = analysis['D_decorr_bs'][0,:,1]
      D_decorr = D_decorr.reshape((meta_steps_array[0],meta_steps_array[1]))
      D_decorr_var = D_decorr_var.reshape((meta_steps_array[0],meta_steps_array[1]))

    
    #print analysis['D_decorr_bs']
    plt.figure(figsize=(5,4))
    
    pltRange = -pltDat['range'][:steps_cplg,1]
    ax0 = plt.axes([0.23,0.23,0.7,0.7])
    
    for j in range(steps_rate):
      col = j/float(steps_rate)
      #col = 1-(j+1)/float(steps_rate)
      c = (col,col,col)
      print j, steps_rate
      if net.topo == 'p':
	
	nu_p = pltDat['range'][j*steps_cplg,0]
	#nu_p_array[j] = pltDat['range'][j*steps_cplg,1]
	
	#mask_val = np.where(pltDat['range'][:,1]==J_p)[0]
	D_decorr = analysis['D_decorr_bs'][0,j*steps_cplg:(j+1)*steps_cplg]
	
	
	#mask = np.invert(np.isnan(D_decorr[:,0])) & (D_decorr[:,0] >= 1)
	mask = np.invert(np.isnan(D_decorr[:,0]))
	mask_sigma = D_decorr[:,2]-D_decorr[:,1] <= 0.5
	
	#plt_bootstrap(pltRange,D_decorr,ax0,c,label=r'$\displaystyle J_p=%g$'%J_p,mask=mask&mask_sigma)
	#popt,perr = plt_bootstrap(pltRange,D_decorr,ax0,c,'-',label=r'$\displaystyle J_p=%g$'%J_p,fit_func=func_exp,p0=[-1,10],mask=mask&mask_sigma)
	print pltRange
	print D_decorr[mask]
	#popt,perr = 
	plt_bootstrap(pltRange,D_decorr,ax0,c,'-',ms=None,mask=mask&mask_sigma,label=r'$\displaystyle \nu_p = %g \,K \bar{\nu}$'%(pltDat['range'][j*steps_cplg,0]/float(net.Const['K'])))
	
	#idx = 3
	#if np.isfinite(D_decorr[idx,0]):
	  #ax0.text(pltRange[idx], D_decorr[idx,0], r'$\displaystyle \nu_p = %g K \bar{\nu}$'%(pltDat['range'][j*steps_cplg,0]/float(net.Const['K'])),fontsize=20)
	
	#ax1.errorbar(-J_p,-popt[0],yerr=perr[0],fmt='.',color=c,markersize=10,ecolor='r',label='$\displaystyle J_p = %g$' % J_p)
	
	#slope_J[j] = popt[0]
	#slope_J_err[j] = perr[0]
	
	#ax0.set_xlim([10**(0),10**(3.1)])
	
    #ax0.set_yticks(10**np.linspace(0,1,6))
    #ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
    ax0.set_yscale('log')
    ax0.set_ylim([10**0,10**1.5])
    ax0.set_xlabel(r'$\displaystyle J_p$')
    ax0.set_ylabel(r'$\displaystyle D_{\phi}^{dc}$')
    
    ax0.set_xticks(np.linspace(0,20,5))
    ax0.set_xticklabels(np.linspace(0,-20,5).astype('int'))
    
    ax0.spines['right'].set_color('none')
    ax0.yaxis.set_ticks_position('left')
    ax0.spines['top'].set_color('none')
    ax0.xaxis.set_ticks_position('bottom')
    
    ax0.legend(prop={'size':25},frameon=False,loc=[0.37,0.5],numpoints=1,handlelength=1)
    
    if save:
      save_str = './pics/%s/samples/LIF_D_decorr.svg' % (pltDat['topo'])
      print 'figure saved as: %s ' % save_str
      plt.savefig(save_str)
      
    plt.show()
	
    #for j in range(steps_rate):
      #col = 1-(j+1)/float(steps_rate)
      #c = (col,col,col)
      
      #if net.topo == 'p':
	  
	#ax0 = plt.axes([0.1,0.23,0.5,0.7])
	#ax1 = plt.axes([0.625,0.58,0.25,0.3])
	
	#J_p = pltDat['range'][j][1]
	#J_p_array[j] = pltDat['range'][j][1]
	
	#mask_val = np.where(pltDat['range'][:,1]==J_p)[0]
	#D_decorr = analysis['D_decorr_bs'][0,mask_val]
	#pltRange = pltDat['range'][mask_val,0]
	
	#mask = np.invert(np.isnan(D_decorr[:,0])) & (D_decorr[:,0] >= 1)
	#mask_sigma = D_decorr[:,2]-D_decorr[:,1] <= 0.5
	
	##print D_decorr
	##print D_decorr[mask&mask_sigma]
	
	##plt_bootstrap(pltRange,D_decorr,ax0,c,label=r'$\displaystyle J_p=%g$'%J_p,mask=mask&mask_sigma)
	##popt,perr = plt_bootstrap(pltRange,D_decorr,ax0,c,'-',label=r'$\displaystyle J_p=%g$'%J_p,fit_func=func_exp,p0=[-1,10],mask=mask&mask_sigma)

	#popt,perr = plt_bootstrap(pltRange,D_decorr,ax0,c,'-',ms='.',label=r'$\displaystyle J_p=%g$'%J_p,mask=mask&mask_sigma,fit_func=func_exp_exp,p0=[-1,10])
	
	#ax1.errorbar(-J_p,-popt[0],yerr=perr[0],fmt='.',color=c,markersize=10,ecolor='r',label='$\displaystyle J_p = %g$' % J_p)
	
	#slope_J[j] = popt[0]
	#slope_J_err[j] = perr[0]
	
	#ax0.set_xlim([10**(0),10**(3.1)])
	
	#ax0.set_yticks(10**np.linspace(0,1,6))
	#ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
	#ax0.set_ylim([10**0,10**1.1])
    
      #else:
	#ax0 = plt.axes([0.2,0.2,0.7,0.7])
	#pltRange = pltDat['range']
	#print pltRange
	##pltRange = [500,1000,2000,5000,10000]
	#pltRange = [0.2,0.5,1,2,5,10,20]
	##print D
	#print meta_steps_array
	#for i in range(meta_steps_array[1]):
	  #col = i/float(meta_steps_array[1])
	  #c = (col,col,col)
	  
	  #ax0.errorbar(pltRange,D_decorr[:,i],D_decorr_var[:,i],color=c,ecolor='r',fmt='-o',label=r'$\displaystyle J_0 = %d$' %meta_array[i]['J0'])
	  ##ax0.set_xlabel(r'$\displaystyle \bar{\nu}$ in Hz',fontsize=12)
	  ##ax0.set_xlabel('N',fontsize=12)
	  ##ax0.set_xlim([0.1,50])
	  ##ax1.plot(0,np.exp(intercept),'s',color=c)
	  ##ax1.plot(1,slope,'o',color=c)
	  ##ax1.set_xticks([0,1])
	  ##ax1.set_xticklabels(['offset','slope'],fontsize=10)
	  ##ax1.set_yticks(np.linspace(0,1,3))
	  ##ax1.set_yticklabels(np.linspace(0,1,3))
	  #ax0.set_xlim([10**(-1),10**(1.5)])
	  #ax0.set_ylim([10**0.5,10**2])

	  ##ax1.set_xlim([-0.5,1.5])
	  ##ax1.set_ylim([0,1])
	  
	##ax0.errorbar(pltDat['range'],D_decorr[1],D_decorr_var[1],color=c,ecolor='r',fmt='-o')
      
      
      
      ##plt.errorbar(pltDat['range'],analysis['D_decorr'][0],analysis['D_decorr_var'][0],ecolor='k',fmt='-ok')
        #norm = pltRange[0]*10.
    
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    
    ax0.set_xlabel(r'$\displaystyle \nu_p$ in Hz',fontsize=12)
    
    
    #if net.topo == 'p':
      #if len(J_p_array) > 1:
	#popt,pcov = curve_fit(func_log,-J_p_array,-slope_J,sigma=slope_J_err)
	
	#perr = np.sqrt(np.diag(pcov))
	#print "2. popt: ", popt
	#print "2. perr: ", perr
	
	#ax1.xaxis.set_label_position('top')
	#ax1.yaxis.set_ticks_position('right')
	#ax1.yaxis.set_label_position('right')
	#ax1.plot(-J_p_array,popt[0]*np.log(-J_p_array)+popt[1],'--r',linewidth=0.5)
      
	#plt.setp(ax1,xscale='log',xlim=[-J_p_array[0],-J_p_array[-1]],xticks=-J_p_array,xlabel=r'$\displaystyle -J_p$')
	#ax1.set_xticklabels(-J_p_array,fontsize=10)
	#ax1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: '%g' % x))
	
	#ax1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: '%g' % x))
	#plt.setp(ax1,ylim=[0,0.01],yticks=np.linspace(0,0.01,3),ylabel=r'$\displaystyle -\gamma/$ms')
	##ax1.set_yticklabels(np.linspace(0,0.01,3),fontsize=10)
	
	#ax1.text(0.5, -0.015, r'$\displaystyle \gamma = m\cdot \log{(J_p)} + b$ \newline $\displaystyle m=(%5.3g\pm%5.3g)10^3$ \newline $\displaystyle b=(%5.3g\pm%5.3g)10^3$ '%(popt[0]*10**3,perr[0]*10**3,popt[1]*10**3,perr[1]*10**3), bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=12)
      
    
    
    ax0.set_ylabel(r'$\displaystyle D_{\phi}^{dc}$')
    #ax0.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0,prop={'size':10})
    
    if save:
      save_str = './pics/%s/samples/D_decorr_N=%d_K=%d_nu=%d.pdf' % (pltDat['topo'],meta_array[0]['N'],meta_array[0]['K'],meta_array[0]['rateWnt'])
      print 'figure saved as: %s ' % save_str
      plt.savefig(save_str)
    plt.show(block=False)
    
    if mode_return == 2:
      return popt, perr
    
    
  if net.Const['special'] in [3,30]:
    print 'plot re-convergence statistics'
    
    #print analysis['tau_RC']
    #print pltDat['range']
    
    #pltRange = pltDat['range'][:steps,1]
    ##print "he"
    ##print pltRange
    
    #plt.figure(figsize=(5,4))
    
    
    #ax = plt.axes([0.23,0.23,0.7,0.7])
    
    #for i in range(steps_rate):
      #col = i/float(steps_rate)
      #c = (col,col,col)
      #mask= np.invert(np.isnan(analysis['tau_RC'][0,i*steps_cplg:(i+1)*steps_cplg]))
      #ax.plot(-pltRange[mask],analysis['tau_RC'][0,i*steps_cplg:(i+1)*steps_cplg][mask],'-',color=c,label=r'$\displaystyle \nu_p = %g K \bar{\nu}$'%(pltDat['range'][i*steps_cplg,0]/100.),linewidth=2)
      #idx = 3
      
      #if np.isfinite(analysis['tau_RC'][0,i*steps_cplg:(i+1)*steps_cplg][mask][idx]):
	#print analysis['tau_RC'][0,i*steps_cplg:(i+1)*steps_cplg][mask][4]
	#ax.text(-pltRange[mask][idx], 1.2*analysis['tau_RC'][0,i*steps_cplg:(i+1)*steps_cplg][mask][idx], r'$\displaystyle \nu_p = %g K \bar{\nu}$'%(pltDat['range'][i*steps_cplg,0]/100.),fontsize=20)
    
    #ax.set_yscale('log')
    #ax.set_xlabel(r'$\displaystyle J_p$')
    #ax.set_ylabel(r'$\displaystyle \tau_{RC} / s$')
    #ax.set_ylim([10**(-2),10**2])
    
    #ax.set_xticks(np.linspace(0,20,5))
    #ax.set_xticklabels(np.linspace(0,-20,5).astype('int'))
    
    #ax.spines['right'].set_color('none')
    #ax.yaxis.set_ticks_position('left')
    #ax.spines['top'].set_color('none')
    #ax.xaxis.set_ticks_position('bottom')
    #if save:
      #save_str = './pics/%s/LIF_tau_RC.svg' % (pltDat['topo'])
      #print 'figure saved as: %s ' % save_str
      #plt.savefig(save_str)
      
    ##ax.legend(prop={'size':10})
    #plt.show()
    
    
    plt.figure(figsize=(2,2))
    pltRange = pltDat['range'][::steps,0]
    ##print pltRange
    ax = plt.axes([0.34,0.25,0.63,0.7])
    ##ax = plt.axes([0.32,0.25,0.65,0.7])
    ##ax1 = plt.axes([0.175,0.125,0.8,0.4])
    
    for i in range(steps):
      col = i/float(steps)
      c = (col,col,col)
      ax.plot(pltRange,analysis['tau_RC'][0][i::steps],'-',color=c)

      #ax1.plot(pltRange,analysis['T_min'][0][i::steps],color=c)
    
    
    ax.set_xscale('log')
    
    ax.set_yscale('log')
    ax.set_ylim([10**(-2),10**(2)])
    ax.set_xlim([10**1.9,10**3.1])
    ax.set_xlabel(r'$\displaystyle \nu_p$ in Hz',fontsize=14)
    #plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_ylabel(r'$\displaystyle \tau_{RC}$ in s',fontsize=14)
    
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    #ax1.set_ylim([10**(-3),10**(0)])
    #ax1.set_xlim([10**1.9,10**3.1])
        
    #ax1.set_xlabel(r'$\displaystyle \nu_p/Hz$',fontsize=14)
    #ax1.set_ylabel(r'$\displaystyle \tau^{offset}/s$',fontsize=14)
    
    
    if save:
      save_str = './pics/%s/J_nu_phase_N=%d.pdf' % (pltDat['topo'],meta_array[0]['N'])
      print 'figure saved as: %s ' % save_str
      plt.savefig(save_str)

    plt.show(block=False)
    
  #steps = len(pltDat['range'])
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='serif')
  if 'drop_sensitivity' in analysis.keys():	# sensitivity horizon!
    print 'plot drop sensitivity (whats that again?)'
    print analysis['drop_sensitivity']
    col = [(x,x,x) for x in np.linspace(0,0.8,steps)]
    plt.figure()
    for i in range(len(analysis['drop_sensitivity'][0])):
      plt.plot(analysis['measureTimes'][0][0],analysis['div_ratio'][0][i],color=col[i],linewidth=2,label=r'$\displaystyle \alpha_{conv}=%g$'%pltDat['range'][i])
      #plt.plot(analysis['measureTimes'][0][0],np.exp(analysis['drop_sensitivity'][0][i][0]*analysis['measureTimes'][0][0])*np.exp(analysis['drop_sensitivity'][0][i][1]),'--r')
    plt.xlim([0,1])
    plt.ylim([0.01,1])
    plt.ylabel('fraction of non-diverged trajectories',fontsize=14)
    plt.xlabel('time in s',fontsize=16)
    #plt.yscale('log')
    plt.legend(prop={'size':16})
        
    plt.show(block=False)
  
  if net.Const['special'] == 4:
    print 'plotting phase correlations:'
    print analysis.keys()
    
    J_p = pltDat['range'][0,1]
    mask_val = np.where(pltDat['range'][:,1]==J_p)[0]
    pltRange = pltDat['range'][mask_val,0]
    
    col = [(x,x,x) for x in np.linspace(0,0.8,steps)]
    
    fig = plt.figure(figsize=(6,2))
    ax_corr = plt.axes([0.02,0.22,0.3,0.73])
    ax0 = plt.axes([0.35,0.22,0.25,0.73])
    ax1 = plt.axes([0.65,0.22,0.25,0.33])
    ax2 = plt.axes([0.65,0.6,0.25,0.35],sharex=ax1)
    
    
    ax1.plot(pltRange,analysis['assigned'][0,:,1]/analysis['assigned'][0,:,0],'-k',linewidth=0.5)
    ax1.set_xscale('log')
    ax1.set_xticks([])
    ax1.yaxis.set_ticks_position('right')
    ax1.yaxis.set_label_position('right')
    ax1.set_xlabel(r'$\displaystyle \nu_p$ in Hz')
    ax1.set_xlim([10**(-0.1),10**(2.6)])
    ax1.set_ylabel(r'$\displaystyle f_{assigned}$')
    ax1.set_yticks(np.linspace(0.1,0.9,5))
    ax1.set_ylim([0.1,0.6])
    
    
    ax2.plot(pltRange,analysis['mean_dt'][0],'-k',linewidth=0.5)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim([10**(-0.1),10**(2.6)])
    ax2.set_ylabel(r'$\displaystyle \langle |\Delta t| \rangle$/s')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    ax2.set_yticks(np.logspace(-3,-1,5))
    ax2.set_ylim([10**(-3),10**(-1)])
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    phi_corr_network_mean = np.zeros((steps,2))
    
    for i in range(steps):
      phi_corr_network_mean[i,0] = np.mean(analysis['phi_corr'][0,i])
      phi_corr_network_mean[i,1] = np.sqrt(np.var(analysis['phi_corr'][0,i]))
      
      print 'mean: ', phi_corr_network_mean[i]
      
      ax_corr.hist(analysis['phi_corr'][0,i].ravel(),bins=501,range=[-1,1],histtype='step',color=col[i],label=r'$\displaystyle \nu_p=%d$'%int(pltDat['range'][i][0]))
      
      ax0.plot(analysis['deltaT_hist_range'][0,i,:-1],analysis['deltaT_hist'][0,i],color=col[i])#,label=r'$\displaystyle \nu_p=%d$'%int(pltDat['range'][i][0]))
      ax1.plot(pltRange[i],analysis['assigned'][0,i,1]/analysis['assigned'][0,i,0],'o',markersize=4,color=col[i])
      ax2.plot(pltRange[i],analysis['mean_dt'][0,i],'o',markersize=4,color=col[i])
    
    ax_corr.set_xlabel(r'corr($\displaystyle \phi_n^i,\phi_n^j$)')
    ax_corr.set_xticks(np.linspace(0,0.8,5))
    ax_corr.set_xlim([-0.15,0.85])
    ax_corr.set_yticks([])
    ax_corr.spines['top'].set_color('none')
    ax_corr.spines['left'].set_color('none')
    ax_corr.spines['right'].set_color('none')
    #ax_corr.legend(bbox_to_anchor=(0.8,0.95),loc='upper left',borderaxespad=0,prop={'size':10},ncol=2)
    ax_corr.xaxis.set_ticks_position('bottom')
    ax_corr.yaxis.set_ticks_position('none')
    
    ax0.set_xscale('log')
    ax0.set_xlabel(r'$\displaystyle |\Delta| t$ in s')
    ax0.set_yticks([])
    ax0.spines['top'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('none')
    #ax0.spines['top'].set_color('none')
    
    #plt.legend(loc='upper center',prop={'size':10},ncol=3)
    
    if save:
      save_str = './pics/%s/assigned_spikes.pdf' % (pltDat['topo'])
      print 'figure saved as: %s ' % save_str
      plt.savefig(save_str)
    plt.show(block=False)
    
    
    
    
    #phi_corr_network_mean = np.zeros((steps,2))
    
    #plt.figure(figsize=(2,2))
    #ax0 = plt.axes([0.05,0.2,0.9,0.75])
    ##ax0 = plt.axes([0.05,0.525,0.9,0.45])
    ##ax1 = plt.axes([0.2,0.15,0.75,0.2])
    #for i in range(steps):
      #print i,steps
      
      #phi_corr_network_mean[i,0] = np.mean(analysis['phi_corr'][0,i])
      #phi_corr_network_mean[i,1] = np.sqrt(np.var(analysis['phi_corr'][0,i]))
      
      #print 'mean: ', phi_corr_network_mean[i]
      
      #ax0.hist(analysis['phi_corr'][0,i].ravel(),bins=501,range=[-1,1],histtype='step',color=col[i],label=r'$\displaystyle \nu_p=%d$'%int(pltDat['range'][i][0]))
      ##ax1.errorbar(pltRange[i],phi_corr_network_mean[i,0],yerr=phi_corr_network_mean[i,1],fmt='-o',color=col[i])
      
      ##print 'len: ', len(phi_corr_mean)
    
    ##ax0.set_xlabel(r'corr($\displaystyle \phi_n^i,\phi_n^j$)')
    ##ax0.set_xticks(np.linspace(0,0.8,5))
    ##ax0.set_xlim([-0.15,0.85])
    ##ax0.set_yticks([])
    
    ##ax1.errorbar(pltRange,phi_corr_network_mean[:,0],yerr=phi_corr_network_mean[:,1],fmt='-o',color='k')
    ##ax1.set_xscale('log')
    ##ax1.set_xlim([10**(-0.1),10**(2.1)])
    ##ax1.set_ylim([-0.1,0.6])
    ##ax1.set_yticks([0,0.5])
    ##ax1.set_xlabel(r'$\displaystyle \nu_p$ in Hz')
    ##ax1.set_ylabel(r'$\displaystyle \langle corr(\phi_n^i,\phi_n^j) \rangle$')
    ##ax1.yaxis.set_ticks_position('right')
    ##ax0.legend(bbox_to_anchor=(0.1,0.95),loc='upper left',borderaxespad=0,prop={'size':10},ncol=3)
    #if save:
      #save_str = './pics/%s/phase_corr.pdf' % (pltDat['topo'])
      #print 'figure saved as: %s ' % save_str
      #plt.savefig(save_str)
    #plt.show(block=False)
    


  
def plt_drive_cur(analysis,pltDat,str_dict,meta_array,save,simulation_paras):
  
  J = simulation_paras.transpose()[1]
 
  plt.figure()
  plt.plot(pltDat['range'],analysis['I_bal'][0],'-k',label='balance equation')
  plt.plot(pltDat['range'],analysis['I_selfcon'][0],'or',label='selfconsistency')
  plt.plot(pltDat['range'],analysis['I_sim'][0],'ob',label='simulations')
  plt.xscale('log')
  plt.yscale('log')
  plt.ylim([10**(-1.5),10**2])
  plt.xlabel(str_dict['plot'])
  plt.xticks(pltDat['range'])
  plt.legend(loc=2)
  plt.show(block=False)
  
  #print analysis
  #print pltDat
  #print simulation_paras
  contraction(I_ext=analysis['I_sim'][0],J=J)


def scriptSONet(alpha_steps,alpha_intervalls,steps=None,mode_sim=['eft'],mean_idx=None,TC=0.2,mode=None,N=1000,save=0):
  execfile('py_code/tubesize.py')
	      
    #elif (('pert_spread' in mode_sim) and ('pert_spread_para' in dat)):
      #analysis['time_div'][sim_idx], analysis['vel_div'][sim_idx], analysis['time_mean'][sim_idx], analysis['vel_mean'][sim_idx] = analyze(mode_sim,fileName_in=sessionPath + 'info/' + dat,scriptName_in='SONET',call=1)
      #sim_idx += 1
	      
    
    ### spikecrossing statistics
    #t_post = analysis['t_post'].transpose(0,2,1)
    #t_pre = analysis['t_pre'].transpose(0,2,1)

    #gs = gridspec.GridSpec(2, 1)

    #ax0 = plt.subplot(gs[0,0])
    #ax1 = plt.subplot(gs[1,0])

    #for i in range(alpha_steps[alpha_idx]):
      #mask = np.invert(np.isnan(t_post[i][0]))
      #ax0.hist(t_post[i][0][mask],bins=10,range=[0,0.025])
      #mask = np.invert(np.isnan(t_pre[i][0]))
      #ax1.hist(-t_pre[i][0][mask],bins=10,range=[0,0.025])
    
    #plt.show(block=False)
    
    

  #if 'pert_spread' in mode_sim:
    
    #mask = np.array(analysis['time_div']==0)
    
    #analysis['time_div'] = np.ma.array(analysis['time_div'],mask=mask)
    
    #data_box_time = [[y for y in row if y] for row in analysis['time_div']]
    #data_box_vel = [[y for y in row if y] for row in analysis['vel_div']]
    
    ##print data_box_time
    
    #fig,axes = plt.subplots(nrows=1,ncols=2)
    ##print analysis['time_div'][mask].transpose().shape
    #axes[0].boxplot(data_box_time,notch=1,positions=alpha_range[alpha_idx])
    #axes[0].plot(alpha_range[alpha_idx],analysis['time_mean'],'ro',markersize=5,label='mean')
    ##axes[0].plot(alpha_range[alpha_idx],analysis['peaks'][0][0],'o',color='grey',markersize=10)
    #axes[0].set_ylabel('Time to divergence',fontsize=18)
    #axes[0].set_ylim([0,0.12])
    #axes[0].set_xlabel(r'$\displaystyle\alpha_{' + str_label[alpha_idx] + '}$',fontsize=18)
    #axes[0].legend()

    #axes[1].boxplot(data_box_vel,notch=1,positions=alpha_range[alpha_idx])
    #axes[1].plot(alpha_range[alpha_idx],analysis['vel_mean'],'ro',markersize=5,label='mean') #axes[0].plot(alpha_range[alpha_idx],analysis['peaks'][0][0],'o',color='grey',markersize=10)
    #axes[1].yaxis.set_label_position('right')
    #axes[1].yaxis.set_ticks_position('right')
    #axes[1].set_ylabel('Velocity of divergence',fontsize=18)
    #axes[1].set_ylim([0,3000])
    #axes[1].set_xlabel(r'$\displaystyle\alpha_{' + str_label[alpha_idx] + '}$',fontsize=18)
    #axes[1].legend()
    
    #if save:
      #save_str = './pics/SONET/divergence_' + str_label[alpha_idx] + '.pdf'
      #plt.savefig(save_str)
      
    #plt.show()
    
    
    #fig,axes = plt.subplots(nrows=2,ncols=1)
    
    #for i in range(valid_idx):
      #mask = np.array(analysis['time_div'][i]!=0)
      #axes[0].hist(analysis['time_div'][i][mask],bins=41,range=[0,0.1],alpha=0.5,label='%d'%i)
      #axes[1].hist(analysis['vel_div'][i][mask],bins=41,range=[0,3000],alpha=0.5,label='%d'%i)
    #axes[0].legend()
    #axes[1].legend()
    #plt.show()
	
	
def dummy_plot(ccmap,border,steps):
  plt.figure()
  # Using contourf to provide my colorbar info, then clearing the figure
  Z = [[0,0],[0,0]]
  levels = np.linspace(0,border,steps)
  cbar_dummy = plt.contourf(Z, levels, cmap=ccmap)
  plt.clf()
  return cbar_dummy
  
  

def plot_distributions(moments,str_label,str_value,plot,save):
  # find peaks of the moments of the distributions
  peak_search_range = 4

  idx_mom = 0
  peaks = np.zeros((2,2,4))	#store up to 4 peaks (moments,axes,peaks)
  peaks[:] = np.nan

  num_bin = 21
  moment = np.zeros((2,2,num_bin))	#store moments (moments,axes,bins)
  
  moment_bars = [50,1.5]
  
  for i in range(len(moment)):
    moment[i][0] = np.histogram(moments[i],bins=num_bin,range=[0,moment_bars[i]])[0]/1000.
    moment[i][1] = np.histogram(moments[i],bins=num_bin,range=[0,moment_bars[i]])[1][:-1]
    
    peak_idx = 0
    
    for idx_dat in range(len(moment[i][0])):
      peak = 1
      
      if idx_dat < peak_search_range:		# define range for search of local maximum
	moment_search = moment[i][0][:idx_dat+peak_search_range+1]	
      elif idx_dat > (num_bin-peak_search_range+1):
	moment_search = moment[i][0][idx_dat-peak_search_range:]
      else:
	moment_search = moment[i][0][idx_dat-peak_search_range:idx_dat+peak_search_range+1]
      
      if np.sum(moment[i][0][idx_dat] <= moment_search)==1:
	peaks[i][0][peak_idx] = moment[i][1][idx_dat]
	peaks[i][1][peak_idx] = moment[i][0][idx_dat]
	peak_idx += 1
      
    idx_mom += 1
    
  if plot:
    # plot the moment distributions for each simulation
    
    plt.figure()
    gs = gridspec.GridSpec(4, 1)
      
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
      
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[1:4,0])
    mask = np.invert(np.isnan(moments[0]))
    ax0.boxplot(moments[0][mask],notch=1,widths=0.75,vert=False)
    ax0.plot(peaks[0][0],np.ones(len(peaks[0][0])),'o',color='grey',markersize=10)
    ax0.set_xlim([0,50])
    ax0.set_xticks(np.linspace(0,50,6))
    ax0.set_xticklabels([])
    ax0.set_yticklabels([])
    
    ax1.plot(moment[0][1],moment[0][0],'.-')
    ax1.plot(peaks[0][0],peaks[0][1],'o',color='grey',markersize=10)
    ax1.set_ylim([0,0.2])
    ax1.set_xticks(np.linspace(0,50,6))
    ax1.set_xticklabels(np.linspace(0,50,6),fontsize=16)
    ax1.set_yticks(np.linspace(0,0.2,5))
    ax1.set_yticklabels(np.linspace(0,0.2,5),fontsize=16)
    ax1.set_xlabel(r'$\displaystyle \bar{\nu} \rightarrow$',fontsize=24)
    ax1.set_ylabel(r'$\displaystyle p(\bar{\nu}) \rightarrow$',fontsize=24)
    plt.suptitle(r'Firing Rates of $\displaystyle\alpha_{' + str_label + '} = ' + str_value +'$',fontsize=24)
    #plt.suptitle(r'Firing Rates of an Erd\H{o}s-R\'{e}nyi-network',fontsize=24)
    
    if save:
      save_str = './pics/SONET_statistics/rate_' + str_label + '=' + str_value + '.pdf'
      plt.savefig(save_str)
    plt.show(block=False)
    
    
    plt.figure()
    gs = gridspec.GridSpec(4, 1)

    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[1:4,0])
    
    mask = np.invert(np.isnan(moments[1]))
    ax0.boxplot(moments[1][mask],notch=1,widths=0.75,vert=False)
    ax0.plot(peaks[1][0],np.ones(len(peaks[0][1])),'o',color='grey',markersize=10)
    ax0.set_xlim([0,1.4])
    ax0.set_xticks(np.linspace(0,1.4,8))
    ax0.set_xticklabels([])
    ax0.set_yticklabels([])
    
    ax1.plot(moment[1][1],moment[1][0],'.-')
    ax1.plot(peaks[1][0],peaks[1][1],'o',color='grey',markersize=10)
    ax1.set_xlim([0,1.4])
    ax1.set_ylim([0,0.2])
    ax1.set_xticks(np.linspace(0,1.4,8))
    ax1.set_xticklabels(np.linspace(0,1.4,8),fontsize=16)
    ax1.set_yticks(np.linspace(0,0.2,5))
    ax1.set_yticklabels(np.linspace(0,0.2,5),fontsize=16)
    ax1.set_xlabel(r'CV $\displaystyle \rightarrow$',fontsize=24)
    plt.suptitle(r'CV of $\displaystyle\alpha_{' + str_label + '} = ' + str_value +'$',fontsize=24)
    #plt.suptitle(r'CV of an Erd\H{o}s-R\'{e}nyi-network',fontsize=24)
    if save:
      save_str = './pics/SONET_statistics/CV_' + str_label + '=' + str_value + '.pdf'
      plt.savefig(save_str)
    plt.show()
    
    ### 3rd & 4th moment - not so useful
    
    #axes[1,0].plot(moments[2][0],moments[2][1],'.-')
    #axes[1,0].plot(moments[2][0][peaks[2]],moments[2][1][peaks[2]],'r.',markersize=10)
    #axes[1,0].set_xlim([-1,4])
    #axes[1,0].set_ylim([0,2])
    #axes[1,0].set_title('skewness')
    
    #axes[1,1].plot(moments[3][0],moments[3][1],'.-')
    #axes[1,1].plot(moments[3][0][peaks[3]],moments[3][1][peaks[3]],'r.',markersize=10)
    #axes[1,1].set_xlim([0,20])
    #axes[1,1].set_ylim([0,0.5])
    #axes[1,1].set_title('kurtosis')
    
    #save_str = './pics/SONET_moments/' + str_label + '=' + str_value + '.pdf'
    
    #plt.suptitle(r'$\displaystyle\alpha_{' + str_label + '} = ' + str_value +'$',fontsize=18)
    #plt.savefig(save_str)
  
    ### end 3rd & 4th moment
  return peaks
  
  
  
######################################## POISSON SCRIPT ################################################### 
def scriptPoisson(hP_steps,hP_intervall=None,mode_sim=['eft'],hPconst=None,constVal=None,scaling='pow',TC=0.2,mode=None,N=1000,save=0):
  
  execfile('py_code/tubesize.py')
  
  sim_min, sim_max = hP_intervall[0], hP_intervall[1]
  if scaling == 'pow':
    hP_range = 10**np.linspace(sim_min,sim_max,hP_steps)
  elif scaling == 'lin':
    hP_range = np.linspace(sim_min,sim_max,hP_steps)

  if not mode:
    mode = raw_input('(r)ead or (s)imulate data? ')
  
  hPplot = []
  sessionPath = 'data/'
  
  analysis = {}
  if 'eft' in mode_sim:
    analysis = {'eft':np.zeros(hP_steps), 'var_eft':np.zeros(hP_steps)}
  
  if 'statistics' in mode_sim:
    analysis = {'moments':np.zeros((hP_steps,2,N)),'chi':np.zeros(hP_steps),'peaks':np.zeros((np.prod(alpha_steps),2,2,4)),'in_degree':np.zeros((np.prod(alpha_steps),N)),'in_degreeNet':np.zeros(np.prod(alpha_steps)),'out_degree':np.zeros((np.prod(alpha_steps),N)),'out_degreeNet':np.zeros(np.prod(alpha_steps)),'corr_KNu':np.zeros(np.prod(alpha_steps)),'corr_KK':np.zeros(np.prod(alpha_steps)),'t_post':np.zeros((np.prod(alpha_steps),N,3)),'t_pre':np.zeros((np.prod(alpha_steps),N,3))}
    
    analysis = {'moments':np.zeros((hP_steps,2,N)),'chi_tmp':np.zeros(hP_steps),'peaks':np.zeros((hP_steps,2,2,4))}
    
  if 'LE' in mode_sim:
    if TC < 1.:
      TC = 1.
    analysis = {'LE':np.zeros((hP_steps,N))}
  
  idx = 0
  
  for hP_step in hP_range:
    
    hPrate, hPcplg, hPplot, string = set_paras(hP_step,hPconst,constVal,hPplot)
    
    if mode == 's':
      
      net = network(mode='t',topo='p',hPrate=hPrate,hPcplg=hPcplg,N=N)
      sim = simulation_data(net,mode_sim=mode_sim,TC=TC)
      print '\nNow simulating network with hPrate=%5.3f, hPcplg=%5.3f.' % (hPrate,hPcplg)#,plus_string)
      break_it = net.setup(sim,hPconst=hPconst,constVal=constVal,check=0)
  
      if break_it == 1:
	erase(net.Path['script'],net.Path['info_file'],net.Hash)
      elif ((break_it == 0) and ('eft' in mode_sim)):
	net.simulation(sim)
	
    elif mode == 'r':
      
      str_hPrate = ('hPrate=%5.3f' % hPrate).rstrip('0').rstrip('.') + '_'
      str_hPcplg = ('hPcplg=%5.3f' % hPcplg).rstrip('0').rstrip('.') + '_'
      
      str_hPconst = {'hP_variance':r'\sigma_p^2','hP_current':r'I_p','hP_rate':r'\bar{\nu}_p','hP_cplg':r'J_p'}
      
      for dat in os.listdir(sessionPath + 'info/'):
	if 'poisson' in dat:

	  if ((str_hPrate in dat) and (str_hPcplg in dat) and (string['rateWnt'] in dat) and (('_N=%d_'%N) in dat)):
	    fileName=sessionPath + 'info/' + dat
		    
	    if (('statistics' in mode_sim) and ('statistics' in dat)):
	      print 'For %s, %s, reading statistics.' % (str_hPrate.split('_')[0],str_hPcplg.split('_')[0])
	      
	      Data = analyze(mode_sim,fileName_in=fileName,scriptName_in='poisson',call=1)
	      
	      analysis['moments'][idx][0] = Data['rateNeurons']
	      analysis['moments'][idx][1] = Data['cvNeurons']

	      analysis['chi_tmp'][idx] = Data['chi']
	      
	      # plot moments and calculate peaks
	      analysis['peaks'][idx] = plot_distributions(analysis['moments'][idx],'label string','irgendein string',plot=0,save=save)
	      
	    
	    elif ('eft' in mode_sim) and ('poisson_para' in dat):
	      print 'For %s, %s, reading flux tube size.' % (str_hPrate.split('_')[0],str_hPcplg.split('_')[0])
	      #try:
	      analysis['eft'][idx],analysis['var_eft'][idx],numsim = analyze(mode_sim,fileName_in=fileName,scriptName_in='poisson',call=1)
	      #except:
		#print 'No data for %s, %s.' % (str_hPrate,str_hPcplg)
		#analysis['eft'] = None
		
	    
	    elif ('LE' in mode_sim) and ('LE' in dat):
	      print 'For %s, %s, reading LE spectrum.' % (str_hPrate.split('_')[0],str_hPcplg.split('_')[0])
	      Data = analyze(mode_sim,fileName_in=fileName,scriptName_in='poisson',call=1)
	      analysis['LE'][idx] = Data['subLyapunovExponents']
	      
      idx += 1
  
  if mode == 'r':
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    if 'eft' in mode_sim:
      for i in range(hP_steps):
	if not analysis['eft'][i]:
	  analysis['eft'][i] = None
	  analysis['var_eft'][i] = None
      
      plt.rc('text', usetex=True)
      plt.rc('font', family='serif')
      plt.figure()
	
      plt.plot(hPplot,analysis['eft'],'s',markersize=10)
      plt.errorbar(hPplot,analysis['eft'],yerr=analysis['var_eft'],ecolor='r',elinewidth=2,fmt=None)
      plt.xlabel(r'$\displaystyle' + string['plot'] + '$',fontsize=18)
      plt.ylabel(r'$\displaystyle\varepsilon_{FT}$',fontsize=18)
      plt.xlim([min(hPplot),max(hPplot)])
      if scaling == 'pow':
	plt.xscale('log')
      plt.title(r'$\displaystyle' + string['title'] + '$',fontsize=20)
      plt.savefig('./pics/poisson_eft/' + string['save'] + '.pdf')
      plt.show(block=False)
  
  
    if 'statistics' in mode_sim:
      
      data_box = np.zeros((4,N,hP_steps))
      
      mean_plot = np.zeros((4,hP_steps))
      
      analysis['peaks'] = analysis['peaks'].transpose(1,2,0,3)	#rearrange, such that (moment,axes,simulation,peaks)

      data_box = analysis['moments'].transpose((1,2,0))	# rearrange, such that (moment,data,simulation)
      
      mask = np.isnan(data_box)
      
      data_box = np.ma.array(data_box,mask=mask)
      mean_plot = np.mean(data_box,axis=1)	# calculate mean value of each simulation
      
      data_box_rate = [[y for y in row if y] for row in data_box[0].T]
      data_box_CV = [[y for y in row if y] for row in data_box[1].T]
      
      #print mean_plot
      
      chi = analysis['chi_tmp'][:idx]
      
      gs = gridspec.GridSpec(3, 2)
      
      ax0 = plt.subplot(gs[0,0])
      ax1 = plt.subplot(gs[1:3,0])
      ax2 = plt.subplot(gs[0:3,1])
      
      #offset = (alpha_intervalls[alpha_idx][1]-alpha_intervalls[alpha_idx][0])/float(alpha_steps[alpha_idx])
      
      ax0.plot(np.log10(hP_range),chi,'r--',linewidth=4)
      #ax0.set_xlim([alpha_intervalls[alpha_idx][0]-offset,alpha_intervalls[alpha_idx][1]+offset])
      ax0.set_xticks(np.log10(hP_range))
      ax0.set_xticklabels(hP_range,fontsize=14)
      
      ax0.set_ylabel(r'Synchrony $\displaystyle\chi$',fontsize=18)
      ax0.set_ylim([0,0.4])
      ax0.set_yticks(np.linspace(0,0.4,5))
      ax0.set_yticklabels(np.linspace(0,0.4,5),fontsize=14)
      
      ax1.boxplot(data_box_rate,notch=1,positions=np.log10(hP_range))
      ax1.plot(np.log10(hP_range),analysis['peaks'][0][0],'o',color='grey',markersize=10)
      #ax1.set_xlabel(r'$\displaystyle\alpha_{' + str_label[alpha_idx] + '}$',fontsize=18)
      #ax1.set_xlim([alpha_intervalls[alpha_idx][0]-offset,alpha_intervalls[alpha_idx][1]+offset])
      ax1.set_xticks(np.log10(hP_range))
      ax1.set_xticklabels(hP_range,fontsize=14)
      
      ax1.set_ylabel(r'Firing rate $\displaystyle \nu$',fontsize=18)
      ax1.set_ylim([0,100])
      ax1.set_yticks(np.linspace(0,100,6))
      ax1.set_yticklabels(np.linspace(0,100,6),fontsize=14)
      
      
      ax2.boxplot(data_box_CV,notch=1,positions=np.log10(hP_range))
      ax2.plot(np.log10(hP_range),mean_plot[1],'s',color='k',markersize=10,label='mean')
      ax2.plot(-10,-10,'o',color='grey',label='peaks')
      ax2.plot(np.log10(hP_range),analysis['peaks'][1][0],'o',color='grey',markersize=10)
      #ax2.set_xlabel(r'$\displaystyle\alpha_{' + str_label[alpha_idx] + '}$',fontsize=18)
      #ax2.set_xlim([alpha_intervalls[alpha_idx][0]-offset,alpha_intervalls[alpha_idx][1]+offset])
      ax2.set_xticks(np.log10(hP_range))
      ax2.set_xticklabels(hP_range,fontsize=14)
      
      ax2.set_ylabel('Coefficient of variation',fontsize=18)
      ax2.set_ylim([0,2])
      ax2.set_yticks(np.linspace(0,2,5))
      ax2.set_yticklabels(np.linspace(0,2,5),fontsize=14)
      ax2.yaxis.set_label_position('right')
      ax2.yaxis.set_ticks_position('right')
      ax2.legend()
      
      plt.suptitle(r'Constant $\displaystyle %s=%d$' % (str_hPconst[hPconst],constVal),fontsize=20)
      if save:
	save_str = './pics/poisson_statistics/const_%s=%d.pdf' % (hPconst,constVal)
	plt.savefig(save_str)
      plt.show(block=False)
      
      
    if 'LE' in mode_sim:
      plt.figure()
      mp = plt.axes()
      sp = plt.axes([.55, .5, .3, .3])
      if hPconst == 'hP_rate':
	label_str = r'$\displaystyle J_p$'
      else:
	label_str = r'$\displaystyle \bar{\nu}_p$'

      for i in range(hP_steps):
	plt_color = (i/float(hP_steps),i/float(hP_steps),i/float(hP_steps))
	mp.plot(range(N),analysis['LE'][i],color=plt_color)#,label='%s=%d' % (label_str,hP_range[i]))
	sp.plot(np.log10(hP_range)[i],analysis['LE'].transpose()[1][i],'o',color=plt_color,markersize=10)
      mp.set_xlabel('Neuron Index i/N',fontsize=18)
      mp.set_xticks(np.linspace(0,N,6))
      mp.set_xticklabels(np.linspace(0,1,6),fontsize=14)
      mp.set_ylabel('Lyapunov Exponents',fontsize=18)
      mp.set_yticks(np.linspace(-100,0,6))
      mp.set_yticklabels(np.linspace(-100,0,6),fontsize=14)
      
      # try of inplot-plot
      #n, bins, patches = hist(s, 400, normed=1)
      sp.set_xlim([np.log10(hP_range)[0]-0.5,np.log10(hP_range)[-1]+0.5])
      sp.set_xticks(np.log10(hP_range))
      sp.set_xticklabels(np.log10(hP_range).astype('int'),fontsize=14)
      sp.set_yticks(np.linspace(-100,0,3))
      sp.set_yticklabels(np.linspace(-100,0,3),fontsize=14)
      sp_str = r'\displaystyle \log_{10}(%s)$' % label_str
      sp.set_xlabel(sp_str,fontsize=18)
      #str_title = r'$\displaystyle 1^{st}$ LEs'
      #sp.set_title(str_title,fontsize=18)
      
      plt.suptitle(r'Constant $\displaystyle %s=%d$' % (str_hPconst[hPconst],constVal),fontsize=20)
      if save:
	save_str = './pics/poisson/LE_const_%s=%d.pdf' % (hPconst,constVal)
	plt.savefig(save_str)
      plt.show()
  #return 0







####def scriptRate(steps,intervall,mode=None):
  ####nu_range = np.linspace(intervall[0],intervall[1],steps)
  
  ####if not mode:
    ####mode = raw_input('(r)ead or (s)imulate data? ')
  
  ####if mode == 's':
    ####for nu_stp in nu_range:
      ####simu = fluxtubesize(mode='t',topo='n')
      ####print ('\nNow simulating network with recurrent firing rate=%5.3f' % nu_stp).rstrip('0').rstrip('.')
      ####simu.const['rateWnt'] = nu_stp
      ####break_it = simu.setup(check=0)
  
      ####if break_it == 1:
	####erase(simu.Path,simu.Hash)
      ####elif break_it == 0:
	####simu.simulation()
	
  ####elif mode == 'r':
    ####sessionPath = 'data/'
    ####eft = np.zeros(steps)
    ####var_eft = np.zeros(steps)
    
    ####idx = 0
    ####for nu in nu_range:
      #####if nu%1:
      ####str_nu = ('rateWnt=%5.3f_' % nu).rstrip('0').rstrip('.')
      #####else:
	#####str_nu = 'hPrate=%d_' % nu
      ####for dat in os.listdir(sessionPath + 'info/'):
	####if 'normal' in dat:
	  ####if str_nu in dat:
	    ####fileName=sessionPath + 'info/' + dat
	    ####print 'Reading %s.' % (str_nu.split('_')[0])
	    ####try:
	      ####eft[idx],var_eft[idx],numsim = analyze('read',fileName_in=fileName,scriptName_in='normal',call=1)
	      ####print "hey"
	    ####except:
	      ####print 'No data for %s.' % (str_nu)
	      ####eft[idx] = None
      ####idx += 1
  
    ####for i in range(steps):
      ####if not eft[i]:
	####eft[i] = None
    
    ####nu_plt, = plt.plot(nu_range,eft,label='no poisson')
    ####plt.errorbar(nu_range,eft,yerr=var_eft,ecolor='r',elinewidth=2,fmt=None)
    ####plt.xlabel('rateWnt')
    ####plt.ylabel('eft')
    ####plt.legend(['poisson','no poisson'])
    ####plt.xscale('log')
    ####plt.xlim([0,max(nu_range)])
    ####plt.show(block=False) 
  
    ####return
