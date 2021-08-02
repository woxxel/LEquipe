import os, imp, hashlib, inspect, time, sys, shutil
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.special as spec
from numpy import complex
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

assert os.getcwd().split('/')[-1] == 'Program', 'Please run the program in the directory of LEquipe to properly set up the paths.'

imp.load_source('read_write_code','Code/py_code/read_write.py')
imp.load_source('file_management_code','Code/py_code/file_management.py')
imp.load_source('support_code', 'Code/py_code/support.py')
#imp.load_source('analysis_code', 'py_code/analysis.py')
from read_write_code import *
from file_management_code import *
from support_code import *
#from analysis_code import *

# Rainers Mail: Analyze several differences in kind of simulations:
#    perturbations in which representation? (what happens to threshold-crossing perturbations?)
#    simulation start at spike or time? (warmup stop when?) check spike times, 1st one when?
#    -> when is perturbation applied?
#    simulation stop at spike or time?
#
#
# next: understand stochastic driving neurons, implement (2 populations)

# connectivities / weights from learning networks -> stabilize?
#    -> run for some time, then take as connectivity
#    -> or implement learning in LEquipe to have it while simulating?


class simulation_data:
  
  def __init__(self,mode_input='t',infoName=None,net=None,mode='eft',TC=0.2):
    
    self.mode = mode
    
    self.cluster = {'scratch':0,'par':0,'q':0}		# clusterstuff
    
    self.Paras = {}					# bib to write simulation parameters to
    
    if self.mode == 'LE':
      # calculate Lypanunov Spectra
      self.TC = 2
      self.Paras['LyapunovExp'] = net.Const['N']	# Lyapunov Exponents
      self.Paras['LyapunovExpConvergence'] = 1
      if net.topo == 'p':
	self.Paras['subLyapunovExp'] = net.Const['N']
	self.Paras['subLyapunovExpConvergence'] = 1
      self.Paras['SWONS'] = 1
      self.Paras['ONstep'] = 10*float(net.Const['N'])/net.Const['K']
      self.Paras['pLE'] = 1
      
      #self.Paras['CLV'] = 1				# Covariant Lyapunov Exponents / local Lyapunov Exp
      #self.Paras['SWCLV'] = 1
      
    
    if self.mode == 'statistics':			# no warmup, thus all measures are made in one simulation
      # calculate cv, skewness, etc
      self.TC = 10.
      self.Paras['train'] = range(net.Const['N'])
      self.Paras['measures'] = 2
      self.Paras['measureTimes'] = np.linspace(0,self.TC,self.TC*1000.)
      self.Paras['synchrony'] = 1
      self.Paras['ISIneurons'] = range(net.Const['N'])			# neurons whose spike statistics are calculated
      self.Paras['ISIstats'] = 2					# moments that are being calculated
      self.Paras['ISIdecorr'] = range(net.Const['N'])					# ISI of potential decorrelation events are calculated for x neurons
    
    
    if self.mode == 'eft':
      # calculate mean flux tube radius
      self.TC = TC
      self.Paras['CDB'] = 1						# only valid for self.simulation
      self.Paras['measureTimes'] = np.linspace(0,TC,101)		# not valid for warmup
      #self.Paras['measures'] = 2
      #self.Paras['train'] = range(net.Const['N'])

      
    if self.mode == 'mosaic':
      # calculate "mosaic" cross-section of fluxtubes
      self.TC = 0.4
      self.Paras['train'] = range(net.Const['N'])
      self.Paras['measureTimes'] = np.linspace(0,TC,101)

    
    if self.mode == 'samples':
      # simulate sample trajectories to e.g. track "cloud" behaviour
      self.TC = 5-np.log10(net.Const['N'])
      self.Paras['train'] = range(net.Const['N'])
      self.Paras['measures'] = 2
      self.Paras['measureTimes'] = np.linspace(0,self.TC,1001)
      #self.Paras['measureSpikes'] = range(
      
    if self.mode == 'drive_cur':
      self.TC = 0
    #if self.mode == 'pert_spread':			# only one perturbed trajectory is calculated
      #self.TC = TC
      #self.Paras['train'] = range(net.Const['N'])
      #self.Paras['measures'] = 2
      #self.Paras['measureTimes'] = np.linspace(0,TC,1001)
      ##if 'plot' in mode_sim:
	##self.steps = {'ps':1, 'co':1,'in':1,'pd':1}
      ##else:
      #self.Paras['CDB'] = 1

      
    if mode_input == 'r':
      tmp_para, tmp_topo, self.steps, self.Hash = read_from_info(infoName)
      
    else:
      if self.mode == 'LE':
	self.steps = {'ps':1, 'co':1, 'in':1, 'pd':1}
      
      if self.mode == 'statistics':      
	self.steps = {'ps':1, 'co':10,'in':1,'pd':1}			# get 10 different sims
      
      if self.mode == 'eft':
	#self.steps = {'ps':11, 'co':5,'in':10,'pd':10}			# iteration steps
	self.steps = {'ps':21, 'co':5,'in':20,'pd':30}			# iteration steps
	#self.steps = {'ps':31, 'co':10,'in':25,'pd':55}		# iteration steps (dont change!)    
      
      if self.mode == 'mosaic':
	n=201				# decent value to get nice image quickly
	#if not (mode_input == 'r'):
	  #n = int(raw_input('How many samples across one dimension? '))	# number of bins across one dimension
	self.steps = {'ps':1, 'co':1,'in':1,'pd':n**2}
      
      if self.mode == 'samples':
	self.steps = {'ps':3, 'co':3,'in':3,'pd':5}
      
      if self.mode == 'drive_cur':
	self.steps = {'ps':1, 'co':3,'in':1,'pd':1}
      #if self.mode == 'pert_spread':
	#self.steps = {'ps':1, 'co':5,'in':10,'pd':20}
      
    
      self.steps['total'] = self.steps['pd']*self.steps['in']*self.steps['co']
    

class network(read_write):
  
  # Liederbuch
  def __init__(self,mode_input='t',infoName=None,alpha=None,hPVars=None,puppet=None,netConst={'N':1000,'K':100,'rateWnt':10,'J0':-1,'NeuronType':2}):
    
    self.Const = {}	# write constants of the network into this
    self.topoConst = {}
    self.puppet = 0
    
    if mode_input == 'r':
      self.infoName = infoName
      self.Const, self.topoConst, tmp_steps, self.Hash = read_from_info(infoName)
      if 'drive_rate' in self.topoConst.keys():
	self.topo = 'p'
	self.scriptName = 'poisson'
      elif 'drive_type' in self.topoConst.keys():
	self.topo = self.topoConst['drive_type']
	self.scriptName = self.topoConst['drive_type']
	self.puppet = 1
      else:
	self.topo = 'S'
	self.scriptName = 'SONET'
	
    else:
      if not all(np.isnan(alpha)):		# mode 'topo': n = normal, p = poisson, s = SONet
	self.topo = 'S'
	self.scriptName = 'SONET'
	self.topoConst['alpha_recip'] = alpha[0]	# 0:recip, 1:conv, 2:div, 3:chain
	self.topoConst['alpha_conv'] = alpha[1]
	self.topoConst['alpha_div'] = alpha[2]
	self.topoConst['alpha_chain'] = alpha[3]
      else:
	self.topoConst['alpha_recip'] = 0
	self.topoConst['alpha_conv'] = 0
	self.topoConst['alpha_div'] = 0
	self.topoConst['alpha_chain'] = 0

      if not all(np.isnan(hPVars)):
	self.topo = 'p'
	self.scriptName = 'poisson'
	self.topoConst['drive_rate'] = hPVars[0]
	self.topoConst['drive_cplg'] = hPVars[1]
      
      if type(puppet[0]) == str:
	self.topo = puppet[0][0]
	self.puppet = 1
	self.scriptName = puppet[0]
	self.topoConst['drive_type'] = puppet[0]
	for key in puppet[1].keys():
	  self.topoConst['drive_'+key] = puppet[1][key]

      for key in netConst.keys():
	self.Const[key] = netConst[key]
      
      self.Const['tauM'] = 10.
      
    
    
    self.rateWnt = self.Const['rateWnt']


  def setup(self,sim,perturbation=None,sv_file=[1,1,1,1,1,1],hide=1,call=0):
    
    # write to simulation dictionaries
    self.ParaNet, self.ParaSim = {}, {}
    self.break_it = 0
    
    for items in self.Const.items():
      self.ParaNet[items[0]] = items[1]
    
    self.perturbation = perturbation
    
    self.ParaNet['Iext'] = -(self.Const['J0']*self.rateWnt*np.sqrt(self.Const['K']))*self.Const['tauM']/1000.
	
    if self.topo == 'p':
      self.ParaNet['Iext'] -= min(self.Const['rateWnt']*self.Const['K'],self.topoConst['drive_rate'])*max(-10,self.topoConst['drive_cplg'])*self.Const['tauM']/1000.
      self.ParaNet['Iext'] = np.kron([self.ParaNet['Iext'],self.topoConst['drive_rate']],np.ones(self.Const['N']))
      
      self.ParaSim['rateWnt'] = 0				# unfix rate of whole network
      self.ParaSim['rateWntSubN'] = self.rateWnt
      self.ParaSim['subN'] = self.Const['N']
      self.ParaNet['N'] = self.Const['N']*2			# extend network by N poissonneurons
      if 'train' in sim.Paras.keys():
	  sim.Paras['train'] = range(self.ParaNet['N'])
	  
      if self.puppet:
	self.ParaNet['NeuronType'] = np.kron(np.array([self.Const['NeuronType'],7]),np.ones(self.Const['N']))
      else:
	self.ParaNet['NeuronType'] = np.kron((self.Const['NeuronType'],5),np.ones(self.Const['N']))
	self.ParaNet['poissonrate'] = np.kron([0,self.topoConst['drive_rate']],np.ones(self.Const['N']))
	self.initPoisson = np.random.randint(1,2**32-1,size=self.Const['N'])	# fix seed for all simulations
    else:
      self.ParaSim['rateWnt'] = self.rateWnt
      
    if not (sim.mode == 'eft'):				# so that directory name includes script name
      self.scriptName += '_' + sim.mode
    self.data_structure(sim,call=call)				# set data structure (names, files, folders)
    if self.break_it or call:					# if something went wrong
      return self.break_it
   
    if sim.mode == 'eft':
      print 'Starting ratefinding and network warmups...'
    
    # set up storing file
    datinit_sv = netcdf.netcdf_file(self.Path['inpara_links'],'w')
    datinit_sv.createDimension('connectivities',sim.steps['co'])
    datinit_sv.createDimension('initial conditions',sim.steps['in'])
    datinit_sv.createDimension('Simulations',sim.steps['co']*sim.steps['in'])
    datinit_sv.createDimension('links',6)
    datinit_sv.createDimension('linksize',len(self.Path['inpara'])+51)
    saveLinks = datinit_sv.createVariable('initial data','S1', ('connectivities','initial conditions','links','linksize'))
    
    sim.steps['ps_ct'], sim.steps['pd_ct'] = None, None
    
    time_start = time.time()
    
    # start getting different topologies with their respective external current
    for sim.steps['co_ct'] in range(sim.steps['co']):
      self.ParaTopo = {}
      self.create_topology(hide)
      if self.break_it:			# if topology could not be set up
	return self.break_it
      
      self.writeTopo()
      
      self.ParaSim['TR'] = 2		# rate finding only when connectivity changes
      
      self.ParaSim['saveFinalState'] = 1
      
      
      
      
	
      for sim.steps['in_ct'] in range(sim.steps['in']):	# change initial conditions
	#print self.ParaNet['Iext'][0]
	
	if not (sim.mode == 'drive_cur'):
	  self.ParaSim['TW'] = 1
	
	if sim.mode in ['LE','statistics']:	# just one simulation
	  for items in sim.Paras.items():
	    self.ParaSim[items[0]] = items[1]
	  self.ParaSim['TC'] = sim.TC
	#else:					# if more simulations needed, this is only warmup
	self.ParaSim['saveFinalState'] = 1
	  
	if (sim.steps['co'] > 1) or (sim.steps['in'] > 1):
	  print_status(sim,time_start)
	
	if self.puppet:
	  try:
	    time_limit = self.ParaSim['TR'][0]+self.ParaSim['TW']
	  except:
	    time_limit = self.ParaSim['TR']+self.ParaSim['TW']
	  if sim.mode in ['LE','statistics']:
	    time_limit += sim.TC
	  self.create_puppetTrain(time_limit)
	
	# set Paths and create .nc files for simulations
	self.writeNet()
	self.writePuppet()
	self.writeSim()
	
	#print self.ParaSim.items()
	HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo'],self.Path['Puppet']])).hexdigest()	# construct and ...
	self.Path['Out'] = self.Path['results'] + 'DataIni-' + HashDataOut + '.nc'	# ... assign Out-Name
	#print self.Path['Out']
	# run of simulation to either calculate specified values or ratefinding and warmup
	run_script(self.Path,sim.cluster,CDB=0,hide=hide)
	
	# remove unnecessary files to save disk space (if specified so)
	if not sv_file[0]:
	  os.remove(self.Path['Net'])
	if not sv_file[1]:
	  os.remove(self.Path['Sim'])
	
	PathIni = self.Path['Out']
	
	Data = readDataOut(self.Path['Out'])	# read results
	
	if (not sim.steps['in_ct']) and (sim.mode in ['eft','mosaic','samples']):	# after ratefinding, only warump is needed for this topology
	  self.ParaSim['TR'] = 0	# disable ratefinding
	  
	  if self.topo == 'p' and not self.puppet:	# set external drive for different initial conditions of this topology to found value
	    self.ParaNet['Iext'] = np.kron([Data['finalCurrents'][0],self.topoConst['drive_rate']],np.ones(self.Const['N']))
	    self.ParaSim['rateWntSub'] = 0
	  else:
	    self.ParaNet['Iext'] = Data['finalCurrents'][0]
	
	if sim.mode in ['eft','pert_spread','mosaic','samples']:	# then run first simulation to have reference trajectory

	  trash = self.ParaSim.pop('TW')
	  
	  for items in sim.Paras.items():
	    self.ParaSim[items[0]] = items[1]
	  
	  self.ParaSim['CDB'] = 0	# no comparing to ref trajectory here
	  if 'CDB' in sim.Paras.keys() or (sim.mode == 'samples'):
	    self.ParaSim['measures'] = 2
	    trash = self.ParaSim.pop('saveFinalState')
	    
	  self.ParaSim['TC'] = sim.TC
	  
	  if self.topo == 'p' and not self.puppet:	# get warmed up initial conditions with fixed poisson seed for all simulations
	    self.ParaNet['init'] = np.append(Data['finalStates'][:self.Const['N']],self.initPoisson)
	  else:
	    self.ParaNet['init'] = Data['finalStates']
	  
	  if self.puppet:
	    self.create_puppetTrain(self.ParaSim['TC'])
	  else:
	    self.Path['Puppet'] = ''
	  
	  self.writeNet()
	  self.writePuppet()
	  self.writeSim()
	  
	  HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()
	  self.Path['Out'] = self.Path['results'] + 'RefTraj-' + HashDataOut + '.nc'

	  run_script(self.Path,sim.cluster,CDB=0,hide=hide)		## run the C++ simulation for calculation of unperturbed trajectory
	  
	  trash = self.ParaSim.pop('measures')
	  trash = self.ParaSim.pop('measureTimes')
	  trash = self.ParaSim.pop('TC')
	  
	  # remove unnecessary files to save disk space (if specified so)
	  if (not sv_file[0]) and not (self.puppet):		# puppetTrain saved in here
	    os.remove(self.Path['Net'])
	  if not sv_file[1]:
	    os.remove(self.Path['Sim'])
	    
	trash = self.ParaNet.pop('init')			# remove initial conditions to generate new ones
	  
	# save links in ncfile for setting references
	saveLinks[sim.steps['co_ct'],sim.steps['in_ct'],:,:] = np.array([list(self.Path['Net']),list(self.Path['Sim']),list(self.Path['Topo']),list(PathIni),list(self.Path['Puppet']),list(self.Path['Out'])])
      
    datinit_sv.close()
    if sim.mode == 'eft':
      print ']'
      print 'time elapsed: %g sec\n' % (time.time()-time_start)
    
    # set parameters for simulation
    self.ParaSim['TC'] = sim.TC
    for items in sim.Paras.items():
      self.ParaSim[items[0]] = items[1]
    if 'CDB' in sim.Paras.keys():
      self.ParaSim['saveFinalState'] = 1
    
    if sim.mode in ['eft','pert_spread','mosaic','samples']:
      self.get_pert_range(sim)
    
    return self.break_it
    
    
  
# ----------------------------------- simulation stuff -----------------------------------------      
  
  def simulation(self,sim,sv_file=[0,0,1,1,1],check=0,hide=1):
    
    # initiate save arrays
    self.initDataSv(sim)
    
    if sim.mode == 'mosaic':  # generate one random vector as direction
      pert_vect = np.zeros((2,self.Const['N']))
      i=0
      while i<2:
	pert_vect_tmp = OrthoVector(np.random.uniform(-1,1,self.Const['N']))
	pert_vect_tmp /= np.linalg.norm(pert_vect_tmp)
	pert_vect[i] = pert_vect_tmp
	if (not i) or (i==1 and np.dot(pert_vect[0],pert_vect[1]) < 10**(-5)):
	  i += 1  
    else:
      pert_vect = None
    
    datinit_rd = netcdf.netcdf_file(self.Path['inpara_links'],'r')	# read input data from this
    
    idx_pert = 0
    
    print 'Now the real simulations start'
    
    for sim.steps['ps_ct'] in range(sim.steps['ps']):
      # setup time and save-stuff for this perturbation size
      time_start = time.time()
      datout_sv,svLinks = sv_links(self.Path,sim.steps,sim.steps['ps_ct'],self.pert_size[idx_pert])
      
      for sim.steps['co_ct'] in range(sim.steps['co']):
	
	for sim.steps['in_ct'] in range(sim.steps['in']):
	  
	  [self.Path['Net'],self.Path['Sim'],self.Path['Topo'],self.Path['Data'],self.Path['Puppet'],self.Path['Ref']] = [''.join(fileNames) for fileNames in datinit_rd.variables['initial data'][sim.steps['co_ct'],sim.steps['in_ct'],:,:]]
	  
	  Data = readDataOut(self.Path['Data'])
	  
	  if self.topo == 'p' and not self.puppet:
	    init = np.append(Data['finalStates'][:self.Const['N']],self.initPoisson)
	    self.ParaNet['Iext'] = np.kron([Data['finalCurrents'][0],self.topoConst['drive_rate']],np.ones(self.Const['N']))
	  else:
	    init = Data['finalStates']
	    self.ParaNet['Iext'] = Data['finalCurrents'][0]
	  
	  for sim.steps['pd_ct'] in range(sim.steps['pd']):

	    print_status(sim,time_start,self.pert_size[idx_pert])
	    pert_fail = self.perturb_state(sim,init,self.pert_size[idx_pert],pert_vect)
	    
	    # set Paths and create .nc files for simulations (need to do that before breaking, otherwise problems with reading
	    self.writeNet()
	    self.writeSim()
	    
	    HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()	# construct and ...
	    self.Path['Out'] = self.Path['results'] + 'DataOut-' + HashDataOut + '.nc'			# ... assign Out-Name
	    
	    svLinks[sim.steps['ps_ct'],sim.steps['co_ct'],sim.steps['in_ct'],sim.steps['pd_ct'],:] = [list(self.Path['Net']),list(self.Path['Sim']),list(self.Path['Topo']),list(self.Path['Out']),list(self.Path['Ref'])]	# archive file names

	    if not pert_fail:
	      run_script(self.Path,sim.cluster,CDB=self.ParaSim['CDB'][0],hide=hide)	## run the C++ simulation
	    
	    trash = self.ParaNet['init']
	    
	    # remove unnecessary files to save disk space (if specified so)
	    if not sv_file[0]:
	      os.remove(self.Path['Net'])
	    if not sv_file[1]:
	      os.remove(self.Path['Sim'])
	    #queue_control(sim.steps,self.Const['N'],max_sim=1000,wait_time=5)
	    # end of pd_stp
	  # end of in_stp
	# end of co_stp
      print ']'
      idx_pert+=1	# maybe doesn't work when running "mosaic"
      datout_sv.close()
      if sim.mode in ['eft','pert_spread','mosaic']:
	# directly read results to promote fast finding
	self.prob_read(sim,steps_rd=[sim.steps['ps_ct']])
	
	
	print 'time elapsed: %g sec\n' % (time.time()-time_start)
	
	if sim.mode == 'eft':
	  p = self.DataSv['p'][sim.steps['ps_ct']]
	  
	  if sim.mode == 'eft':
	    if ((p < 0.1) or (p > 0.7)):	# fast forward if far from 1-1/e
	      idx_pert+=3
	    if ((idx_pert >= sim.steps['ps']) or (p==1) or (not self.DataSv['correct_ct'][sim.steps['ps_ct']])):# or (self.pert_size[idx_pert] > 1)):
	      break
	
    datinit_rd.close()
    
    if sim.mode == 'eft':
      self.fluxSimSv(sim.steps)
      for dat in os.listdir(self.Path['results']):	# remove all data but the processed one
	if not 'processed' in dat:
	  os.remove(self.Path['results'] + dat)
    
    
    #------------------------------ support code -----------------------------
    
  def initDataSv(self,sim):
    self.DataSv = {}
    self.DataSv['pert_size'] = np.zeros(sim.steps['ps'])
    self.DataSv['p'] = np.zeros(sim.steps['ps'])
    self.DataSv['correct_ct'] = np.zeros(sim.steps['ps'])
    self.DataSv['CDBTimes'] = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    self.DataSv['div_track'] =  np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    self.DataSv['distances'] = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
  
  
  
  def perturb_state(self,sim,init,pert_size,pert_vect):
      # add perturbation of single neurons?
      # add perturbation in voltage instead of phases?
    
    if sim.mode == 'mosaic':
      
      n = int(np.sqrt(sim.steps['pd']))
      #construct initial conditions in square around unperturbed init cond
      frac0 = 2*(sim.steps['pd_ct']%n-(n-1)/2)/float(n-1)
      frac1 = 2*(sim.steps['pd_ct']/n-(n-1)/2)/float(n-1)
      init_add = frac0*pert_size*pert_vect[0] + frac1*pert_size*pert_vect[1]
      
      if self.topo == 'p' and not self.puppet:
	init_add = np.append(init_add,np.zeros(self.Const['N']))
      self.ParaNet['init'] = init + init_add
      
      return 0
      
    #elif sim.mode == 'samples':	# needed?
      #self.ParaNet['init'] = init + pert_size*pert_vect
      #return 0
    elif 'poisson_pert' in self.scriptName:
      
      ncidPuppet = netcdf.netcdf_file(self.Path['Puppet'],'r')
      puppetTrain = ncidPuppet.variables['puppetTrain'][:]
      self.puppetTrainSz = ncidPuppet.dimensions['puppetTrainSz']
      ncidPuppet.close()
      
      sigma = 0.1
      puppetTrain_tmp = puppetTrain + np.random.normal(scale=sigma,size=(self.Const['N'],self.puppetTrainSz))
      np.copyto(puppetTrain_tmp,puppetTrain,where=puppetTrain_tmp<0)
      
      self.ParaNet['puppetTrain'] = puppetTrain_tmp
      self.ParaNet['init'] = init
      return 0
      
    else:
      for i in range(10):	# loop to ensure not perturbing over threshold
	pert_vect = OrthoVector(np.random.uniform(-1,1,self.Const['N']))	# get orthogonal vector
	
	pert_vect *= pert_size/np.linalg.norm(pert_vect)				# ... normalize and scale it ...
	
	if self.topo == 'p':	# when poisson neurons present
	  pert_vect = np.append(pert_vect,np.zeros(self.Const['N']))		# (... extend by 0s for poisson neurons ...)
	
	#elif self.mode == 'single_neuron_perturbation':				# this mode is not at all added, yet
	  #pert_vect = np.zeros(self.Const['N'])
	  #pert_vect[np.random.randint(self.Const['N'])] = pert_size
	init_pert = init + pert_vect						# ...and add it to initial phase point

	if max(init_pert[:self.Const['N']]) < 1:
	  self.ParaNet['init'] = init_pert
	  return 0
    return 1	# if no proper perturbation could be found
  
  
  
  def get_pert_range(self,sim):
    self.eft_exp()
    
    if sim.steps['ps'] == 1:	# perturbation of expected fluxtube size (should look at others as well?
      self.pert_size = np.array([self.perturbation])
    else:
      if sim.mode == 'samples':
	self.pert_size = self.perturbation*10**(np.linspace(-1,0,sim.steps['ps']))
      else:
	self.pert_size = self.eft_exp*10**np.linspace(-2,2,sim.steps['ps'])	# pertSizes logscaled around eft_exp
  
  
  def eft_exp(self):
    self.eft_exp = 1./(math.sqrt(self.Const['N']*self.Const['K'])*self.Const['rateWnt']*self.Const['tauM']/float(1000))		# expected tube size from paper as medium guess
  
  
  
  def calc_max_dist(self,I_0=None):
    
    # define and process inputs
    V_T = I_T = 1.
    V_R = 0.
    
    #nu_array = [1,2,5,10,20]
    
    plt.figure()
    #for j in range(5):
    #self.Const['rateWnt'] = nu_array[j]
    
    #if not I_0:		# can be taken from simulation after warmup
      #I_0 = I_T - self.Const['J0']*np.sqrt(self.Const['K'])*self.Const['rateWnt']*self.Const['tauM']/1000.
    #print I_0
      #if self.topo == 'p':
	#I_0 += - self.topoConst['drive_cplg']*self.topoConst['drive_rate']*self.Const['tauM']/1000.	
      #print I_0

    sigma2 = self.Const['J0']**2*self.Const['rateWnt']*self.Const['tauM']/1000.	# could also be taken, but easier like that ?!
    mu = I_T + np.sqrt(self.Const['K'])*(I_0 + self.Const['J0']*self.Const['rateWnt']*self.Const['tauM']/1000.)
    
    #if self.topo == 'p':
      #sigma2 += self.topoConst['drive_cplg']**2*self.topoConst['drive_rate']*self.Const['tauM']/1000.
      #mu += self.topoConst['drive_cplg']*self.topoConst['drive_rate']*self.Const['tauM']/1000.
    
    sigma = np.sqrt(sigma2)
    
    
    erfi = lambda z: complex(0.,-1.)*spec.erf(complex(0.,1.)*z)	#define imaginary error function
    
    steps = 1001
    V = np.linspace(-0.5,1.5,steps)
    
    erfi_tmp = np.zeros(steps)
    
    for i in range(steps):
      if V[i] < 0:
	#print erfi((V_R-mu)/sigma)
	erfi_tmp[i] = erfi((V_R-mu)/sigma)
      else:
	#print erfi((V[i]-mu)/sigma)
	erfi_tmp[i] = erfi((V[i]-mu)/sigma)
      #print erfi_tmp[i]
    
    P_V = np.sqrt(math.pi)*self.Const['rateWnt']/sigma*np.exp(-(V-mu)**2/sigma2)*(erfi((V_T-mu)/sigma)-erfi_tmp)
    
    #plt.figure()
    plt.plot(V,P_V)
    plt.ylim([0,5])
    plt.xlim([-0.5,1.5])
  
    plt.show(block=False)
    
    
  
  
  def create_topology(self,hide):
    
    N = self.Const['N']
    K = self.Const['K']
    
    run_str = './Code/SONETs/run_secorder '
    str_save_tmp = self.Path['inpara'] + 'results_topo.nc'
    #alpha_str = '%5.3f %5.3f %5.3f %5.3g' % (self.topoConst['alpha_recip'],self.topoConst['alpha_conv'],self.topoConst['alpha_div'],self.topoConst['alpha_chain'])
    seed = np.random.randint(0,2**32-1)
    str_paras = '%d %5.3f %5.3f %5.3f %5.3f %5.3f' % (N,K/float(N),self.topoConst['alpha_recip'],self.topoConst['alpha_conv'],self.topoConst['alpha_div'],self.topoConst['alpha_chain'])
    run_str +=  '%s %s %d' % (str_paras,str_save_tmp,seed)
    fileParas = '_'.join(str_paras.split(' '))
    
    found = 0
    i = 0
    while i <= 10:
      runIt(run_str,hide)		# get in- and out-degree from here if needed for sims
      
      try:
	# read in the results of topology creation
	ncid = netcdf.netcdf_file(str_save_tmp,'r')
	for item in ncid.variables.items():
	  try:
	    self.ParaTopo[item[0]] = item[1][:]
	  except:
	    self.ParaTopo[item[0]] = item[1].getValue()
	ncid.close()
	
	os.remove(str_save_tmp)
	for dat in os.listdir(self.Path['session']):
	  if fileParas in dat:
	    os.remove(self.Path['session'] + dat)
	found = 1
	break
      
      except:
	print "No topology found - try again!"
	i += 1
      if i >= 10:
	print "No topology could be found!"
	if not self.topoConst['alpha_recip'] and not self.topoConst['alpha_conv'] and not self.topoConst['alpha_div'] and not self.topoConst['alpha_chain']:
	  print "Creating random matrix without SONET-control (could yield unwanted in/out degree correlations!)"
	  graph = random_graph(N,K)
	  for item in graph.items():
	    self.ParaTopo[item[0]] = item[1]
	found = 1
	
    if (len(np.unique(np.array(self.ParaNet['J0']))) > 1) or (self.topo == 'p'):
	self.ParaTopo['J'] = np.ones(len(self.ParaTopo['postSyn']))*self.ParaNet['J0']/np.sqrt(K)
	if self.topo == 'p':
	  self.ParaTopo['J'] = np.append(self.ParaTopo['J'],np.ones(N)*self.topoConst['drive_cplg'])
	  preSyn_tmp = np.zeros(len(self.ParaTopo['J']))
	  # adjust in/out degree and synapses for poisson networks
	  idx_deg = 0
	  idx_add = 0
	  for inDeg in self.ParaTopo['inDegree']:
	    preSyn_tmp[idx_deg+idx_add:idx_deg+idx_add+inDeg] = self.ParaTopo['preSyn'][idx_deg:idx_deg+inDeg]
	    preSyn_tmp[idx_deg+idx_add+inDeg] = idx_add+self.Const['N']
	    idx_deg += inDeg
	    idx_add += 1
	  self.ParaTopo['preSyn'] = preSyn_tmp
	  self.ParaTopo['postSyn'] = np.append(self.ParaTopo['postSyn'],np.arange(self.Const['N']))
	  self.ParaTopo['inDegree'] = np.append(self.ParaTopo['inDegree']+1,np.zeros(self.Const['N']))
	  self.ParaTopo['outDegree'] = np.append(self.ParaTopo['outDegree'],np.ones(self.Const['N']))
    else:
      self.ParaTopo['J'] = self.ParaNet['J0']/np.sqrt(K)
  


    if not found:
      print "how?"
      self.break_it = 1
  
  
  def create_puppetTrain(self,time_limit):
    #print "time: %g" % time_limit
    if self.topoConst['drive_type'] == 'poisson':
      self.puppetTrainSz = 0
      puppets = self.Const['N']
      puppetTrain = []
      for n in range(puppets):
	puppetTrain.append([])
	dt = 0
	while dt <= time_limit*1000:	# time measured in ms
	  dt_add = -np.log(np.random.uniform(0,1))/self.topoConst['drive_rate']*1000
	  puppetTrain[n].append(dt_add)
	  dt += dt_add
	self.puppetTrainSz = max(self.puppetTrainSz,len(puppetTrain[n]))

      self.ParaNet['puppetTrain'] = (np.ones((puppets,self.puppetTrainSz)) + time_limit)*1000	# fill with spiketimes > TC
      
      for n in range(puppets):
	self.ParaNet['puppetTrain'][n][:len(puppetTrain[n])] = puppetTrain[n]

  
  # ---------------------------------------------- data structure stuff ----------------------------
  
  def data_structure(self,sim,call=0):
    
    #generate output data paths
    self.Path = {}
    
    if sim.cluster['scratch'] == 1:
      self.Path['session'] = '/scratch01/aschmidt/data/'	#path for storage on scratch
    else:
      self.Path['session'] = 'data/'			#path for local storage
    if (not os.path.exists(self.Path['session']) and not call):
      os.mkdir(self.Path['session'])
      
    self.Path['script'] = self.Path['session'] + self.scriptName + '/'
    if (not os.path.exists(self.Path['script']) and not call):
      print 'creating new directory for this script: %s' % self.Path['script']
      os.mkdir(self.Path['script'])

    self.Path['info'] = self.Path['session'] + 'info/'
    if (not os.path.exists(self.Path['info']) and not call):
      os.mkdir(self.Path['info'])
    
    self.InfoHash(sim,call)
    if self.break_it:
      return
	
    self.Path['inpara'] = self.Path['script'] + 'inparas_' + self.Hash + '/'
    if (not os.path.exists(self.Path['inpara']) and not call):
      os.mkdir(self.Path['inpara'])
    
    self.Path['inpara_links'] = self.Path['inpara'] + 'datInit.nc'
    
    self.Path['results'] = self.Path['script'] + 'results_' + self.Hash + '/'      
    if (not os.path.exists(self.Path['results']) and not call):
      os.mkdir(self.Path['results'])
    self.Path['results_links'] = self.Path['results'] + 'filenames.nc'
    
    self.Path['txt'] = self.Path['script'] + 'txt_' + self.Hash + '/' 
    if (not os.path.exists(self.Path['txt']) and not call):
      os.mkdir(self.Path['txt'])



  def InfoHash(self,sim,call=0):
    # creating name for identification of simulation in info file
    self.Hash_array = []
    self.break_it = 0
    #print self.Path
    if not call:
      para_str_name, para_str_write = self.Info_generate_lines(self.Const,'float','para')
      topo_str_name, topo_str_write = self.Info_generate_lines(self.topoConst,'float','topo')
      step_str_name, step_str_write = self.Info_generate_lines(sim.steps,'int','step')
      saves_str_name, saves_str_write = self.Info_generate_lines(sim.Paras,None,'save')
      self.Path['info_file'] = self.Path['info'] + 'SimInfo_%s%s%s%s.txt' % (self.scriptName,para_str_name,topo_str_name,step_str_name)
    else:
      self.Path['info_file'] = self.infoName
    
    if call:
      try:
	tmp_para, tmp_topo, tmp_step, self.Hash = read_from_info(self.Path['info_file'])
      except:
	self.break_it = 1
      return
    elif os.path.exists(self.Path['info_file']):
      print 'The simulation %s exists already. Skipping...' % self.Path['info_file']
      self.break_it = 2
      return
    else:
      self.Hash = hashlib.sha1(np.array(self.Hash_array)).hexdigest()
      print 'This simulation can be accessed by its Hash: "%s".' % self.Hash
      
    # write info into information file
    info_sv = open(self.Path['info_file'],'w')
    info_sv.write('# Simulation script %s, executed on %s with a total of %d simulations per perturbation size\n' % (self.scriptName,time.strftime('%c'),sim.steps['total']) + \
    '#\n' + \
    para_str_write + topo_str_write + step_str_write + saves_str_write + \
    '#\nHash is <%s>' % self.Hash)
    info_sv.close()


  def Info_generate_lines(self,data,data_type,str_name):
    str_write = str_name + ' Variables: \t<'
    str_name = '_' + str_name
    
    if data_type == 'int':	# avoiding '.' for integers (better readout, better name_construction)
      dat_print = '%d'
    elif data_type == 'float':
      dat_print = '%5.3f'
      
    for item in data.items():
      if not data_type:
	string = str(item[0])
      else:
	if type(item[1])==str:
	  string = '%s=%s' % (item[0],item[1])
	else:
	  string = '%s=%s' % (item[0],dat_print) % item[1]
	if not(data_type == 'int'):
	  string = string.rstrip('0').rstrip('.')
	self.Hash_array.append(item[1])
      str_write += '\t' + string
      str_name += '_' + string
      
    str_write += '>\n'
    
    return str_name, str_write


  def fluxSimSv(self,steps):
    print "Saving results..."
    
    for item in self.DataSv.items():
      self.DataSv[item[0]] = item[1][:steps['ps_ct']+1]
    
    processedPath = self.Path['results'] + 'results_processed.nc'
    if os.path.exists(processedPath):
      os.remove(processedPath)
    ncid = netcdf.netcdf_file(processedPath,'w')
    
    ncid.createDimension('one',1)
    ncid.createDimension('ps_steps',steps['ps_ct']+1)
    ncid.createDimension('co_steps',steps['co'])
    ncid.createDimension('in_steps',steps['in'])
    ncid.createDimension('pd_steps',steps['pd'])
    
    VarN = ncid.createVariable('N','i', ('one',))
    VarK = ncid.createVariable('K','i', ('one',))
    Varf = ncid.createVariable('f','i', ('one',))
    VartauM = ncid.createVariable('tauM','i', ('one',))
    VarPertSz = ncid.createVariable('pert_size', 'd', ('ps_steps',))
    VarProbs = ncid.createVariable('p','d', ('ps_steps',))
    VarDivTrack = ncid.createVariable('div_track','d',('ps_steps','co_steps','in_steps','pd_steps'))
    VarTimes = ncid.createVariable('CDBTimes','d',('ps_steps','co_steps','in_steps','pd_steps'))
    VarDist = ncid.createVariable('distances','d',('ps_steps','co_steps','in_steps','pd_steps'))
    VarCorr = ncid.createVariable('correct_ct','i',('ps_steps',))
    
    VarN[:] = self.Const['N']
    VarK[:] = self.Const['K']
    Varf[:] = self.rateWnt
    VartauM[:] = self.Const['tauM']
  
    VarPertSz[:] = self.DataSv['pert_size']
    VarProbs[:] = self.DataSv['p']
    VarDivTrack[:] = self.DataSv['div_track']
    VarTimes[:] = self.DataSv['CDBTimes']
    VarDist[:] = self.DataSv['distances']
    VarCorr[:] = self.DataSv['correct_ct']

    if steps['ps'] > 1:
      eft, var_eft, var_p = self.eft_calc()
      VarTubeRadius = ncid.createVariable('flux_tube_radius','d', ('one',))
      VarTubeRadius[:] = eft
    
    ncid.close()


  def prob_read(self,sim,steps_rd=None,errorMsg=0):
    
    sv_all = 0
    
    ncid = netcdf.netcdf_file(self.Path['results_links'],'r')
    pert_size = ncid.variables['pert_size'][:]
    
    if not steps_rd:			# if not specified, read all
      steps_rd = range(len(pert_size))
    ps_steps = len(steps_rd)
    
    if sim.mode == 'eft':
      div_track = np.zeros((ps_steps,sim.steps['co'],sim.steps['in'],sim.steps['pd']))
      dist = np.ones((ps_steps,sim.steps['co'],sim.steps['in'],sim.steps['pd']))
      times = np.zeros((ps_steps,sim.steps['co'],sim.steps['in'],sim.steps['pd']))
      div_track[:] = np.nan
    
    print 'reading %d files from list %s...' % (ps_steps*sim.steps['total'],self.Path['results_links'])
    
    idx_ps = 0
    idx_ct = 0
    if sim.mode == 'samples':
      idx_cut = self.Const['N']*self.Const['rateWnt']	# number of spikes in the recurrent neurons after about 1 sec
      idx_cut_drive = self.Const['N']*self.topoConst['drive_rate']
      
      self.analysis['trainTimeDrive'] = np.zeros((sim.steps['in']*sim.steps['co']*ps_steps,idx_cut_drive))
      self.analysis['trainNeuronDrive'] = np.zeros((sim.steps['in']*sim.steps['co']*ps_steps,idx_cut_drive))
      
      self.analysis['trainTime'] = np.zeros((2,sim.steps['total']*ps_steps,idx_cut))
      self.analysis['trainNeuron'] = np.zeros((2,sim.steps['total']*ps_steps,idx_cut))
      self.analysis['trainDeltaTime'] = np.zeros((sim.steps['total']*ps_steps,idx_cut))
    
    for ps_ct in steps_rd:
      
      for co_ct in range(sim.steps['co']):
	
	for in_ct in range(sim.steps['in']):
	  link_unperturbed = ''.join(ncid.variables['links'][ps_ct,co_ct,in_ct,0,4])
	  try:
	    Data_unpert = readDataOut(link_unperturbed,1)
	  except:
	    print 'Data is not here anymore! Reading stopped!'
	    ncid.close()
	    return
	  if sim.mode == 'mosaic':
	    ct_fluxtube = 0
	    FluxtubeState = np.zeros((sim.steps['pd'],self.Const['N']))	# initialize array of maximum possible number of fluxtubes
	    FluxtubeState[ct_fluxtube] = Data_unpert['finalStates'][:self.Const['N']]
	    
	    n = int(np.sqrt(sim.steps['pd']))
	    fluxtube = np.zeros((n,n))
	    fluxtube[:] = np.nan
	  
	  if sim.mode == 'samples':
	    mask = np.where(Data_unpert['trainNeuron'] >= self.Const['N'])[0]
	    self.analysis['trainTimeDrive'][idx_ct/sim.steps['pd']] = Data_unpert['trainTime'][mask][:idx_cut_drive]
	    self.analysis['trainNeuronDrive'][idx_ct/sim.steps['pd']] = Data_unpert['trainNeuron'][mask][:idx_cut_drive]
	  for pd_ct in range(sim.steps['pd']):
	    link_perturbed = ''.join(ncid.variables['links'][ps_ct,co_ct,in_ct,pd_ct,3])

	    try:
	      Data_pert = readDataOut(link_perturbed,1)
	      
	      if sim.mode == 'eft':
		if sim.Paras['CDB']:
		  idx_time = len(Data_pert['measureTimes'])-1
		  times[idx_ps][co_ct][in_ct][pd_ct] = Data_pert['measureTimes'][idx_time]
		  dist_vec = (Data_unpert['measurePhases'][idx_time] - Data_pert['finalStates'])[:self.Const['N']]
		else:
		  times[idx_ps][co_ct][in_ct][pd_ct] = np.nan
		  dist_vec = (Data_unpert['finalStates'] - Data_pert['finalStates'])[:self.Const['N']]

		dist[idx_ps][co_ct][in_ct][pd_ct] = np.linalg.norm(OrthoVector(dist_vec))	# distance orthogonal to trajectory
		
		if dist[idx_ps][co_ct][in_ct][pd_ct] >= 0.01:	# if not converged, treat as diverged
		  div_track[idx_ps][co_ct][in_ct][pd_ct] = 1
		else:
		  div_track[idx_ps][co_ct][in_ct][pd_ct] = 0
	      
	      elif sim.mode == 'mosaic':
		for idx_fluxtube in range(ct_fluxtube+1):
		  dist_vec = (FluxtubeState[idx_fluxtube] - Data_pert['finalStates'][:self.Const['N']])
		  dist = np.linalg.norm(OrthoVector(dist_vec))
		  if dist < 0.01:
		    fluxtube[pd_ct%n][pd_ct/n] = idx_fluxtube
		    break
		if np.isnan(fluxtube[pd_ct%n][pd_ct/n]):
		  ct_fluxtube += 1
		  FluxtubeState[ct_fluxtube] = Data_pert['finalStates'][:self.Const['N']]
		  fluxtube[pd_ct%n][pd_ct/n] = ct_fluxtube
		  
	      elif sim.mode == 'samples':
		print idx_ct
		
		if not 'measureTimes' in self.analysis.keys():
		  self.analysis['measureTimes'] = Data_unpert['measureTimes']
		  time_len = Data_unpert['measurePhases'].shape[0]
		  
		  ### keep commmented stuff for video creation!
		  if sv_all:
		    self.analysis['distanceVector'] = np.zeros((sim.steps['total']*ps_steps,time_len,self.Const['N']))
		    self.analysis['distanceOrthVector'] = np.zeros((sim.steps['total']*ps_steps,time_len,self.Const['N']))
		    self.analysis['distance2'] = np.zeros((sim.steps['total']*ps_steps,time_len))
		  else:
		    distanceOrthVector = np.zeros((time_len,self.Const['N']))
		  self.analysis['distanceOrth'] = np.zeros((sim.steps['total']*ps_steps,time_len))
		  #self.analysis['distance1'] = np.zeros((sim.steps['total']*ps_steps,time_len))
		  self.analysis['convergeVelocity'] = np.zeros((2,sim.steps['total']*ps_steps))
		  self.analysis['kinkTime'] = np.zeros(sim.steps['total']*ps_steps)
		  #self.analysis['convergeVelocity2'] = np.zeros((sim.steps['total']*ps_steps,2)
		  
		if sv_all:
		  self.analysis['distanceVector'][idx_ct] = (Data_unpert['measurePhases']-Data_pert['measurePhases']).transpose()[:self.Const['N']].transpose()
		  self.analysis['distance2'][idx_ct] = np.sqrt(np.sum(self.analysis['distanceVector'][idx_ct]**2,axis=1))
		  #self.analysis['distance1'][idx_ct] = np.sum(self.analysis['distanceVector'][idx_ct],axis=1)
		else:
		  distanceVector = (Data_unpert['measurePhases']-Data_pert['measurePhases']).transpose()[:self.Const['N']].transpose()
		
		for t in range(time_len):
		  if sv_all:
		    self.analysis['distanceOrthVector'][idx_ct][t] = OrthoVector(self.analysis['distanceVector'][idx_ct][t])
		    self.analysis['distanceOrth'][idx_ct] = np.sqrt(np.sum(self.analysis['distanceOrthVector']**2,axis=1))
		  else:
		    distanceOrthVector[t] = OrthoVector(distanceVector[t])
		    self.analysis['distanceOrth'][idx_ct] = np.sqrt(np.sum(distanceOrthVector**2,axis=1))
		    #self.analysis['distance1'][idx_ct] = np.sum(np.abs(distanceOrthVector),axis=1)
		
		if self.analysis['distanceOrth'][idx_ct][-1] < 0.01:# only get convergence speed for converged runs
		  fit_idx = self.get_convergence(time_len,idx_ct)  
		else:
		  self.analysis['convergeVelocity'][0][idx_ct] = np.nan
		  self.analysis['convergeVelocity'][1][idx_ct] = np.nan
		
		
		mask = np.where(Data_unpert['trainNeuron'] < self.Const['N'])[0]
		  
		self.analysis['trainTime'][0][idx_ct] = Data_unpert['trainTime'][mask][:idx_cut]
		self.analysis['trainNeuron'][0][idx_ct] = Data_unpert['trainNeuron'][mask][:idx_cut]
		
		self.analysis['trainTime'][1][idx_ct] = Data_pert['trainTime'][mask][:idx_cut]
		self.analysis['trainNeuron'][1][idx_ct] = Data_pert['trainNeuron'][mask][:idx_cut]
		
		# and now, get "cloud" of initial conditions
		
		#if self.analysis['convergeVelocity'][0][idx_ct]:
		assignSpikes(self.analysis,idx_ct,idx_cut)	# just until 1 sec
		
		if self.analysis['convergeVelocity'][1][idx_ct] > 0:
		  levs = range(160)
		  rb = mcolors.LinearSegmentedColormap.from_list(name='red_black',colors=[(0,(0,0,0)),(1,(1,0,0))],N=len(levs)-1,)
		
		  fig,axes = plt.subplots(nrows=2,ncols=1)
		  
		  axes[0].scatter(self.analysis['trainTime'][0][idx_ct],self.analysis['trainNeuron'][0][idx_ct],c=np.log(1+self.analysis['trainDeltaTime'][idx_ct]),cmap=rb,vmin=0,vmax=max(np.log(1+self.analysis['trainDeltaTime'][idx_ct]))/100,marker='D')
		  axes[0].set_xlim([0,1])#self.analysis['measureTimes'][2*fit_idx[1]-1]])
		  axes[0].set_ylim([0,np.max(self.analysis['trainNeuron'])])
		  
		  axes[1].plot(self.analysis['measureTimes'],self.analysis['distanceOrth'][idx_ct],'-k')
		  axes[1].plot(self.analysis['measureTimes'][fit_idx],self.analysis['distanceOrth'][idx_ct][fit_idx],'or')
		  axes[1].set_xlim([0,1])#self.analysis['measureTimes'][2*fit_idx[1]-1]])
		  #axes[1].set_xlim([0,self.analysis['measureTimes'][-1]])
		  axes[1].set_ylim([10**(-16),10**1])
		  axes[1].set_yscale('log')
		  
		  plt.show()
		
		#return
	      # clean after reading
	    
	      os.remove(link_perturbed)
	    
	    except:
	      if errorMsg:
		print linkName + ' does not exist. Something seems to have gone wrong'
	      
	    idx_ct += 1
	    
      
      if sim.mode == 'eft':
	mask = np.invert(np.isnan(div_track[idx_ps]))
	correct_ct = np.sum(mask)
	
	if correct_ct:
	  p = np.sum(div_track[idx_ps][mask])/correct_ct
	else:
	  p = 0
	print p
	
	if sim.steps['ps'] > 1:
	  print 'probability of divergence at perturbation size ps = %5.3g is %5.3g (%d measurements)' % (pert_size[ps_ct],p,correct_ct)
      
	self.DataSv['pert_size'][ps_ct] = pert_size[ps_ct]
	self.DataSv['div_track'][ps_ct] = div_track[idx_ps]
	self.DataSv['p'][ps_ct] = p
	self.DataSv['CDBTimes'][ps_ct] = times[idx_ps]
	self.DataSv['distances'][ps_ct] = dist[idx_ps]
	self.DataSv['correct_ct'][ps_ct] = correct_ct
	
	idx_ps += 1
    
    if sim.mode == 'mosaic':
      
      ncidwrite = netcdf.netcdf_file(self.Path['results'] + 'results_processed.nc','w')
      ncidwrite.createDimension('n',n)
      ncidwrite.createDimension('one',1)
      ncidVarFluxSv = ncidwrite.createVariable('fluxtubes','i',('n','n'))
      ncidVarPertSz = ncidwrite.createVariable('pert_size','d',('one',))
      ncidVarFluxSv[:] = fluxtube
      ncidVarPertSz[:] = pert_size
      ncidwrite.close()
    
    elif sim.mode == 'samples' and ('convergeVelocity' in self.analysis.keys()):
      # remove these - too costly!!
      
      self.writeSamples(sim,time_len,ps_steps)
      
      for dat in os.listdir(self.Path['results']):
	if not ('results' in dat) and not ('filenames' in dat):
	  os.remove(self.Path['results'] + dat)
      #self.analysis.pop('distance2')
      #self.analysis.pop('distanceVector')
      #self.analysis.pop('distanceOrthVector')

    ncid.close()
  
  
  
  def writeSamples(self,sim,time_len,ps_steps):
    print self.Path['results'] + 'results_processed.nc'
    ncidwrite = netcdf.netcdf_file(self.Path['results'] + 'results_processed.nc','w')
    
    ncidwrite.createDimension('time',time_len)
    ncidwrite.createDimension('N',self.Const['N'])
    ncidwrite.createDimension('steps',sim.steps['total']*ps_steps)
    ncidwrite.createDimension('two',2)
    ncidwrite.createDimension('spikes',self.analysis['trainTime'].shape[-1])
    ncidwrite.createDimension('drive_steps',self.analysis['trainNeuronDrive'].shape[0])
    
    ncidVarTime = ncidwrite.createVariable('measureTimes','d',('time',))
    ncidVarDistOrth = ncidwrite.createVariable('distanceOrth','d',('steps','time'))
    ncidVarConvVel1 = ncidwrite.createVariable('convergeVelocity','d',('two','steps'))
    #ncidVarConvVel2 = ncidwrite.createVariable('convergeVelocity2','d',('steps',))
    ncidVarTrainTime = ncidwrite.createVariable('trainTime','d',('two','steps','spikes'))
    ncidVarTrainNeuron = ncidwrite.createVariable('trainNeuron','i',('two','steps','spikes'))
    ncidVardeltaTime = ncidwrite.createVariable('trainDeltaTime','d',('steps','spikes'))
    ncidVarTrainNeuronDrive = ncidwrite.createVariable('trainNeuronDrive','i',('drive_steps','spikes'))
    ncidVarTrainTimeDrive = ncidwrite.createVariable('trainTimeDrive','d',('drive_steps','spikes'))
    ncidVarKink = ncidwrite.createVariable('kinkTime','d',('steps',))
    
    ncidVarTime[:] = self.analysis['measureTimes']
    ncidVarDistOrth[:] = self.analysis['distanceOrth']
    ncidVarConvVel1[:] = self.analysis['convergeVelocity']
    #ncidVarConvVel2[:] = self.analysis['convergeVelocity2']
    ncidVarTrainTime[:] = self.analysis['trainTime']
    ncidVarTrainNeuron[:] = self.analysis['trainNeuron']
    ncidVardeltaTime[:] = self.analysis['trainDeltaTime']
    ncidVarTrainNeuronDrive[:] = self.analysis['trainNeuronDrive']
    ncidVarTrainTimeDrive[:] = self.analysis['trainTimeDrive']
    ncidVarKink[:] = self.analysis['kinkTime']
    
    #ncidVarDist1 = ncidwrite.createVariable('distance1','d',('steps','time'))
    #ncidVarDist1[:] = self.analysis['distance1']
    
    #ncidVarDistVec = ncidwrite.createVariable('distanceVector','d',('steps','time','N'))
    #ncidVarDist2 = ncidwrite.createVariable('distance2','d',('steps','time'))
    
    #ncidVarDistOrthVec = ncidwrite.createVariable('distanceOrthVector','d',('steps','time','N'))
    #ncidVarDistVec[:] = self.analysis['distanceVector']
    #ncidVarDist2[:] = self.analysis['distance2']
    #ncidVarDistOrthVec[:] = self.analysis['distanceOrthVector']
      
    ncidwrite.close()
    

  def get_convergence(self,time_len,idx_ct):
    
    if all(self.analysis['distanceOrth'][idx_ct] == 0):
      self.analysis['convergeVelocity'][0][idx_ct] = np.nan
      self.analysis['convergeVelocity'][1][idx_ct] = np.nan
      self.analysis['kinkTime'][idx_ct] = np.nan
      return

    min_idx_search = np.where(self.analysis['distanceOrth'][idx_ct] < 10**(-10))[0]
    
    if not len(min_idx_search):
      idx_est = time_len - 1
    else:
      idx_est = min_idx_search[0]
    
    m2,b2,idx_est2 = get_slope(self.analysis['measureTimes'],np.log(self.analysis['distanceOrth'][idx_ct]),idx_est,0)
    self.analysis['convergeVelocity'][1][idx_ct] = m2
    
    if self.analysis['measureTimes'][idx_est2] > 0.05:
      m1,b1,idx_est1 = get_slope(self.analysis['measureTimes'],np.log(self.analysis['distanceOrth'][idx_ct]),0,idx_est2)    
      self.analysis['convergeVelocity'][0][idx_ct] = m1
      
      # get kink
      x = np.linspace(0,self.analysis['measureTimes'].max(),1000)
      f1 = m1*x+b1
      f2 = m2*x+b2
      deltaf = np.abs(f1-f2)
      kink_idx = np.where(deltaf == min(deltaf))[0]
      #print kink_idx
      self.analysis['kinkTime'][idx_ct] = self.analysis['measureTimes'][kink_idx]
      #print self.analysis['kinkTime'][idx_ct]
    else:
      self.analysis['convergeVelocity'][0][idx_ct] = np.nan
      self.analysis['kinkTime'][idx_ct] = np.nan
      idx_est1 = 0
    
    return [idx_est1,idx_est2]
    
  
  def eft_calc(self,var_p=None):	# add proper error estimation (from eft calcs for different steps)
  
    if all(self.DataSv['p']==0):	# add here: maximum distance, calculated by c++ solver
      eft = np.nan
      var_eft = np.nan
      var_p = np.zeros(len(self.DataSv['p']))
      
    else:
      if not type(var_p)==np.ndarray:	# rather take variance error
	var_p = np.sqrt(self.DataSv['correct_ct']*self.DataSv['p']*(1-self.DataSv['p']))/self.DataSv['correct_ct']

      try:
	idx = np.where(self.DataSv['p'] > 0.9)[0][0]
      except:
	idx = len(self.DataSv['p'])
      mask = np.where(self.DataSv['p'] > 0)
      p_fit = self.DataSv['p'][mask][:idx]
      if len(p_fit):
	print p_fit
      else:
	print "All stable!"

      try:
	m, b, r_value, p_value, std_err = stats.linregress(np.log(self.DataSv['pert_size'][mask][:idx]),np.log(p_fit))
	
	print m, std_err
	eft = np.exp((np.log(1-1./np.exp(1)) - b)/m)
	var_eft = (np.log(1-1./np.exp(1)) - b)/(m**2)*np.exp((np.log(1-1./np.exp(1)) - b)/m)*std_err		# get proper error for eft! only slope-error taken into account here
      #idx = np.where(self.DataSv['p'] - (1-1/np.exp(1)) < 0)[0].max()
      
      #try:
	#delta_m, delta_p = self.DataSv['pert_size'][idx], self.DataSv['pert_size'][idx+1]
	#p_m, p_p = self.DataSv['p'][idx], self.DataSv['p'][idx+1]
	#var_p_m, var_p_p = var_p[idx], var_p[idx+1]
	#gamma_p, gamma_m = (1-1/np.exp(1)-p_p)/(p_p-p_m), (1-1/np.exp(1)-p_m)/(p_p-p_m)
	
	#eft = delta_p**gamma_m/(delta_m**gamma_p)
	#var_eft = np.sqrt(((eft*np.log(delta_p/delta_m))/(p_p-p_m))**2*((var_p_m*gamma_p)**2+(var_p_p*gamma_m)**2)) # why delta_m/delta_p? signmistake somewhere yields delta_m*delta_p, but this gives realistic results.
	
      except:
	eft = None
	var_eft = None
	print 'perturbation sizes are not in range of simulated data'
      
    return eft, var_eft, var_p
  
  
def get_slope(x,y,fix_idx,test_idx):
  ### rather iterate over errors, and change idx_cut accordingly - break, when precision high enough
  
  first_idx = min(test_idx,fix_idx)
  final_idx = max(test_idx,fix_idx)
  sign = np.sign(test_idx-fix_idx)
  
  #get first estimate:
  steps_est = np.where(x>(x[fix_idx]+sign*0.1))[0][0] - fix_idx
  
  m_est, b_est, r_value, p_value, std_err = stats.linregress(x[fix_idx:fix_idx+steps_est:sign],y[fix_idx:fix_idx+steps_est:sign])
    
  #plt_range = range(first_idx,final_idx)
  
  fit_err = np.zeros(final_idx)
  min_err = 1
  min_ct = 0
  min_idx = 0
  j = test_idx
  
  while ((final_idx-first_idx) > 10):
    if sign > 0:
      final_idx -= 1
      j -= 1
    else:
      first_idx += 1
      j += 1
    m, b, r_value, p_value, std_err = stats.linregress(x[first_idx:final_idx],y[first_idx:final_idx])
    #get error of fit
    
    fit_curve = m*x+b
    fit_err = np.var((fit_curve-y)[first_idx:final_idx])/(final_idx-first_idx)**2	# should cover the most idx as possible -> therefor var/N
    
    fit_err += 2*fit_err*((m-m_est)/m_est)**2	# quadratic bias towards expected slope value
    #fit_err[j] += fit_err[j]*((b-b_est)/b_est)**2	# linear bias towards intercept
    if fit_err < min_err:
      min_err = fit_err
      min_idx = j
      min_ct = 0
    else:
      min_ct += 1
    if min_ct == 10:
      break
  #print m
  #min_idx = np.where(fit_err == np.min(fit_err[j:test_idx:sign]))[0][0]
  #plt.figure()
  #plt.plot(plt_range,fit_err)
  #plt.plot(min_idx,fit_err[min_idx],'or')
  #plt.show(block=False)
  
  return m, b, min_idx

  
  
def sv_links(Path,steps,ps_ct,pert):
  
  if os.path.exists(Path['results_links']):
    tmp_datout = netcdf.netcdf_file(Path['results_links'],'r')
    tmp_links = tmp_datout.variables['links'][:]
    tmp_pert = tmp_datout.variables['pert_size'][:]
    #tmp_time = tmp_datout.variables['time'][:]
    tmp_datout.close()
  try:
    os.remove(Path['results_links'])
  except:
    print "start!"

  datout_sv = netcdf.netcdf_file(Path['results_links'],'w')	# save filenames to this

  datout_sv.createDimension('perturbation sizes',ps_ct+1)
  datout_sv.createDimension('connectivities',steps['co'])
  datout_sv.createDimension('initial conditions',steps['in'])
  datout_sv.createDimension('perturbation directions',steps['pd'])
  datout_sv.createDimension('Simulations',steps['co']*steps['in']*steps['pd'])
  datout_sv.createDimension('links',5)
  datout_sv.createDimension('linksizes',len(Path['results'])+51)	#hashnames are 51 long
  
  svLinks = datout_sv.createVariable('links','S1',('perturbation sizes','connectivities','initial conditions','perturbation directions','links','linksizes'))
  svPert = datout_sv.createVariable('pert_size', 'f4', ('perturbation sizes',))
  
  if ps_ct:
    svLinks[:ps_ct] = tmp_links
    svPert[:ps_ct] = tmp_pert
    
  svPert[-1] = pert
  
  return datout_sv, svLinks
  
    
def OrthoVector(vector):
  N = len(vector)
  norm_traj = np.ones(N)/np.linalg.norm(np.ones(N))	# norm vector of trajectory
  vector -= np.dot(vector,norm_traj)*norm_traj		# ... substract projection on trajectory
  return vector


  

def analyze(sim_mode,net=None,sim=None,keys=None,restrictions=None,plot=0,call=0,save=0,reread=0):
  
  scratch = 0
  if not net:
    scratch = int(raw_input('You want to access data from local (0) or from scratch (1)? '))
    if scratch==1:
      sessionPath = '/scratch01/aschmidt/data/'	#path for storage on scratch
    else:
      sessionPath = 'data/'			#path for local storage
    
    scriptName = pick_script(sessionPath)
    
  repeat = 'y'
  while repeat == 'y':
    if not net:
      pickName = pick_file(scriptName,sessionPath + 'info/',restrictions=restrictions)
      fileName = sessionPath + 'info/' + pickName
    
      net = network('r',fileName)
      net.scriptName = scriptName
    
      #if not (sim_mode == 'eft'):
	#net.scriptName += '_' + sim_mode
      
      sim = simulation_data('r',fileName,net)
      sim.cluster['scratch'] = scratch
    # was hiervon brauche ich wirklich?
    net.data_structure(sim,call=1)
    net.initDataSv(sim)
    
    if sim_mode == 'eft':
      out = analyze_eft(net,sim,plot,save,reread)
    elif sim_mode == 'LE':
      out = analyze_LE(net,sim,plot,save)
    elif sim_mode == 'statistics':
      out = analyze_statistics(net,sim,plot,save)
    elif sim_mode == 'mosaic':
      out = analyze_mosaic(net,sim,plot,save)
    elif sim_mode == 'samples':
      out = analyze_samples(net,sim,plot,save,reread)
    elif sim_mode == 'drive_cur':
      out = analyze_drive_cur(net,sim,plot,save)
    #elif sim_mode == 'pert_spread':
      #out = pert_spread(fileName=net.Path['results_links'],steps=sim.steps,net=net)
            
    if call:
      repeat = 'n'
    else:
      net = None
      sim = None
      repeat = raw_input('You want to analyze more from this simulation? (y/n) ')

  if call:
    return out


def analyze_eft(net,sim,plot,save,reread):
  
  analysis = {}
  
  if 'SONET' in net.scriptName:
    para_str = ''
    title_str = r'$'
    for key in net.topoConst.keys():
      if 'alpha' in key:
	para_str += '_%s=%g' % (key,net.topoConst[key])
	str_tmp = key.split('_')
	title_str += r'\displaystyle \%s_{%s}=%g,\;' % (str_tmp[0],str_tmp[1],net.topoConst[key])
    title_str += r'$'
  else:
    title_paras = ['N','K','J0','rateWnt']
    para_str = 'drive_rate=%g_drive_cplg=%g' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
    title_str = r'$\displaystyle\bar{\nu}_p=%g, J_p=%g$, ' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
      
  # read from already processed data    
  ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r')

  sim.steps['ps'] = ncid.dimensions['ps_steps']
  try:
    net.DataSv['pert_size'] = ncid.variables['pert_size'][:]
  except:
    net.DataSv['pert_size'] = ncid.variables['perturbation_sizes'][:]
    
  net.DataSv['div_track'] = ncid.variables['div_track'][:]
  net.DataSv['p'] = ncid.variables['p'][:]
  net.DataSv['CDBTimes'] = ncid.variables['CDBTimes'][:]
  net.DataSv['distances'] = ncid.variables['distances'][:]
  net.DataSv['correct_ct'] = ncid.variables['correct_ct'][:]
  #net.DataSv['flux_tube_radius'] = ncid.variables['flux_tube_radius'].getValue()
  ncid.close()
  
  sim.steps['ps_ct'] = sim.steps['ps'] - (1 + np.sum(net.DataSv['correct_ct'] == 0))
  
  net.eft_exp()		# mean fluxtube radii (from analysis & from finding)
  if reread:
    net.DataSv['p'] = np.zeros(sim.steps['ps'])
    net.DataSv['div_track'] = (net.DataSv['distances'] > 0.05).astype('int')
    for ps_ct in range(sim.steps['ps']):
      mask = np.invert(np.isnan(net.DataSv['div_track'][ps_ct]))
      if net.DataSv['correct_ct'][ps_ct]:
	net.DataSv['p'][ps_ct] = np.sum(net.DataSv['div_track'][ps_ct][mask])/float(net.DataSv['correct_ct'][ps_ct])
      else:
	net.DataSv['p'][ps_ct] = 0
    #print net.DataSv['p']
    net.fluxSimSv(sim.steps)
    
  div_track = net.DataSv['div_track'].astype('float')
  
  #net.DataSv.pop('flux_tube_radius')
  np.place(div_track,net.DataSv['div_track']==-1,np.nan)	# does this work?
  
  net.DataSv['div_track'] = div_track
  
  #if np.sum(np.isnan(net.DataSv['div_track'])) or not net.DataSv['p'][-1]:
    #os.remove(net.Path['results'] + 'results_processed.nc')
    #net.fluxSimSv(sim.steps)
  
  var_p = np.zeros(sim.steps['ps'])
  p_conv = np.zeros((sim.steps['ps'],sim.steps['total']))
  p_conv[:] = np.nan
  
  for ps_ct in range(sim.steps['ps']):
    mask = np.invert(np.isnan(net.DataSv['div_track'][ps_ct]))
    p_conv[ps_ct][:np.sum(mask)] = np.cumsum(net.DataSv['div_track'][ps_ct][mask])/(np.arange(np.sum(mask))+1).astype('float')
    # maximum - minimum value in last half
    #print p_conv[ps_ct]
    mask = np.invert(np.isnan(p_conv[ps_ct][sim.steps['total']/2:]))
    if len(p_conv[ps_ct][sim.steps['total']/2:][mask]):
      var_p[ps_ct] = np.max(p_conv[ps_ct][sim.steps['total']/2:][mask])-np.min(p_conv[ps_ct][sim.steps['total']/2:][mask])
    else:
      var_p[ps_ct] = np.nan
  p_conv = p_conv.transpose()
  analysis['eft'], analysis['var_eft'], var_p = net.eft_calc(var_p)
  #net.DataSv['flux_tube_radius'] = analysis['eft']
  
  #for ps_ct in range(len(net.DataSv['pert_size'])):
    #print np.min(net.DataSv['distances'][ps_ct])
  
  if plot and analysis['eft']:

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    plt.figure()
    gs = gridspec.GridSpec(3, 2)
  
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[1:3,0])
    ax2 = plt.subplot(gs[1:3,1])
    #fig = plt.figure()
    #ax1 = plt.subplot(121)
    #ax2 = plt.subplot(122)
    #print p_conv.shape
    
    net.DataSv['distances'] = net.DataSv['distances'].reshape(np.prod(net.DataSv['distances'].shape))
    ax0.hist(np.log10(net.DataSv['distances']),bins=100,range=[-8,1])
    #print min(np.log10(net.DataSv['distances']))
    #print max(np.log10(net.DataSv['distances']))
    
    ax0.set_xlim([-8,1])
    #ax0.set_xscale('log')
    
    ax1.plot(np.arange(sim.steps['total']),p_conv)
    ax1.set_xlim([0,sim.steps['total']])
    ax1.set_xlabel('simulation step n')
    ax1.set_ylim([0,1])
    ax1.set_ylabel('Probability of divergence',fontsize=18)
    
    # curves
    appr_plot_range = analysis['eft']*10**np.linspace(-1.5,1.5,1000)
    appr_exp_curve = 1-np.exp(-appr_plot_range/net.eft_exp)
    appr_curve = 1-np.exp(-appr_plot_range/analysis['eft'])
    
    ax2.plot(net.DataSv['pert_size'],net.DataSv['p'],'.',color='b',markersize=10)      
    ax2.plot([net.eft_exp,net.eft_exp],[0,1],color='0.75')
    ax2.plot([analysis['eft'],analysis['eft']],[0,1],color='k')
    ax2.plot(appr_plot_range,appr_exp_curve,'--',color='0.75')
    ax2.plot(appr_plot_range,appr_curve,'--',color='k')
    ax2.errorbar(net.DataSv['pert_size'],net.DataSv['p'],yerr=var_p,ecolor='r',elinewidth=2,fmt=None)
    ax2.errorbar([analysis['eft'],analysis['eft']],[1-1/np.exp(1),1-1/np.exp(1)],xerr=analysis['var_eft'],ecolor='r',fmt=None)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    #ax2.set_xlim([net.DataSv['pert_size'].min(),net.DataSv['pert_size'].max()])
    ax2.set_xlim([10**(-4),10])
    ax2.set_xlabel(r'$\displaystyle \varepsilon_{FT}$',fontsize=18)
    ax2.set_ylim([10**(-2),10**(0)])
    
    plt.tick_params(labelsize=14)	# adjust to get proper labelsizes
    plt.title(title_str)
    plt.show(block=False)
  try:
    print 'found flux tube radius: %5.3f+/-%5.3f' %(analysis['eft'],analysis['var_eft'])
  except:
    print "No fluxtube radius found."
    print analysis

  return analysis


def analyze_LE(net,sim,plot,save):
  
  analysis = {}
  
  # read data
  for dat in os.listdir(net.Path['results']):
    Data = readDataOut(net.Path['results']+dat)
  
  # 
  if 'subLyapunovExponents' in Data.keys():
    analysis['LyapunovExponents'] = Data['subLyapunovExponents']
    Data['LEconvergence'] = Data['subLEconvergence']
  else:
    analysis['LyapunovExponents'] = Data['LyapunovExponents']
  
  if plot:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    N = len(analysis['LyapunovExponents'])
    
    # plot CLV, calculate "Chaos index"
    fig = plt.figure()

    #ax1 = plt.subplot(121)
    #axcb = plt.subplot(1,30,16)
    ax2 = plt.subplot(121)
    ax3 = plt.subplot(122)
    
    #levs = range(80)
    #assert len(levs) % 2 == 0, 'N levels must be even.'
    #zero_val = abs(Data['localLE'].min())/(Data['localLE'].max()-Data['localLE'].min())
    #rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors=[(0,(0,0,1)),(zero_val,(1,1.,1)),(1,(1,0,0))],N=len(levs)-1,)
    #dim0 = Data['localLE'].shape[0]
    #dim1 = N#Data['localLE'].shape[1]
    #times = np.tile(Data['LEtimes'][0:dim0,np.newaxis], (1,dim1))	#times = np.tile(Data['LEtimes'][0:dim0], (dim1,1))
    #index = np.tile(np.linspace(1/float(dim1),1,dim1), (dim0,1))	#index = np.tile(np.linspace(1/float(dim1),1,dim1)[:,np.newaxis], (1,dim0))
    
    #if Data['localLE'].shape[1] > N:
      #N_idx = N
    #else:
      #N_idx = 0
    #Data['localLE'] = Data['localLE'].transpose()[N_idx:].transpose()
    
    str_add = ''
    for item in net.topoConst.items():
      if 'alpha' in item[0] and net.topo == 'S':
	str_add += '_%s=%g' % (item[0],item[1])
      elif not ('alpha' in item[0]):
	if not (item[0] == 'drive_type'):
	  str_add += '_%s=%g' % (item[0],item[1])
    
    ### get covariance of localLE with indegree/outdegre and/or firing rate
    #pc = ax1.pcolor(times, index, Data['localLE'],cmap=rwb)
    #ax1.set_xlabel('time (ms)')
    #ax1.set_ylabel('i / N')
    #ax1.set_xlim([0,np.max(Data['LEtimes'][0:dim0])])
    #ax1.set_title('local Lyapunov exponents per spike')
    #plt.colorbar(pc, cax=axcb)
    
    ax2.plot(Data['LEtimes'],Data['LEconvergence'].transpose()[1])
    ax2.plot(Data['LEtimes'],Data['LEconvergence'].transpose()[N/2])
    ax2.plot(Data['LEtimes'],Data['LEconvergence'].transpose()[-1])
    ax2.set_xlabel('time')
    ax2.set_ylabel('Lyapunov Exponent')
    ax2.set_xlim([0, max(Data['LEtimes'])])
    ax2.set_title('LE convergence')
    ax2.set_ylim([-120,0])
    ax3.plot(np.linspace(1/float(N),1,N), analysis['LyapunovExponents'])
    #ax3.plot(np.linspace(1/float(N),1,N), Data['LEclv'][N_idx:])
    ax3.set_ylabel(r'$\displaystyle \lambda_i (s^{-1})$')
    ax3.set_xlabel('i/N')
    ax3.set_ylim([-120,0])
    
    #ax3.set_ylim([np.min(analysis['LyapunovExponents']),np.max(analysis['LyapunovExponents'])])
    plt.suptitle(r'Lyapunov spectra for N=%d, K=%d, $\displaystyle\nu=%g %s$'%(net.Const['N'],net.Const['K'],net.Const['rateWnt'],',\; '.join(str_add.split('_'))),fontsize=16)
    
    if save:
      save_str = './pics/%s/LE/single_N=%d_K=%d_f=%3.1f%s.png' % (net.scriptName.split('_')[0],net.Const['N'],net.Const['K'],net.Const['rateWnt'],str_add)
      plt.savefig(save_str)
      print "Figure saved to %s." % save_str
    else:
      plt.show(block=False)
  
  return analysis


def analyze_statistics(net,sim,plot,save):
  stats = ['ISIdecorr','moments','subthreshold','topology']		  # stats = ['topology','moments','ISIdecorr','subthreshold'], choose which ones to plot
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
    
  analysis = {}
  
  N = net.Const['N']
  
  str_add = []
  for item in net.topoConst.items():
    if 'alpha' in item[0] and not (net.topo == 'p'):
      str_add.append('%s=%g' % (item[0],item[1]))
    elif not ('alpha' in item[0]):
      str_add.append('%s=%g' % (item[0],item[1]))
      
  # specify values that should be read out
  read_list = ['rateNeurons', 'cvNeurons', 'mean_ISIdecorr', 'cv_ISIdecorr', 'chi','finalCurrents','measure1stStateVarNeurons','measureTimes']
  read_list_topo = ['p_hat', 'alpha_recip_hat', 'alpha_conv_hat', 'alpha_div_hat', 'alpha_chain_hat', 'inDegree', 'outDegree']
  
  Data = {}
  numsim = len(os.listdir(net.Path['results']))
  
  analysis['inDegree'] = np.zeros((numsim,N))
  analysis['outDegree'] = np.zeros((numsim,N))
  
  idx_file = 0
  for dat in os.listdir(net.Path['results']):
    Data, ncidTopo = readDataStats(net.Path['results'] + dat,Data,read_list,numsim,idx_file)	# read results file
    Data = readDataStats(ncidTopo,Data,read_list_topo,numsim,idx_file)				# read topology file

    analysis['inDegree'][idx_file] = Data['inDegree'][idx_file][:N]		# truncate to size of recurrent population
    analysis['outDegree'][idx_file] = Data['outDegree'][idx_file][:N]
    
    # resize measures to proper shape and population
    if not idx_file:
      measureSz = Data['measureTimes'].shape[-1]
      tmp_measureSz = Data['measure1stStateVarNeurons'].shape[-1]
      
      analysis['measures'] = np.zeros((numsim,N,measureSz))
      analysis['measureTimes'] = Data['measureTimes'][idx_file]

    analysis['measures'][idx_file] = Data['measure1stStateVarNeurons'][idx_file].reshape(measureSz,tmp_measureSz/measureSz).transpose()[:N]
    
    idx_file += 1
  
  np.place(Data['rateNeurons'],np.isnan(Data['rateNeurons']),0)	# nan-entries == silent neurons
  
  # compute statistics...
  print Data['finalCurrents']
  # ... topology ...
  if 'topology' in stats:
    Data['rateNetwork'] = np.zeros(numsim)
    Data['inDegreeNet'] = np.zeros(numsim)
    Data['outDegreeNet'] = np.zeros(numsim)
    Data['corr_KNu'] = np.zeros(numsim)
    Data['corr_KK'] = np.zeros(numsim)
    Data['rateNetwork'] = np.zeros(numsim)
    Data['rateNetwork'] = np.zeros(numsim)
    
    for idx_file in range(numsim):
      Data['rateNetwork'][idx_file] = np.mean(Data['rateNeurons'][idx_file])
      Data['inDegreeNet'][idx_file] = np.mean(analysis['inDegree'][idx_file])
      Data['outDegreeNet'][idx_file] = np.mean(analysis['outDegree'][idx_file])
      
      # correlation of in-degree and firing rate
      Data['corr_KNu'][idx_file] = np.mean((Data['rateNeurons'][idx_file]-Data['rateNetwork'][idx_file])*(analysis['inDegree'][idx_file]-Data['inDegreeNet'][idx_file]))/np.sqrt(np.var(Data['rateNeurons'][idx_file])*np.var(analysis['inDegree'][idx_file]))

      # correlation of in- and out-degree
      Data['corr_KK'][idx_file] = np.mean((analysis['outDegree'][idx_file]-Data['outDegreeNet'][idx_file])*(analysis['inDegree'][idx_file]-Data['inDegreeNet'][idx_file])) /np.sqrt(np.var(analysis['outDegree'][idx_file])*np.var(analysis['inDegree'][idx_file]))
    
    analysis['corr_KNu'] = np.mean(Data['corr_KNu'])
    analysis['corr_KK'] = np.mean(Data['corr_KK'])
    
    analysis['inDegree'] = np.reshape(analysis['inDegree'],N*numsim)
    analysis['outDegree'] = np.reshape(analysis['outDegree'],N*numsim)
    
    if plot:
      K = net.Const['K']
      y_border = 0.1
      ax_border = [2*K,2*K]
      if max(analysis['inDegree']) > 2*K:
	ax_border[0] = N
      if max(analysis['outDegree']) > 2*K:
	ax_border[1] = N
      
      fig, axes = plt.subplots(nrows=1,ncols=2)
      
      axes[0].hist(analysis['inDegree'],bins=min(ax_border[0],100),range=[0,ax_border[0]],color='k')
      axes[1].hist(analysis['outDegree'],bins=min(ax_border[1],100),range=[0,ax_border[1]],color='k')
      
      axes[0].set_xlabel('In-degree',fontsize=18)
      axes[1].set_xlabel('Out-degree',fontsize=18)
      axes[0].set_ylabel('Probability',fontsize=18)
      axes[0].set_xlim([0,ax_border[0]])
      axes[1].set_xlim([0,ax_border[1]])
      axes[0].set_ylim([0,y_border*N*numsim])
      axes[1].set_ylim([0,y_border*N*numsim])
      
      axes[0].set_xticks(np.linspace(0,ax_border[0],5))
      axes[0].set_xticklabels(np.linspace(0,ax_border[0],5).astype('int'),fontsize=14)
      axes[1].set_xticks(np.linspace(0,ax_border[1],5))
      axes[1].set_xticklabels(np.linspace(0,ax_border[1],5).astype('int'),fontsize=14)
      
      axes[0].set_yticks(np.linspace(0,y_border*N*numsim,3))
      axes[0].set_yticklabels(np.linspace(0,y_border,3),fontsize=14)
      axes[1].set_yticks(np.linspace(0,y_border*N*numsim,3))
      axes[1].set_yticklabels(np.linspace(0,y_border,3),fontsize=14)
      
      #axes[0].legend(prop={'size':14})
      #axes[1].legend(prop={'size':14})
      #plt.suptitle(r'Erd\H{o}s-R\'{e}nyi network',fontsize=20)
      if save:
	save_str = './pics/%s/statistics/degrees_%s.png' % (net.scriptName.split('_')[0],'_'.join(str_add))
	plt.savefig(save_str)
	print 'Figure saved to ' + save_str
      plt.show(block=False)
  
  # ... synchrony ...
  analysis['chi'] = np.mean(Data['chi'])
  
  
  # ... moments ...
  if 'moments' in stats:
    m = 0
    search_range = 4
    num_bin = 21

    moment_bars = [50,1.5]
    moment_labels = ['Firing rate','CV']
    y_bars = [0.2,0.1]
    
    analysis['peaks'] = np.zeros((2,2,4))	#store up to 4 peaks (moments,axes,peaks)
    analysis['peaks'][:] = np.nan      
    moment = np.zeros((2,2,num_bin))	#store moments (moments,axes,bins)
	
    for key_moment in ['rateNeurons','cvNeurons']:	#iterate over moments
      analysis[key_moment] = np.reshape(Data[key_moment],np.prod(Data[key_moment].shape))  # flatten array
      
      # convert to histograms
      tmp_hist = np.histogram(analysis[key_moment],bins=num_bin,range=[0,moment_bars[m]])
      moment[m][0] += tmp_hist[0]/float(numsim*N)
      moment[m][1] = tmp_hist[1][:-1]
      
      # search for peaks in histograms
      peak_idx = 0
      
      for idx_dat in range(len(moment[m][0])):
	peak = 1
	moment_search = moment[m][0][max(0,idx_dat-search_range):min(num_bin,idx_dat+search_range+1)]	
	
	if np.sum(moment[m][0][idx_dat] <= moment_search)==1:
	  analysis['peaks'][m][0][peak_idx] = moment[m][1][idx_dat]
	  analysis['peaks'][m][1][peak_idx] = moment[m][0][idx_dat]
	  peak_idx += 1
      
      if plot:	# convert data to histograms and find peaks
      
	plt.figure()
	gs = gridspec.GridSpec(3, 1)

	ax0 = plt.subplot(gs[0,0])
	ax1 = plt.subplot(gs[1:3,0])

	mask = np.invert(np.isnan(analysis[key_moment]))
	
	ax0.boxplot(analysis[key_moment][mask],notch=1,vert=False,widths=0.5)
	ax0.plot(analysis['peaks'][m][0],np.ones(4),'o',color='grey',markersize=10,label=moment_labels[m] + 'peak')
	ax0.set_xlim([0,moment_bars[m]])
	ax0.set_xticklabels([])
	ax0.set_yticklabels([])
	
	ax1.plot(moment[m][1],moment[m][0])
	ax1.plot(analysis['peaks'][m][0],analysis['peaks'][m][1],'o',color='grey',markersize=10,label='peak')
	ax1.set_xlim([0,moment_bars[m]])
	ax1.set_ylim([0,0.2])
	ax1.set_ylabel('Probability',fontsize=18)
	ax1.set_xlabel(moment_labels[m],fontsize=18)
	
	ax1.legend(prop={'size':14})
	#plt.suptitle(r'Erd\H{o}s-R\'{e}nyi network',fontsize=20)
	if save:
	  save_str = './pics/%s/statistics/moments%d_%s.png' % (net.scriptName.split('_')[0],m,'_'.join(str_add))
	  plt.savefig(save_str)
	  print 'Figure saved to ' + save_str

	plt.show(block=False)
      
      m += 1
  
  
  # ISI stats of potential decorrelation events
  if 'ISIdecorr' in stats:
    # put them all together
    analysis['mean_ISIdecorr'] = np.reshape(Data['mean_ISIdecorr'].transpose(2,0,1),(2,N*numsim))
    analysis['cv_ISIdecorr'] = np.reshape(Data['cv_ISIdecorr'].transpose(2,0,1),(2,N*numsim))
    analysis['finalCurrents'] = np.mean(Data['finalCurrents'],axis=0)[0]
    
    if plot:
      mask = np.invert(np.isnan(analysis['mean_ISIdecorr']))
      N_eff = np.sum(mask[1])
      
      #x_max = np.ceil(np.max(Data['ISIdecorr'][mask])/10)*10
      x_max = 20
      
      plt.figure()
      
      ax0 = plt.subplot(211)
      
      ax0.hist(analysis['mean_ISIdecorr'][0],range=[0,x_max],bins=41,alpha=0.7,label='post-prespike time')
      ax0.hist(analysis['mean_ISIdecorr'][1],range=[0,x_max],bins=41,alpha=0.7,label='pre-postspike time')
      
      ax0.set_xlim([0,x_max])
      ax0.set_xlabel(r'$\displaystyle \Delta t$ in ms')
      ax0.set_ylim([0,N_eff])
      ax0.set_yticks(np.linspace(0,N_eff,3))
      ax0.set_yticklabels(np.linspace(0,1,3),fontsize=14)
      ax0.legend()
      ax0.set_xlabel(r'$\displaystyle ISI_{decorr}$',fontsize=18)
      
      ax1 = plt.subplot(212)
      ax1.hist(analysis['cv_ISIdecorr'][0],range=[0.8,1.8],bins=21,alpha=0.7,label='cv post-prespike time')
      ax1.hist(analysis['cv_ISIdecorr'][1],range=[0.8,1.8],bins=21,alpha=0.7,label='cv pre-postspike time')

      ax1.set_xlim([0.8,1.8])
      #ax1.set_xlabel('CV')
      
      ax1.set_xlabel(r'CV of $\displaystyle ISI_{decorr}$',fontsize=18)
      ax1.set_ylim([0,N_eff])
      ax1.set_yticks(np.linspace(0,N_eff,3))
      ax1.set_yticklabels(np.linspace(0,1,3),fontsize=14)
      ax1.legend()
      
      plt.suptitle(r'Erd\H{o}s-R\'{e}nyi network',fontsize=20)
      if save:
	save_str = './pics/%s/statistics/ISIdecorr_%s.png' % (net.scriptName.split('_')[0],'_'.join(str_add))
	plt.savefig(save_str)
	print 'Figure saved to ' + save_str
      #plt.suptitle(r'ISI of potential decorr. events (%d of %d neurons active)' % (N_eff/numsim,N),fontsize=20)
      plt.show(block=False)
  
  
  # subthreshold statistics: phase distributions
  
  #single trajectories (randomly drawn from first simulation)
  if 'subthreshold' in stats:
    if plot:
      plt.figure()
      plt.plot(analysis['measureTimes'][:1000],analysis['measures'][0][np.random.randint(N)][:1000])
      plt.ylim([-1,1])
      plt.show(block=False)
      
      # distribution of phases in network (assuming ergodicity, thus distribution is the same at all times)
      phase_distr, phase_bins = np.histogram(analysis['measures'].reshape(np.prod(analysis['measures'].shape)),bins=100,range=[-1,1])
      plt.figure()
      plt.hist(analysis['measures'].reshape(np.prod(analysis['measures'].shape)),bins=100,range=[-1,1])
      print "Neuron fraction prior to spiking: %4.3f" % (phase_distr[-1]/float(np.sum(phase_distr))*100)
      print "phase mean: %g" % np.mean(analysis['measures'])
      print "phase CV: %g" % np.var(analysis['measures'])
      plt.show(block=False)
  
  
  return analysis
    

def analyze_mosaic(net,sim,plot,save):
  
  ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r')
  n = ncid.dimensions['n']
  pert_size = ncid.variables['pert_size'].getValue()
  fluxtube = ncid.variables['fluxtubes'][:]
  ncid.close()
  
  print pert_size**2
  # number of different fluxtubes
  ct_fluxtube = np.max(fluxtube)
  
  area = np.zeros(ct_fluxtube+1)
  # area of fluxtubes
  for idx_fluxtube in range(ct_fluxtube+1):
    area[idx_fluxtube] = np.sum(fluxtube==idx_fluxtube)/float(n**2)*pert_size**2
  
  #print fluxtube
  #print ct_fluxtube
  print "area: %g" % np.mean(area)
  print "area CV: %g" % (np.sqrt(np.var(area)/np.mean(area)**2))
  #counter = np.bincount(fluxtube)
  #print counter
  
  plt.matshow(fluxtube)
  plt.xticks([])
  plt.yticks([])
  plt.show()
  
  return fluxtube


def analyze_samples(net,sim,plot,save,reread):
  
  net.analysis = {}
  if not os.path.exists(net.Path['results'] + 'results_processed.nc'):
    print "read it"
    net.prob_read(sim)		# get net.analysis
  if reread:
    ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r')
    #print ncid.variables.keys()
    # hand over all variables
    for key in ncid.variables.keys():
      net.analysis[key] = np.copy(ncid.variables[key][:])
    if reread == 2:
      time_len = len(net.analysis['measureTimes'])
      for idx_ct in range(net.analysis['distanceOrth'].shape[0]):
	if net.analysis['distanceOrth'][idx_ct][-1] < 0.1:
	  net.get_convergence(time_len,idx_ct)
	else:
	  try:
	    net.analysis['convergeVelocity'][0][idx_ct] = np.nan
	    net.analysis['convergeVelocity'][1][idx_ct] = np.nan
	    net.analysis['kinkTime'][idx_ct] = np.nan
	  except:
	    net.analysis['convergeVelocity'] = np.zeros((2,net.analysis['distanceOrth'].shape[0]))
	    net.analysis['convergeVelocity'][:] = np.nan
	    
      net.writeSamples(sim,time_len,sim.steps['ps'])
  #print net.analysis.keys()
  ## implement reading from topo-file here to get connections of spike-crossing neurons
  
  if plot:
    # get title and parastring from net.Const/net.topoConst
    if 'SONET' in net.scriptName:
      para_str = ''
      title_str = r'$'
      for key in net.topoConst.keys():
	if 'alpha' in key:
	  para_str += '_%s=%g' % (key,net.topoConst[key])
	  str_tmp = key.split('_')
	  title_str += r'\displaystyle \%s_{%s}=%g,\;' % (str_tmp[0],str_tmp[1],net.topoConst[key])
      title_str += r'$'
    else:
      title_paras = ['N','K','J0','rateWnt']
      para_str = 'drive_rate=%g_drive_cplg=%g' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
      title_str = r'$\displaystyle\bar{\nu}_p=%g, J_p=%g$, ' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
      
      for para in title_paras:
	title_str += '$%s=%g$, ' %(para,net.Const[para])
    print net.analysis.keys()
    sim_idx = np.random.randint(sim.steps['co']*sim.steps['ps'])
    
    spiked = np.zeros(net.analysis['trainNeuronDrive'].shape[-1])
    for t in range(net.analysis['trainNeuronDrive'].shape[-1]):
      spiked[t] = len(np.unique(net.analysis['trainNeuronDrive'][sim_idx][:t]))/float(net.Const['N'])
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    
    fig,axes = plt.subplots(nrows=2,ncols=1)
    #print sim.steps['total']

    axes[0].plot(net.analysis['trainTimeDrive'][sim_idx],spiked,'-k')
    axes[0].set_xlim([0,net.analysis['measureTimes'].max()])
    axes[0].set_ylim([0,1])
    
    axes[1].plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][sim_idx])
    axes[1].set_yscale('log')
    axes[1].set_ylim([10**-12,10**1])
    
    plt.plot()
    
    
    plt.figure()
    plt.show(block=False)
    if plot == 2:	# interactive mode on if video/stream
      plt.ion()
    
    gs = gridspec.GridSpec(2, 2)
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    
    # try estimated decay
    
    mask = np.invert(np.isnan(net.analysis['convergeVelocity'][0]))
    mean_convVel1 = np.mean(net.analysis['convergeVelocity'][0][mask])
    
    mask = np.invert(np.isnan(net.analysis['convergeVelocity'][1]))
    mean_convVel2 = np.mean(net.analysis['convergeVelocity'][1][mask])
    
    decay1 = np.exp(net.analysis['measureTimes']*mean_convVel1)
    decay2 = np.exp(net.analysis['measureTimes']*mean_convVel2)
    
    ax0.plot(net.analysis['measureTimes'],decay1,'--k',linewidth=5,label='approximate initial convergence')
    ax0.plot(net.analysis['measureTimes'],decay2,'--r',linewidth=5,label='approximate later convergence')
    
    #ax0.plot(net.analysis['measureTimes'],np.abs(net.analysis['distance1']).transpose(),'-r')	# orthogonalized distance
    
    ax0.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'].transpose(),'-b')	# orthogonalized distance
    ax0.set_xlim([0,2])
    ax0.set_ylim([10**(-16),6])
    ax0.set_yscale('log')
    ax0.set_xlabel('time in s')
    ax0.set_ylabel('distance from reference')
    ax0.legend(loc=4)
    
    print net.analysis['kinkTime']
    ax1.plot(range(len(net.analysis['kinkTime'])),net.analysis['kinkTime'],'ok')

    #print np.mean(net.analysis['convergeVelocity1'])
    #print np.mean(net.analysis['convergeVelocity2'])
    ax2.plot(range(len(net.analysis['convergeVelocity'][0])),net.analysis['convergeVelocity'][0],'ko')
    ax2.set_ylabel('initial convergence speed')
    ax2.set_xlabel('trial \#')
    
    ax3.plot(range(len(net.analysis['convergeVelocity'][1])),net.analysis['convergeVelocity'][1],'ro')
    ax3.set_ylabel('later convergence speed')
    ax3.set_xlabel('trial \#')
    
    plt.suptitle(title_str)
    plt.show(block=False)
    #for idx_ct in range(sim.steps['total']):#range(sim.steps['pd']):
      #png_name = '%s_%s_%d' % (net.scriptName,para_str,idx_ct)
      ##print png_name

      ##ax0.cla()
      ## the 2-norm distances:
      #ax0.plot(net.analysis['measureTimes'],net.analysis['distance'][idx_ct],'-k')	# distance
      #ax0.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][idx_ct],'-r')	# orthogonalized distance
      ## the 1-norm:
      #ax0.plot(net.analysis['measureTimes'],np.abs(np.sum(net.analysis['distanceVector'][idx_ct],axis=1)),'k')	# distance
      #ax0.plot(net.analysis['measureTimes'],np.abs(np.sum(net.analysis['distanceOrthVector'][idx_ct],axis=1)),'r')	# orthogonalized distance
      
      #ax0.set_xlim([0,net.analysis['measureTimes'][-1]])
      #ax0.set_ylim([10**(-20),6])
      ##ax0.set_ylim([10**(-6),10])
      
      #ax0.set_yscale('log')
      #ax0.set_ylabel('distance to reference')
      #ax0.set_xlabel('time t in s')
      ##str_title = net.scriptName.split('_')[0]
      #ax0.set_title(title_str)
      
      #plt.draw()
      
      #if (net.analysis['distance'][idx_ct][-1] < 0.1) and (plot == 2):
	#t = 0
	#while np.linalg.norm(net.analysis['distanceVector'][idx_ct][t]) > 10**(-10):

	  #ax2.cla()
	  #ax2.plot(range(net.Const['N']),net.analysis['distanceOrthVector'][idx_ct][t],'k',label='orthogonalized distance')
	  #ax2.set_xlabel('Neuron Index')
	  
	  #ax1.cla()
	  #ax1.plot(range(net.Const['N']),dist_tmp[t],'k',label='total distance')
	  #if np.max(dist_tmp[t] > 10**(-5)):
	    #border = 10**(-3)
	  #elif np.max(dist_tmp[t] > 10**(-7)):
	    #border = 10**(-5)
	  #else:
	    #border = 10**(-7)
	  
	  #ax1.set_ylim([-border,border])
	  #ax2.set_ylim([-border,border])

	  #ax1.set_yticks(np.linspace(-border,border,5))
	  #ax1.set_yticklabels(np.linspace(-border,border,5))

	  #ax1.set_xlabel('Neuron Index')
	  
	  #ax0.plot(net.analysis['measureTimes'][t],10**(-16),'Dk') # plot process of video
	  
	  #plt.draw()
	  
	  #if save == 2:
	    #plt.savefig('./pics/vid/%s_%d.png'%(png_name,t))
	  #t += 1
	  #if t >= len(net.analysis['measureTimes']) or (net.analysis['distanceOrth'][idx_ct][t-1] < 10**(-8)):
	    #break
      #if save == 2:
	## convert to video
	#video_string = "ffmpeg -f image2 -r 25 -i './pics/vid/" + png_name + "_%d.png' ./pics/vid/" + png_name + ".mp4"
	#print 'Saving video as "%s.mp4"' % png_name
	#runIt(video_string)
      
	#for dat in os.listdir('./pics/vid/'):	#clean up
	  #if ('%s_%s_%d_' % (net.scriptName,para_str,pd_ct)) in dat:
	    #os.remove('./pics/vid/' + dat)    
    

    
    #Data = []
    #for dat in os.listdir(net.Path['results']):
      #if 'DataOut' in dat:
	#Data.append(readDataOut(net.Path['results'] + dat))
    
    ## plot spike trains
    #plt.figure()
    #idx = 100
    #for i in range(sim.steps['pd']):
      ##dist = np.sqrt(np.sum(
      #train = Data[i]['trainTime'][np.where(Data[i]['trainNeuron'] == 100)]
      
      #plt.plot(train,np.ones(len(train))*i,'ok')
      #plt.plot(Data[i]['trainTime'][:idx],np.ones(idx)*i,'ok')
      ##print Data[i]['trainTime'][idx]
    #plt.ylim([-1,10])
    #plt.show(block=False)
  return net.analysis


def erfi(z):
  #print z<0
  return (complex(0.,-1.)*spec.erf(complex(0.,1.)*z)).real


def analyze_drive_cur(net,sim,plot,save):
  
  K = net.Const['K']
  rateWnt = net.Const['rateWnt']
  tauM = net.Const['tauM']/1000.
  V_T = 1
  V_R = 1
  
  #print net.Const
  analyze = {}
  
  
  I_sim_tmp = np.zeros(sim.steps['co'])
  i=0
  for dat in os.listdir(net.Path['results']):
    ncid = netcdf.netcdf_file(net.Path['results'] + dat,'r')
    I_sim_tmp[i] = ncid.variables['finalCurrents'][:][0]
    ncid.close()
    i+=1
    
  mask = np.where(I_sim_tmp)
  print I_sim_tmp
  analyze['I_sim'] = np.mean(I_sim_tmp[mask])
  
  
  # driving current from balance equation
  analyze['I_bal'] = -rateWnt*net.Const['J0']*tauM
  
  # driving current from selfconsistency equation (Brunel,Hakim)
  V_steps = 1001
  V = np.linspace(-0.5,1,V_steps)
  
  selfcon_sv = 'data/selfcon_K%d_rateWnt%g.nc' % (K,rateWnt)
  runStr = './Code/other/calc_I %d %g %s' % (K,rateWnt,selfcon_sv)
  runIt(runStr)
  #try:
  ncid = netcdf.netcdf_file(selfcon_sv,'r')
  analyze['I_selfcon'] = ncid.variables['I0'].getValue()
  ncid.close()
  
  os.remove(selfcon_sv)
  
  if plot:
    sigma2 = net.Const['J0']**2*rateWnt*tauM
    sigma = np.sqrt(sigma2)
    mu = V_T + np.sqrt(K)*(analyze['I_selfcon'] + net.Const['J0']*rateWnt*tauM)
    print mu
    
    erfi_tmp = np.zeros(V_steps)
    
    for j in range(V_steps):
      if V[j] < 0:
	erfi_tmp[j] = erfi((V_R-mu)/sigma)
      else:
	erfi_tmp[j] = erfi((V[j]-mu)/sigma)

    P_V = np.sqrt(math.pi)*rateWnt*tauM/sigma*np.exp(-(V-mu)**2/sigma2)*(erfi((V_T-mu)/sigma)-erfi_tmp)
    
    plt.figure()
    plt.plot(V,P_V,'-k')
    plt.show()
      
  #except:
    #print "Parameters K=%d and rateWnt=%g gave no solution!" % (K,rateWnt)
    #analyze['I_selfcon'] = np.nan
  
  return analyze

  
def assignSpikes(Data,idx_ct,idx_cut):
  
  assigned = []
  for i in range(idx_cut):
    
    spikeRef = Data['trainNeuron'][0][idx_ct][i]
    spike = Data['trainNeuron'][1][idx_ct][i]
    
    spikeTimeRef = Data['trainTime'][0][idx_ct][i]
    spikeTime = Data['trainTime'][1][idx_ct][i]

    if spikeRef==spike and not (i in assigned):
      Data['trainDeltaTime'][idx_ct][i] = abs(spikeTime-spikeTimeRef)
      assigned.append(i)
    else:
      #if (Topo[DataRef_spike][Data_spike] or Topo[Data_spike][DataRef_spike]):
	#DataRef['decorr'].append(DataRef['trainTime'][i])	#get decorrelation event times of postsyn. neurons
    
      # search for the belonging spike
      test_idx_list = np.where(Data['trainNeuron'][1][idx_ct] == spikeRef)[0]	# get pos. of neurons AP in perturbed trajectory
      test_Data_spiketimes = Data['trainTime'][1][idx_ct][test_idx_list]
      test_Data_spiketimes_idx = abs(spikeTimeRef-test_Data_spiketimes) < 10**(-1)
      
      search_list = list(test_Data_spiketimes[test_Data_spiketimes_idx])
      
      found = 0
      while found == 0:
	if not len(search_list):
	  Data['trainDeltaTime'][idx_ct][i] = 1
	  break
	
	test_time = min(search_list)
	test_idx = test_idx_list[np.where(test_time==test_Data_spiketimes)[0]]
	
	if test_idx in assigned:
	  search_list.remove(test_time)
	else:
	  Data['trainDeltaTime'][idx_ct][i] = test_time
	  assigned.append(test_idx)
	  found = 1
  #return deltaTime


def readDataStats(fileName,Data,read_list,numsim,idx_file):
  
  ncid = netcdf.netcdf_file(fileName,'r')
  for key in read_list:
    if not (key in Data.keys()):
      try:
	Data[key] = np.zeros(sum(((numsim,),ncid.variables[key][:].shape),()))
	assert Data[key].shape[1]>1
      except:
	Data[key] = np.zeros(numsim)
    try:
      Data[key][idx_file] = ncid.variables[key][:]
    except:
      Data[key][idx_file] = ncid.variables[key].getValue()
      
  try:
    ncidTopo = ''.join(ncid.variables['inputfile_topology'][:])
    ncid.close()
    return Data, ncidTopo
  except:
    ncid.close()
    return Data

  
def pert_spread(fileName,steps,net=None):
  
  numread = steps['co']*steps['in']*steps['pd']
  
  if os.path.exists(net.Path['results'] + 'results_processed.nc') and numread > 1:
    ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r')
    
    time_div = ncid.variables['time_div'][:]
    vel_div = ncid.variables['velocity_div'][:]
    
    time_mean = ncid.variables['time_mean'].getValue()
    time_var = ncid.variables['time_var'].getValue()
    
    vel_mean = ncid.variables['velocity_mean'].getValue()
    vel_var = ncid.variables['velocity_var'].getValue()
    
    print 'Already read out'
  else:
    print 'Analyze spread of perturbation...'
    
    #search for indices at wich t >= t_ana
    t_ana = 1
    sort_it = 0
    ncid = netcdf.netcdf_file(fileName,'r')  
    
    N = net.Const['N']
    
    idx_decorr = np.zeros(numread).astype('int')
    time_div = np.zeros(numread)
    vel_div = np.zeros(numread)
    
    idx_sim = 0
    
    print 'reading %d files from list %s...' % (numread,fileName)
    
    for co_ct in range(steps['co']):
      for in_ct in range(steps['in']):
	link_reference = ''.join(ncid.variables['links'][0,co_ct,in_ct,0,4])
	#try:
	DataRef = readDataOut(link_reference,1)
	for pd_ct in range(steps['pd']):
	  link_perturbed = ''.join(ncid.variables['links'][0,co_ct,in_ct,pd_ct,3])
	  link_topology = ''.join(ncid.variables['links'][0,co_ct,in_ct,pd_ct,2])
	  Data = readDataOut(link_perturbed,1)
	  ncidTop = netcdf.netcdf_file(link_topology,'r')
	  Topo = postsynapticvector2adjacencymatrix(ncidTop.variables['postSyn'][:],ncidTop.variables['outDegree'][:])
	  ncidTop.close()
	  
	  if net.mode['topo'] == 'p':
	    DataRef['trainTime'] = DataRef['trainTime'][np.where(DataRef['trainNeuron']<N)]
	    Data['trainTime'] = Data['trainTime'][np.where(Data['trainNeuron']<N)]
	    
	    DataRef['trainNeuron'] = DataRef['trainNeuron'][np.where(DataRef['trainNeuron']<N)]
	    Data['trainNeuron'] = Data['trainNeuron'][np.where(Data['trainNeuron']<N)]
	    
	  #print len(Data['trainTime'])
	  
	  if t_ana > DataRef['measureTimes'][-1]:
	    t_ana = DataRef['measureTimes'][-1]
	  try:
	    idx_phase = min(np.where(DataRef['measureTimes']>=t_ana)[0][0],np.where(Data['measureTimes']>=t_ana)[0][0])
	    idx_train = min(np.where(DataRef['trainTime']>=t_ana)[0][0],np.where(Data['trainTime']>=t_ana)[0][0])
	  except:
	    idx_phase = min(DataRef['measureTimes'].shape[0],Data['measureTimes'].shape[0])
	    idx_train = min(DataRef['trainTime'].shape[0],Data['trainTime'].shape[0])

	  # calculate differences in phases
	  DataRef['deltaPhase'] = abs(DataRef['measurePhases'][:idx_phase] - Data['measurePhases'][:idx_phase])
	  DataRef['spread'] = np.sum(DataRef['deltaPhase'] >= 10**(-2),axis=1)/float(N)	# #-neurons with a certain perturbation
	  
	  DataRef['distance'] = np.zeros(idx_phase)
	  DataRef['distance'] = np.sqrt(np.sum(DataRef['deltaPhase']**2,axis=1))
	  
	  ##if DataRef['distance'][-1] > 6:
	  
	  DataRef['deltaPhase'] = DataRef['deltaPhase'].transpose()

	  
	  ### here comes the spike-time deviation ###
	  DataRef['deltaTime'] = np.zeros(idx_train)
	  # search for corresponding spike of each neuron in train
	  assigned = []
	  DataRef['decorr'] = []
	  for i in range(idx_train):
	    DataRef_spike = DataRef['trainNeuron'][i]
	    Data_spike = Data['trainNeuron'][i]
	    
	    DataRef_spiketime = DataRef['trainTime'][i]
	    
	    if DataRef_spike==Data_spike:
	      DataRef['deltaTime'][i] = abs(DataRef['trainTime'][i]-Data['trainTime'][i])
	      assigned.append(i)
	    else:
	      if (Topo[DataRef_spike][Data_spike] or Topo[Data_spike][DataRef_spike]):
		DataRef['decorr'].append(DataRef['trainTime'][i])	#get decorrelation event times of postsyn. neurons
	    
	      # search for the belonging spike
	      test_idx_list = np.where(Data['trainNeuron'] == DataRef_spike)[0]	# get pos. of neurons AP in perturbed trajectory
	      test_Data_spiketimes = Data['trainTime'][test_idx_list]
	      test_Data_spiketimes_idx = abs(DataRef_spiketime-test_Data_spiketimes) < 10**(-1)
	      
	      search_list = list(test_Data_spiketimes[test_Data_spiketimes_idx])
	      
	      found = 0
	      while found == 0:
		if not len(search_list):
		  DataRef['deltaTime'][i] = 1
		  break
		
		test_time = min(search_list)
		test_idx = test_idx_list[np.where(test_time==test_Data_spiketimes)[0]]
		
		if test_idx in assigned:
		  search_list.remove(test_time)
		else:
		  DataRef['deltaTime'][i] = test_time
		  assigned.append(test_idx)
		  found = 1
		  
	  DataRef['deltaPhase'], DataRef['trainNeuron'] = sortTrain(DataRef,N,sort_it)	  # sort neuron indices according to activity (if wanted - very costly)

	  if DataRef['distance'][-1] >=5 and len(DataRef['decorr']):	# look only at diverged trajectories
	    idx_decorr[idx_sim] = 0
	    found = 0
	    while not found:
	      if any(DataRef['distance'][np.where(DataRef['measureTimes'][:idx_phase]>=DataRef['decorr'][idx_decorr[idx_sim]])[0]] <= 0.1):
		idx_decorr[idx_sim] += 1
	      else:
		try:
		  time_div[idx_sim] = DataRef['decorr'][0]
		  vel_div[idx_sim] = DataRef['distance'][-1]/(DataRef['trainTime'][idx_train]-DataRef['decorr'][idx_decorr[idx_sim]])
		  found = 1
		except:
		  break
	  
	  if numread==1:# or idx_decorr[idx_sim]:
	    plot_pert_spread(DataRef,idx_phase,idx_train,net,sort_it)
	  
	  idx_sim += 1
	  if not (idx_sim%100):
	    print idx_sim
	    
    ncid.close()
    
    mask = np.array(time_div!=0)
    
    time_mean = np.mean(time_div[mask])
    time_var = np.sqrt(np.sum((time_div[mask]-time_mean)**2)/sum(mask))
    
    vel_mean = np.mean(vel_div[mask])
    vel_var = np.sqrt(np.sum((vel_div[mask]-vel_mean)**2)/sum(mask))
    
    ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','w')
    
    ncid.createDimension('one',1)
    ncid.createDimension('sims',numread)

    VarP = ncid.createVariable('p','d', ('one',))
    
    VarTime = ncid.createVariable('time_div','d', ('sims',))
    VarVel = ncid.createVariable('velocity_div','d', ('sims',))
    VarSwitch = ncid.createVariable('switch','i', ('sims',))
    
    VarTime_mean = ncid.createVariable('time_mean','d',('one',))
    VarTime_var = ncid.createVariable('time_var', 'd', ('one',))
    VarVel_mean = ncid.createVariable('velocity_mean', 'd', ('one',))
    VarVel_var = ncid.createVariable('velocity_var', 'd', ('one',))
    VarSwitch_mean = ncid.createVariable('switch_mean', 'd', ('one',))
    
    VarP[:] = np.sum(mask)/float(numread)
    VarTime[:] = time_div
    VarVel[:] = vel_div
    VarSwitch[:] = idx_decorr
    VarTime_mean[:] = time_mean
    VarTime_var[:] = time_var
    VarVel_mean[:] = vel_mean
    VarVel_var[:] = vel_var
    VarSwitch_mean[:] = np.sum(idx_decorr[mask]!=0)
    
    ncid.close()
    
    erase = 1
    if erase and numread>1:
      cleanup_pert(net)
    
  time_CV = time_var/time_mean
  print 'time_div = %5.3f+-%5.3f' % (time_mean,time_var)
  print 'time CV = %5.3f' % time_CV
  
  vel_CV = vel_var/vel_mean
  print 'velocity = %5.3f+-%5.3f' % (vel_mean,vel_var)
  print 'velocity CV = %5.3f' % vel_CV
  
  #mask = np.array(time_div!=0)
  
  #plt.figure()
  #plt.hist(time_div[mask],bins=41,range=[0,0.1])
  #plt.show(block=False)
  
  #plt.figure()
  #plt.hist(vel_div[mask],bins=41,range=[0,3000])
  #plt.show(block=False)
  
  return time_div, vel_div, time_mean, vel_mean



def cleanup_pert(net):
  if erase:
    print "cleaning up..."
    for dat in os.listdir(net.Path['script']):
      if net.Hash in dat:
	for fileName in os.listdir(net.Path['script'] + dat):
	  if not ('processed' in fileName):
	    os.remove(net.Path['script'] + dat + '/' + fileName)



def plot_pert_spread(Data,idx_phase,idx_train,net,sort_it):
  ### here comes the plotting
  
  N = net.const['N']
  # setting up the color map
  levs = range(160)
  rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors=[(0,(0,0,0)),(1,(1,0,0))],N=len(levs)-1,)
  
  fig,axes = plt.subplots(nrows=3,ncols=1)

  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  
  if sort_it and sum(Data['trainNeuron']==N/2) == 0:
    train_max = N/2
  else:
    train_max = N

  #print idx_train
  #print len(Data['trainTime'])
  
  im0 = axes[0].imshow(Data['deltaPhase'], extent=[0,Data['deltaPhase'].shape[1],0,Data['deltaPhase'].shape[0]],aspect='auto',origin='lower',interpolation='none',cmap=rwb,vmin=0,vmax=0.5)
  axes[0].set_xlim([0,Data['trainTime'][idx_train-1]])
  axes[0].set_xticks(np.linspace(0,idx_phase,5))
  axes[0].set_xticklabels([])
  axes[0].tick_params(axis='both', which='major', labelsize=14)
  #axes[0].set_ylabel('Neuron Index',fontsize=18)
  axes[0].set_yticks(np.linspace(0,N,3))
  axes[0].set_ylim([0,train_max])

  divider0 = make_axes_locatable(axes[0])
  cax0 = divider0.append_axes("right", size="5%", pad=0.2)
  cbar0 = plt.colorbar(im0, cax=cax0)
  cbar0.ax.set_ylabel(r'$\displaystyle\Delta\theta$',fontsize=18)
  cbar0.set_ticks(np.linspace(0,0.5,3))
  cbar0.set_ticklabels(np.linspace(0,0.5,3))
  cbar0.ax.tick_params(labelsize=12)
  
  
  time_max = 10**(-1)
  #time_max = 5*10**(-4)
  im1 = axes[1].scatter(Data['trainTime'][:idx_train],Data['trainNeuron'][:idx_train],c=Data['deltaTime'],cmap=rwb,vmin=0,vmax=time_max,marker='D')
  axes[1].set_xlim([0,Data['trainTime'][idx_train-1]])
  axes[1].set_ylim([0,train_max])
  axes[1].tick_params(axis='both', which='major', labelsize=14)
  axes[1].set_ylabel('Neuron Index',fontsize=18)
  axes[1].set_xticks(np.linspace(0,Data['trainTime'][idx_train-1],5))
  axes[1].set_xticklabels([])
  axes[1].set_yticks(np.linspace(0,N,3))
  
  divider1 = make_axes_locatable(axes[1])
  cax1 = divider1.append_axes("right", size="5%", pad=0.2)
  cbar1 = plt.colorbar(im1, cax=cax1)
  cbar1.ax.set_ylabel(r'$\displaystyle\Delta t^{p}$ in ms',fontsize=18)
  cbar1.set_ticks(np.linspace(0,time_max,3))
  cbar1.set_ticklabels(np.linspace(0,time_max*1000,3))
  cbar1.ax.tick_params(labelsize=12)


  axes[2].plot(Data['measureTimes'][:idx_phase],Data['spread'],'k--',label=r'\# disturbed neurons')
  axes[2].plot(Data['measureTimes'][:idx_phase],Data['distance']/np.max(Data['distance']),'k-',label='distance from reference')
  axes[2].plot(Data['decorr'],np.zeros(len(Data['decorr'])),'kD',markersize=8,label='spike crossing events')
  axes[2].set_xlim([0,Data['trainTime'][idx_train-1]])
  axes[2].tick_params(axis='both', which='major', labelsize=14)
  axes[2].set_xlabel('time in s',fontsize=18)
  axes[2].set_xticks(np.linspace(0,round(Data['trainTime'][idx_train-1]*10)/10,5))
  #axes[2].set_xticklabels(np.linspace(0,t_ana,5))
  axes[2].set_ylim([0,1])
  axes[2].set_yticks(np.linspace(0,1,3))
  #axes[2].set_xticks(np.linspace(0,DataRef['trainTime'][idx_train],5))
  #axes[2].set_ylabel(r'\# disturbed neurons / distance')

  leg_loc = 2
  
  if len(Data['decorr']):
    if Data['decorr'][0] < 0.1:
      leg_loc = 4
      
  axes[2].legend(loc=leg_loc,prop={'size':12})
  divider2 = make_axes_locatable(axes[2])
  cax2 = divider2.append_axes("right", size="5%", pad=0.2)
  cax2.axis('off')
  
  #plt.suptitle(r'$\displaystyle \alpha_{recip}=%3.1f,\, \alpha_{conv}=%3.1f,\, \alpha_{div}=%3.1f,\, \alpha_{chain}=%3.1f$' % (net.const['alpha_recip'],net.const['alpha_conv'],net.const['alpha_div'],net.const['alpha_chain']),fontsize=20)

  #save_str = './pics/pert_spread/phase_%3.1f_%3.1f_%3.1f_%3.1f_N=%d.png' % (net.const['alpha_recip'],net.const['alpha_conv'],net.const['alpha_div'],net.const['alpha_chain'],net.const['N'])
  #plt.savefig(save_str)
  plt.show(block=False)



def sortTrain(Data,N,sort_it):

  if sort_it:
  # order trains:
    sort = {}
    train_num = [(n,sum(Data['trainNeuron']==n)) for n in range(N)]	# get number of spikes in window
    train_num = sorted(train_num, key=lambda tup: tup[1])		# order by number of spikes
    deltaPhase_sorted = np.zeros(Data['deltaPhase'].shape)
    for n in range(N):
      sort[train_num[n][0]] = N-1-n		# write replacement spike index into dictionary (highest first)
      deltaPhase_sorted[N-1-n] = Data['deltaPhase'][train_num[n][0]]		# rearrange phase measurement accordingly
    
    train_len = len(Data['trainNeuron'])
    trainNeuron_sorted = np.zeros(train_len)
    for t in range(train_len):
      trainNeuron_sorted[t] = sort[Data['trainNeuron'][t]]
  else:
    deltaPhase_sorted = Data['deltaPhase']
    trainNeuron_sorted = Data['trainNeuron']
  
  return deltaPhase_sorted, trainNeuron_sorted
  #plt.figure()
  #plt.scatter(Data['trainTime'][:idx_train],trainNeuron_sorted[:idx_train])
  #plt.show(block=False)




  
  
################## file accessment and reading ##################

def read_from_info(infoName):
  data_read = open(infoName,'r')
  para, topo, step = None, None, None
  for line in data_read:
    if 'para' in line[:15]:
      para = read_values(line,'=',[0,1])
    elif 'topo' in line[:15]:
      topo = read_values(line,'=',[0,1])
    elif 'step' in line[:15]:
      step = read_values(line,'=',[0,1])
    elif 'iter' in line[:15]:	# old file structure has to be taken into account
      step = read_values(line,'x',[1,0])
    
    elif 'Hash' in line:
      Hash = ''
      read_status = 0
      for char in line:
	if char == '<':
	  read_status = 1
	elif char == '>':
	  read_status = 0
	elif read_status == 1:
	  Hash += char

  data_read.close()
  
  return para, topo, step, Hash
  
  
def read_values(line,splittag,idx):
  vals = {}
  string = ''
  read_status = 0
  for char in line:
    if char == '<':
      read_status = 1
    elif char == '>':
      read_status = 0
    elif read_status == 1:
      string += char
  for val_pair in string.split('\t'):
    if val_pair != '':
      get_val = val_pair.split(splittag)
      try:
	if '.' in get_val[idx[1]]:
	  vals[get_val[idx[0]]] = float(get_val[idx[1]])
	else:
	  vals[get_val[idx[0]]] = int(get_val[idx[1]])
      except:
	vals[get_val[idx[0]]] = get_val[idx[1]]
  return vals


def pick_script(directory):
  idx = 0
  string = []
  for dat in os.listdir(directory):
    if not ('info' in dat):
      print '\t' + str(idx) + ' - ' + dat
      string.append(dat)
      idx += 1
  script_idx = int(raw_input('Pick a script: '))
  
  return string[script_idx]


def pick_file(scriptName, infoPath, file_excl='', restrictions=None):
  idx = 0
  fileName = []
  for dat in os.listdir(infoPath):
    invalid = 0
    if (((scriptName + '_para') in dat) and ('SimInfo' in dat) and (not (dat==file_excl))):
      para, topo, step, Hash = read_from_info(infoPath + dat)
      string = str(idx) + ' - '
      for j in range(len(para)):
	string += '\t' + para.keys()[j] + '=' + str(para.values()[j])
      if topo:	# does not have to be in
	for j in range(len(topo)):
	  string += '\t' + topo.keys()[j] + '=' + str(topo.values()[j])
      for j in range(len(step)):
	string += '\t' + step.keys()[j] + '=' + str(step.values()[j])
      string += '\n'
      #print "he"
      if restrictions:
	for item in restrictions.items():
	  if item[0] in para.keys():
	    if not (para[item[0]]==item[1]):
	      invalid = 1
	      break
	  elif item[0] in topo.keys():
	    if not (topo[item[0]]==item[1]):
	      invalid = 1
	      break
      if not invalid:
	print string
	fileName.append(dat)
	idx += 1
      
  ch_idx = int(raw_input('Pick a file: '))
  return fileName[ch_idx]


def clean_erase():
  scriptPath = 'data/SONET/'
  
  for dat in os.listdir(scriptPath):
    #print dat
    if 'results' in dat:
      for fileName in os.listdir(scriptPath + dat):
	if not ('results_processed' in fileName):
	  #print scriptPath + dat + '/' + fileName
	  os.remove(scriptPath + dat + '/' + fileName)
      #break


def erase(scriptPath=None,infoName=None,Hash=None,restrictions=None):
  if infoName and scriptPath:
    tmp_para, tmp_topo, tmp_stp, Hash = read_from_info(infoName)
    os.remove(infoName)
    for dat in os.listdir(scriptPath):
      if Hash in dat:
	shutil.rmtree(scriptPath + '/'+ dat,ignore_errors=True)
  else:
    scratch = int(raw_input('You want to access data from local (0) or from scratch (1)? '))
    if scratch==1:
      sessionPath = '/scratch01/aschmidt/data/'	#path for storage on scratch
    else:
      sessionPath = 'data/'			#path for local storage
    
    scriptName = pick_script(sessionPath)
    print scriptName
    not_satisfied = 'y'
    while not_satisfied == 'y':
      
      fileName = pick_file(scriptName,sessionPath + 'info/',restrictions=restrictions)
      infoName = sessionPath + 'info/' + fileName
      net = network('r',infoName)
      sim = simulation_data('r',infoName,net)
      #net.setup(sim,check=0,call=1)
      net.scriptName = scriptName
      
      
      os.remove(sessionPath + 'info/' + fileName)
      for dat in os.listdir(sessionPath + scriptName):
	if net.Hash in dat:
	  shutil.rmtree(sessionPath + scriptName + '/' + dat,ignore_errors=True)
      not_satisfied = raw_input('Another one? ')