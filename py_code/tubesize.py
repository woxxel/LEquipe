import os, imp, hashlib, heapq, inspect, time, subprocess, sys, shutil
from scipy.io import netcdf
import numpy as np
#import pylab as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import scipy.special as spec
import scipy as sp
from numpy import complex
from scipy import stats
from scipy import signal
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
#try:
  #from paida import *
#except:
  #1


assert os.getcwd().split('/')[-1] == 'Program', 'Please run the program in the directory of LEquipe to properly set up the paths.'

imp.load_source('read_write_code','Code/py_code/read_write.py')
imp.load_source('file_management_code','Code/py_code/file_management.py')
imp.load_source('support_code', 'Code/py_code/support.py')
imp.load_source('analysis_code', 'Code/py_code/analysis.py')
imp.load_source('fit_func', 'Code/py_code/fit_functions.py')
imp.load_source('bootstrapping', 'Code/py_code/bootstrapping.py')
from read_write_code import *
from file_management_code import *
from support_code import *
from analysis_code import *
from fit_func import *
from bootstrapping import *

# Rainers Mail: Analyze several differences in kind of simulations:
#    perturbations in which representation? (what happens to threshold-crossing perturbations?)
#    simulation start at spike or time? (warmup stop when?) check spike times, 1st one when?
#    -> when is perturbation applied?
#
# next: understand stochastic driving neurons, implement (2 populations)
#
# connectivities / weights from learning networks -> stabilize?
#    -> run for some time, then take as connectivity
#    -> or implement learning in LEquipe to have it while simulating?


class simulation_data:
  
  def __init__(self,mode_input='t',infoName=None,net=None,cluster=None,mode='eft',TC=0.5):
    
    self.mode = mode
    
    if not cluster:
      self.cluster = {'scratch':0,'par':0,'q':0}		# clusterstuff
    else:
      self.cluster = cluster

    self.Paras = {}					# bib to write simulation parameters to
    
    if self.mode == 'LE':
      # calculate Lypanunov Spectra
      self.TC = 20
      self.Paras['SWONS'] = 1
      #if net.topo == 'p':
	#self.Paras['ONstep'] = min(100,10*float(net.Const['N'])/net.Const['K']+net.topoConst['drive_rate']*net.Const['N'])
      #else:
      self.Paras['ONstep'] = 10*float(net.Const['N'])/net.Const['K']
      self.Paras['pLE'] = 1
	
      if net.Const['special'] == 0:
	if net.topo == 'p':
	  self.Paras['subLyapunovExp'] = net.Const['N']
	  self.Paras['subLyapunovExpConvergence'] = 1
	else:  
	  self.Paras['LyapunovExpConvergence'] = 1
	self.Paras['LyapunovExp'] = net.Const['N']	# Lyapunov Exponents
	
      elif net.Const['special'] == 1:
	if net.topo == 'p':
	  self.Paras['subLyapunovExp'] = 1
	  self.Paras['subLyapunovExpConvergence'] = 1
	else:
	  self.Paras['LyapunovExpConvergence'] = 1
	  
	self.Paras['LyapunovExp'] = 2	# Lyapunov Exponents
      elif net.Const['special'] == 2:
	self.TC = 2
	self.Paras['LyapunovExp'] = net.Const['N']
	self.Paras['CLV'] = 1				# Covariant Lyapunov Exponents / local Lyapunov Exp
	self.Paras['SWCLV'] = 2
	self.Paras['ONstep'] = 1
      
    
    if self.mode == 'statistics':			# no warmup, thus all measures are made in one simulation
      # calculate cv, skewness, etc
      
      if net.Const['special'] == 0:
	self.TC = 100
	self.Paras['measures'] = 2
	self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	self.Paras['train'] = range(net.Const['N'])
	self.Paras['synchrony'] = 1
	self.Paras['ISIneurons'] = range(net.Const['N'])			# neurons whose spike statistics are calculated
	self.Paras['ISIdecorr'] = range(net.Const['N'])					# ISI of potential decorrelation events are calculated for x neurons
	self.Paras['ISIstats'] = 2					# moments that are being calculated
	
      elif net.Const['special'] == 1:
	self.TC = 100
	self.Paras['train'] = range(net.Const['N'])
	self.Paras['measures'] = 2
	self.Paras['measureTimes'] = np.linspace(0,1,1000+1)
	self.Paras['synchrony'] = 1
	self.Paras['ISIneurons'] = range(net.Const['N'])
	self.Paras['ISIstats'] = 2
	
      #else:
	#self.TC = 1000.
	##self.Paras['measures'] = 2
	#self.Paras['train'] = range(net.Const['N'])
	#self.Paras['synchrony'] = 1
	#self.Paras['ISIneurons'] = range(net.Const['N'])			# neurons whose spike statistics are calculated
	#self.Paras['ISIstats'] = 2					# moments that are being calculated
	##self.Paras['ISIdecorr'] = range(net.Const['N'])					# ISI of potential decorrelation events are calculated for x neurons
    
    
    if self.mode == 'eft':
      # calculate mean flux tube radius
      self.TC = 0.5
      #self.TC = TC
      self.Paras['CDB'] = 2
      self.Paras['measureTimes'] = np.linspace(0,self.TC,51)		# not valid for warmup
      #self.Paras['measures'] = 2
      #self.Paras['train'] = range(net.Const['N'])
      #self.Paras['CDB'] = 0
    
    
    if self.mode == 'eft_tmp':
      # calculate mean flux tube radius
      self.TC = 2
      if net.topo == 'p':
	self.Paras['CDB'] = 3						# only valid for self.simulation
      else:
	self.Paras['CDB'] = 4
      #self.Paras['measureTimes'] = np.linspace(0,TC,11)		# not valid for warmup
      #self.Paras['measures'] = 2
    
      
    if self.mode == 'mosaic':
      # calculate "mosaic" cross-section of fluxtubes
      self.TC = 0.4
      self.Paras['train'] = range(net.Const['N'])
      self.Paras['measureTimes'] = np.linspace(0,TC,101)

    
    if self.mode == 'samples':
      assert 'special' in net.Const.keys(), 'If you calculate samples, you have to specify a "special"!'
      # simulate sample trajectories to e.g. track "cloud" behaviour
      
      if net.Const['special'] == 0:		# large batch for obtaining conv/div statistics
	self.TC = 2
	self.Paras['measures'] = 2
	if net.Const['N'] > 1000:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	else:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
	#self.Paras['CDB'] = 2
	
      if net.Const['special'] == 1:		# few trajectories for plotting
	self.Paras['measures'] = 2
	if net.Const['N'] > 1000:
	  self.TC = 0.5
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	else:
	  self.TC = 5
	  self.Paras['train'] = range(net.Const['N'])
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)	#self.Paras['measureSpikes'] = range(int(self.TC*net.Const['N']*net.Const['rateWnt']))
      
      if net.Const['special'] == 2:		# few for obtaining D_decorr
	self.TC = 0.2
	self.Paras['measures'] = 2
	if net.Const['N'] > 1000:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	else:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
	
      if net.Const['special'] in [3,30]:		# obtain border to totally stable networks
	if net.Const['special'] == 3:
	  self.TC = 1
	elif net.Const['special'] == 30:
	  self.TC = TC
	#self.Paras['CDB'] = 1			# no breaking, as this is all-to-all comparison
	self.Paras['measures'] = 2		# measures to get phase distribution over time
	self.Paras['measureTimes'] = np.append(0,10**np.linspace(-3,np.log10(self.TC),100))# np.append(np.linspace(0,0.099,100),np.linspace(0.1,self.TC,1000))
	#self.Paras['train'] = range(net.Const['N'])
    
      if net.Const['special'] == 4:		# add single spike failure
	self.TC = 1
	self.Paras['measures'] = 2		# measures to get phase distribution over time
	self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
	self.Paras['train'] = range(net.Const['N'])
      
      if net.Const['special'] == 5:
	self.TC = 10
	self.Paras['measures'] = 2		# measures to get phase distribution over time
	self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
      
      if net.Const['special'] == 6:
	self.TC = 1
	self.Paras['measures'] = 2
	self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
	self.Paras['train'] = range(net.Const['N'])
	
	self.Paras['pertSpike'] = 1
      
      if net.Const['special'] == 7:		# get pseudo Lyapunov Exponents
	self.TC = 0.05
	self.Paras['measures'] = 2
	if net.Const['N'] > 1000:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	else:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
	self.Paras['train'] = range(net.Const['N'])
      
      if net.Const['special'] == 8:
	self.Paras['train'] = range(net.Const['N'])
	self.TC = 1
	self.Paras['measures'] = 2
	if net.Const['N'] > 1000:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,100+1)
	else:
	  self.Paras['measureTimes'] = np.linspace(0,self.TC,1000+1)
      #print "ToDo: singlespike failure"
	
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
	if net.Const['special'] == 0:
	  self.steps = {'ps':1, 'co':10,'in':1,'pd':1}			# get 10 different sims
	if net.Const['special'] == 1:
	  self.steps = {'ps':1, 'co':1,'in':1,'pd':1}
	  
      if self.mode == 'eft':
	#self.steps = {'ps':21, 'co':2,'in':5,'pd':10}			# iteration steps
	if 'special' in net.Const.keys():
	  if net.Const['special'] == -2:	# for test runs
	    self.steps = {'ps':11, 'co':2,'in':5,'pd':100}
	    
	  elif net.Const['special'] == -1:	# for test runs (with larger number)
	    self.steps = {'ps':15, 'co':5,'in':10,'pd':100}
	    
	  elif net.Const['special'] == 0:
	    if net.topoConst['alpha_conv'] > 0:
	      self.steps = {'ps':31, 'co':50,'in':50,'pd':100}
	      #self.TC = 2
	    else:
	      #self.steps = {'ps':31, 'co':10,'in':50,'pd':100}
	      self.steps = {'ps':31, 'co':5,'in':20,'pd':50}
	    
	  elif net.Const['special'] == 1:				# get convergence
	    self.steps = {'ps':31, 'co':20,'in':50,'pd':200}	# iteration steps (dont change!)
	    

	else:
	  #if net.topoConst['alpha_conv'] > 0:
	    #self.steps = {'ps':31, 'co':20,'in':50,'pd':100}	# iteration steps
	  #else:
	  self.steps = {'ps':31, 'co':20,'in':50,'pd':100}	# iteration steps
	
      
      if self.mode == 'eft_tmp':
	if 'special' in net.Const.keys():
	  if net.Const['special'] == -1:	# for test runs
	    self.steps = {'ps':13, 'co':10,'in':5,'pd':50}
	  else:
	    if net.topoConst['alpha_conv'] > 0:
	      self.steps = {'ps':13, 'co':20,'in':50,'pd':100}
	      #self.TC = 2
	    else:
	      self.steps = {'ps':13, 'co':10,'in':50,'pd':100}		# iteration steps
	#self.steps = {'ps':9, 'co':2,'in':5,'pd':10}			# iteration steps
	
      if self.mode == 'mosaic':
	n=201				# decent value to get nice image quickly
	#if not (mode_input == 'r'):
	  #n = int(raw_input('How many samples across one dimension? '))	# number of bins across one dimension
	self.steps = {'ps':1, 'co':1,'in':1,'pd':n**2}
      
      if self.mode == 'samples':
	if net.Const['special'] == 0:			# get convergence statistics
	  self.steps = {'ps':1, 'co':2,'in':10,'pd':10}
	if net.Const['special'] == 1:			# trajectories for plotting
	  self.steps = {'ps':1, 'co':2,'in':5,'pd':10}
	
	if net.Const['special'] == 2:			# trajectories for obtaining decorr distance
	  #if cluster['scratch']:
	  self.steps = {'ps':1, 'co':5,'in':5,'pd':10}
	  #else:
	    #self.steps = {'ps':1, 'co':2,'in':2,'pd':5}
	if net.Const['special'] == 3:			# obtain border to totally stable networks
	  if cluster['scratch']:
	    self.steps = {'ps':1, 'co':2,'in':20,'pd':50}
	  else:
	    self.steps = {'ps':1, 'co':2,'in':5,'pd':20}
	if net.Const['special'] == 30:			# obtain border to totally stable networks
	  if cluster['scratch']:
	    self.steps = {'ps':1, 'co':2,'in':50,'pd':50}
	  else:
	    self.steps = {'ps':1, 'co':2,'in':5,'pd':20}
	#if net.Const['special'] in [2,3]:
	  #assert self.steps['in'] == 1, '"in" should always be 1, as this tests only decorrelated trajectories!'	
	if net.Const['special'] == 4:			# samples for random initial conditions
	  self.steps = {'ps':1, 'co':2,'in':5,'pd':10}
	if net.Const['special'] == 5:			# samples for random initial conditions (understand distance statistics)
	  self.steps = {'ps':1, 'co':10,'in':1,'pd':10}
	if net.Const['special'] == 6:
	  self.steps = {'ps':1, 'co':1,'in':10,'pd':1}
	if net.Const['special'] == 7:
	  self.steps = {'ps':1, 'co':2,'in':5,'pd':10}
	if net.Const['special'] == 8:
	  self.steps = {'ps':3, 'co':5,'in':5,'pd':10}
	
      if self.mode == 'drive_cur':
	self.steps = {'ps':1, 'co':3,'in':1,'pd':1}
      #if self.mode == 'pert_spread':
	#self.steps = {'ps':1, 'co':5,'in':10,'pd':20}
      
    
      self.steps['total'] = self.steps['pd']*self.steps['in']*self.steps['co']
    

class network(read_write):
  
  def __init__(self,mode_input='t',infoName=None,alpha=None,hPVars=None,hPcorr=1,netConst={'N':1000,'K':100,'rateWnt':10,'J0':-1,'NeuronType':2}):
    
    self.Const = {}	# write constants of the network into this
    self.topoConst = {}
    #self.puppet = 0
    
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
	self.topoConst['corr'] = hPcorr
	
      #if type(puppet[0]) == str:
	#self.topo = puppet[0][0]
	#self.puppet = 1
	#self.scriptName = puppet[0]
	#self.topoConst['drive_type'] = puppet[0]
	#for key in puppet[1].keys():
	  #self.topoConst['drive_'+key] = puppet[1][key]

      for key in netConst.keys():
	self.Const[key] = netConst[key]
      if self.topo == 'p':
	self.Const['Ndrive'] = self.Const['N']/self.topoConst['corr']
      
      self.Const['tauM'] = 10.
      
    
    
    self.rateWnt = self.Const['rateWnt']


  def setup(self,sim,perturbation=None,sv_file=[1,1,1,1,1],hide=1,call=0):
    
    if not 'corr' in self.topoConst.keys():
      self.topoConst['corr'] = 1
    # write to simulation dictionaries
    self.ParaNet, self.ParaSim = {}, {}
    self.break_it = 0				# breaking status
    sim.state = 0
    # read-in of given network/neuron parameters
    for items in self.Const.items():
      self.ParaNet[items[0]] = items[1]
    self.perturbation = perturbation
    
    # procession of network/neuron parameters
    self.ParaNet['Iext'] = -(self.Const['J0']*self.rateWnt*np.sqrt(self.Const['K']))*self.Const['tauM']/1000.
    if self.Const['rateWnt'] < 10:	#choose smaller starting guess to 
      self.ParaNet['Iext'] /= 100.
    if self.topo == 'p':
      self.ParaNet['Iext'] -= min(self.Const['rateWnt']*self.Const['K'],self.topoConst['drive_rate'])*max(-10,self.topoConst['drive_cplg'])*self.Const['tauM']/1000.
      
      
      self.ParaSim['rateWnt'] = 0				# unfix rate of whole network
      self.ParaSim['rateWntSubN'] = self.rateWnt
      self.ParaSim['subN'] = self.Const['N']
      if self.topoConst['corr'] > 1:
	self.ParaNet['N'] = self.Const['N'] + self.Const['Ndrive']
	print "network size: ", self.ParaNet['N']
	self.ParaNet['NeuronType'] = np.append(self.Const['NeuronType']*np.ones(self.Const['N']),5*np.ones(self.Const['Ndrive']))
	self.ParaNet['poissonrate'] = np.append(np.zeros(self.Const['N']),self.topoConst['drive_rate']*np.ones(self.Const['Ndrive']))
	self.ParaNet['Iext'] = np.append(self.ParaNet['Iext']*np.ones(self.Const['N']),self.topoConst['drive_rate']*np.ones(self.Const['Ndrive']))
      else:
	self.ParaNet['N'] = self.Const['N']*2			# extend network by N poisson-neurons
	self.ParaNet['NeuronType'] = np.kron((self.Const['NeuronType'],5),np.ones(self.Const['N']))
	self.ParaNet['poissonrate'] = np.kron([0,self.topoConst['drive_rate']],np.ones(self.Const['N']))
	self.ParaNet['Iext'] = np.kron([self.ParaNet['Iext'],self.topoConst['drive_rate']],np.ones(self.Const['N']))
      
      if 'train' in sim.Paras.keys():
	sim.Paras['train'] = range(self.ParaNet['N'])
      
    else:
      self.ParaSim['rateWnt'] = self.rateWnt
    
    #if sim.mode in ['eft']:
      #self.ParaSim['D_decorr'] = 10
    if 'CDB' in sim.Paras:# or (sim.mode == 'samples' and self.Const['special'] == 0):
      try:
	if self.analysis:
	  self.ParaSim['D_decorr'] = self.analysis['D_decorr'][0]
	else:
	  if self.Const['N'] >= 10000:
	    self.ParaSim['D_decorr'] = 10
	  else:
	    self.ParaSim['D_decorr'] = np.nan
	
	if ((self.ParaSim['D_decorr'] < 2) or np.isnan(self.ParaSim['D_decorr'])):
	  if sim.mode == 'eft_tmp':
	    sim.Paras['CDB'] = 3
	  else:
	    sim.Paras['CDB'] = 1
	  print "Not using any upper CDB border."
	else:
	  print "Using %g as upper CDB border." % self.ParaSim['D_decorr']
      except:
	if sim.mode == 'eft_tmp':
	    sim.Paras['CDB'] = 3
	else:
	  sim.Paras['CDB'] = 1
	print "Not using any upper CDB border."
   
    #print self.analysis.keys()
      
    #if sim.mode in ['eft']:
      #if 'tau_conv' in self.analysis.keys():
	##print self.analysis['tau_conv']
	#sim.TC = max(0.5,min(5,np.log(10**(-8))*self.analysis['tau_conv'][1]))
	#print np.log(10**(-8))*self.analysis['tau_conv'][1]
	
	#print "Time duration for simulation (TC = %.2gs) obtained from characteristic convergence time (max = 5s)." % sim.TC
      #elif self.Const['K']*self.Const['rateWnt'] >= 1000:
	#sim.TC = 0.5
      #else:
	#sim.TC = 2
      
      sim.Paras['measureTimes'] = np.linspace(0,sim.TC,51)
    
    
    # setting up data structure
    if not (sim.mode in ['eft']):				# so that directory name includes script name
      self.scriptName += '_' + sim.mode
    self.data_structure(sim,hide=hide,call=call)				# set data structure (names, files, folders)
    if self.break_it or call:					# if something went wrong
      return self.break_it
    #if 'CDB' in sim.Paras.keys():
    
    #try:
      
    ### create topologies and get the respective driving currents to obtain the wanted firing rates ###
   
    if sim.mode in ['eft','eft_tmp']:
      print 'Starting ratefinding for setting external currents...'
    
    time_start = time.time()
    
    # set up storing file for topologies
    datTop_sv = netcdf.netcdf_file(self.Path['inTopo_links'],'w')
    datTop_sv.createDimension('Simulations',sim.steps['co'])
    datTop_sv.createDimension('links',2)
    datTop_sv.createDimension('linksize',len(self.Path['inpara'])+52)
    saveLinks = datTop_sv.createVariable('rateFinding','S1', ('Simulations','links','linksize'))
    
    # get different topologies with their respective external currents
    self.ParaSim['TR'] = 2
    #self.ParaSim['SR'] = 2

    # set Paths and create .nc files for simulations
    self.writeNet()
    self.writePuppet()
    self.writeSim()
    
    os.mkdir(self.Path['results'] + 'init/')
    
    sim.steps['ct'] = 1
    for sim.steps['co_ct'] in range(sim.steps['co']):
      self.ParaTopo = {}
      
      if sim.steps['co'] > 1:
	print_status(sim,time_start)
      
      self.create_topology(hide)
      
      if self.break_it:			# if topology could not be set up
	return self.break_it
      
      self.writeTopo()
      
      HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo'],self.Path['Puppet']])).hexdigest()	# construct and ...
      self.Path['Out'] = self.Path['results'] + 'init/con-' + HashDataOut + '.nc'	# ... assign Out-Name
      
      #if sim.cluster['q']:
      self.Path['ID'] = '%s_R_%02d' % (self.scriptName[0],sim.steps['co_ct'])
      # run of simulation to either calculate specified values or ratefinding and warmup
      run_script(self.Path,sim,CDB=0,hide=hide)
      
      saveLinks[sim.steps['co_ct']] = np.array([list(self.Path['Topo']),list(self.Path['Out'])])
      
      sim.steps['ct'] += 1
    datTop_sv.close()
    
    if not sv_file[0] and not sim.cluster['q']:
      os.remove(self.Path['Net'])
    if not sv_file[1] and not sim.cluster['q']:
      os.remove(self.Path['Sim'])
    
    print ']'
        
    if sim.mode == 'drive_cur':
      return
    
    ### now warm up the network and get a reference trajectory ###
    if sim.mode in ['eft','eft_tmp']:
      print 'Warming up networks and calculating reference trajectories...'
    
    # set up storing file for reference trajectories
    datinit_sv = netcdf.netcdf_file(self.Path['inpara_links'],'w')
    datinit_sv.createDimension('connectivities',sim.steps['co'])
    datinit_sv.createDimension('initial conditions',sim.steps['in'])
    datinit_sv.createDimension('Simulations',sim.steps['co']*sim.steps['in'])
    datinit_sv.createDimension('links',4)
    datinit_sv.createDimension('linksize',len(self.Path['inpara'])+52)
    saveLinks = datinit_sv.createVariable('initial data','S1', ('connectivities','initial conditions','links','linksize'))
    
    datTop_rd = netcdf.netcdf_file(self.Path['inTopo_links'],'r',mmap=False)
    
    # warmup networks and get reference trajectories
    self.ParaSim['TR'] = 0
    self.ParaSim['TW'] = 1
    
    if sim.mode == 'samples':
      #if self.Const['special'] in [0,1]:
	#self.ParaSim['measures'] = 2				# to get reference phases
      if self.Const['special'] in [2,3,4,5,30]:			# get initial conditions from natural distributions
	self.ParaSim['TC'] = sim.steps['pd']
	self.ParaSim['measureTimes'] = np.linspace(0,self.ParaSim['TC'],sim.steps['pd'])
	self.ParaSim['measures'] = 2
	
	#self.ParaSim['measureTimes'] = np.array(list(self.ParaSim['measureTimes']) + list(np.arange(sim.steps['in']*sim.steps['pd'])+self.ParaSim['TC']))
	#self.ParaSim['TC'] += sim.steps['in']*sim.steps['pd']
	#print self.ParaSim['TC']
      
      #if self.topo == 'p':
	#self.ParaSim['train'] = range(self.ParaNet['N'][0])	# to get spike times of poisson neurons (for puppets)

    if sim.mode == 'LE':
      #if self.topo == 'p':
	
      for items in sim.Paras.items():
	print items
	self.ParaSim[items[0]] = items[1]
	#self.ParaSim['LyapunovExp'] = self.ParaNet['N'][0]
      self.ParaSim['TC'] = sim.TC
    
    if sim.mode == 'statistics':
      self.ParaSim['TC'] = sim.TC
      for items in sim.Paras.items():
	self.ParaSim[items[0]] = items[1]
	  #print items
      
    self.writeSim()
    
    sim.steps['ct'] = 1
    
    for sim.steps['co_ct'] in range(sim.steps['co']):
      
      # read in topology and strength of external driving current (from simulation)
      [self.Path['Topo'],rateFinding] = [''.join(dat) for dat in datTop_rd.variables['rateFinding'][sim.steps['co_ct']][:]]
      Data = wait_for_godot(rateFinding)
      
      if self.topo == 'p':
	#self.ParaNet['Iext'] = np.kron([Data['finalCurrents'][0],self.topoConst['drive_rate']],np.ones(self.Const['N']))
	self.ParaNet['Iext'] = np.append(Data['finalCurrents'][0]*np.ones(self.Const['N']),self.topoConst['drive_rate']*np.ones(self.Const['Ndrive']))
      else:
	self.ParaNet['Iext'] = Data['finalCurrents']
      
      for sim.steps['in_ct'] in range(sim.steps['in']):	# change initial conditions
	if (sim.steps['co'] > 1) or (sim.steps['in'] > 1):
	  print_status(sim,time_start)
	
	trash = self.ParaNet.pop('init')		# remove initial conditions to create new ones
	  
	# set Paths and create .nc files for simulations
	self.writeNet()
	self.writePuppet()
	
	HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()	# construct and ...
	self.Path['Out'] = self.Path['results'] + 'init/ini-' + HashDataOut + '.nc'	# ... assign Out-Name
	
	if sim.cluster['q']:
	  self.Path['ID'] = '%s_W_%02d' % (self.scriptName[0],sim.steps['co_ct'])
	  
	# run of simulation to either calculate specified values or ratefinding and warmup
	if sim.mode in ['LE','statistics'] and sim.cluster['par']:
	  run_script_stats(self.Path,sim,CDB=0,hide=hide)
	#elif sim.cluster['q']:
	  #run_script_cluster(self.Path,sim,CDB=self.ParaSim['CDB'][0],loop='in',hide=hide)
	else:
	  run_script(self.Path,sim,CDB=0,hide=hide)
	  
	sim.steps['ct'] += 1
	
	# save links in ncfile for setting references
	saveLinks[sim.steps['co_ct'],sim.steps['in_ct'],:,:] = np.array([list(self.Path['Net']),list(self.Path['Sim']),list(self.Path['Topo']),list(self.Path['Out'])])
    
    datTop_rd.close()
    datinit_sv.close()
    
    if sim.mode in ['eft','eft_tmp']:
      print ']'
      print 'time elapsed: %g sec\n' % (time.time()-time_start)
      
    if sim.mode in ['eft','eft_tmp']:
      self.ParaSim['pd'] = sim.steps['pd']
    elif sim.mode == 'samples':
      if self.Const['special'] in [0,1,6,7,8]:
	self.ParaSim['pd'] = sim.steps['pd']
      
    
    # set parameters for simulation
    self.ParaSim['TW'] = 0
    
    # remove saving of phases and spikes
    trash = self.ParaSim.pop('measures')
    trash = self.ParaSim.pop('train')
    
    #sim.TC = 1
    self.ParaSim['TC'] = sim.TC
    

    for items in sim.Paras.items():
      self.ParaSim[items[0]] = items[1]
    
    if sim.mode == 'eft_tmp':
      self.ParaSim['measureTimes'] = np.linspace(0,sim.TC,int(sim.TC*10)+1)
      #self.ParaSim['measures'] = 2
      print 'measures at '
      print self.ParaSim['measureTimes']
    
    self.writeSim()	# this will not change anymore!
    #print self.ParaSim
    #if self.topo == 'p':
      #self.ParaNet['NeuronType'] = np.kron((self.Const['NeuronType'],7),np.ones(self.Const['N']))
    
    if sim.mode in ['eft','eft_tmp','pert_spread','mosaic','samples']:
      self.get_pert_range(sim)
    #print self.pert_size
      
    return self.break_it
    
  
# ----------------------------------- simulation stuff -----------------------------------------      
  
  def simulation(self,sim,sv_file=[1,1,1,1,1],check=0,hide=1):
    
    # initiate save arrays
    self.initDataSv(sim)
    sim.state = 1
    
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
    
    scaled_down = 0
    
    datinit_rd = netcdf.netcdf_file(self.Path['inpara_links'],'r',mmap=False)	# read input data from this
    
    pert_done = []
    
    idx_pert = 0
    
    puppet_str = []
    
    print 'Now the real simulations start (TC=%5.3g)' %sim.TC
    sim.steps['ct'] = 1
    
    for sim.steps['ps_ct'] in range(sim.steps['ps']):
      # setup time and save-stuff for this perturbation size
      time_start = time.time()
      
      if sim.mode in ['eft','eft_tmp']:
	self.ParaSim['pertSize'] = [self.pert_size[idx_pert]]
	pert_done.append(self.pert_size[idx_pert])
	print 'pertSize: %g' % self.ParaSim['pertSize'][0]
      
      elif sim.mode == 'samples':
	if self.Const['special'] in [0,1,7,8]:
	  self.ParaSim['pertSize'] = [self.pert_size[idx_pert]]
      
      self.writeSim()
      
      if sim.mode in ['eft','eft_tmp']:
	datout_sv,svLinks,svIDs = sv_links(sim,self.Path,self.pert_size[idx_pert])
      elif sim.mode == 'samples':
	if self.Const['special'] in [0,1,6,7,8]:
	  datout_sv,svLinks,svIDs = sv_links(sim,self.Path,self.pert_size[idx_pert])
	elif self.Const['special'] in [2,3,4,5,30]:
	  datout_sv,svLinks,svIDs = sv_links(sim,self.Path,self.pert_size[idx_pert],linkNr=3)
      pathPert = '%sP%02d' % (self.Path['results'],sim.steps['ps_ct'])
      
      if os.path.exists(pathPert):
	shutil.rmtree(pathPert,ignore_errors=True)
      
      os.mkdir(pathPert)
      
      for sim.steps['co_ct'] in range(sim.steps['co']):
	if not sim.steps['ps_ct']:
	  puppet_str.append([])
	
	if sim.mode in ['eft','eft_tmp']:
	  ID = sim.steps['co_ct']
	  os.mkdir('%sP%02d/%03d' % (self.Path['results'],sim.steps['ps_ct'],ID))
	  
	for sim.steps['in_ct'] in range(sim.steps['in']):
	  
	  if not (sim.mode in ['eft','eft_tmp']):
	    ID = sim.steps['in_ct'] + sim.steps['in']*sim.steps['co_ct']
	    os.mkdir('%sP%02d/%03d' % (self.Path['results'],sim.steps['ps_ct'],ID))
	  
	  [netPath_tmp,simPath_tmp,self.Path['Topo'],self.Path['Ref']] = [''.join(fileNames) for fileNames in datinit_rd.variables['initial data'][sim.steps['co_ct'],sim.steps['in_ct'],:,:]]

	  if self.topo == 'p':
	    ncid = netcdf.netcdf_file(netPath_tmp,'r',mmap=False)
	    initPoisson = ncid.variables['initStates'][:][self.Const['N']:]
	    ncid.close()
	  
	  Data = wait_for_godot(self.Path['Ref'])
	  
	  if self.topo == 'p':
	    self.ParaNet['Iext'] = np.append(Data['finalCurrents'][0]*np.ones(self.Const['N']),self.topoConst['drive_rate']*np.ones(self.Const['Ndrive']))
	    #self.ParaNet['Iext'] = np.kron([Data['finalCurrents'][0],self.topoConst['drive_rate']],np.ones(self.Const['N']))
	  else:
	    self.ParaNet['Iext'] = Data['finalCurrents']	# get driving current
	  
	  #if self.topo == 'p':
	    #if not sim.steps['ps_ct']:			# get spiketimes from Ref for puppetTrain
	      #self.get_puppetTrain(Data['trainNeuron'],Data['trainTime'],sim.TC)
	      #self.writePuppet()
	      #puppet_str[sim.steps['co_ct']].append(self.Path['Puppet'])
	    #else:					# and use those for other perturbations
	      #self.Path['Puppet'] = puppet_str[sim.steps['co_ct']][sim.steps['in_ct']]
	  
	  if sim.mode in ['eft','eft_tmp']:	# make reference and perturbed trajectory in one simulation
	    print_status(sim,time_start,pertSize=self.pert_size[idx_pert])
	    
	    self.ParaNet['init'] = Data['finalStates'][:]
	    if self.topo == 'p':
	      self.ParaNet['init'][self.Const['N']:] = initPoisson
	    
	    self.writeNet()
	    
	    HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()	# construct and ...
	    self.Path['Out'] = '%s/Data-%s.nc' % (pathPert,HashDataOut)			# ... assign Out-Name
	    
	    svLinks[sim.steps['co_ct'],sim.steps['in_ct'],:] = [list(self.Path['Topo']),list(self.Path['Out'])]	# archive file names
	    
	    self.Path['ID'] = '%s%s_P%d_%d' % (self.scriptName[0],self.Hash[:2],sim.steps['ps_ct'],ID)
	    svIDs[sim.steps['co_ct'],sim.steps['in_ct']] = ID
	    run_script(self.Path,sim,CDB=self.ParaSim['CDB'][0],hide=hide)	## run the C++ simulation
	    
	    sim.steps['ct'] += sim.steps['pd']
	    
	  elif sim.mode in ['samples']:
	    
	    #init = Data['measureStates'][0]		# get initial conditions
	    
	    if self.Const['special'] in [0,1,6,7,8]:
	      print_status(sim,time_start,pertSize=self.pert_size[idx_pert])
	      
	      self.ParaNet['init'] = Data['finalStates'][:]
	      if self.topo == 'p':
		self.ParaNet['init'][self.Const['N']:] = initPoisson
	      
	      self.writeNet()
	      
	      HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()	# construct and ...
	      self.Path['Out'] = '%sP%02d/%03d/D%s.nc' % (self.Path['results'],sim.steps['ps_ct'],ID,HashDataOut)			# ... assign Out-Name
	      
	      svLinks[sim.steps['co_ct'],sim.steps['in_ct'],:] = [list(self.Path['Topo']),list(self.Path['Out'])]	# archive file names (save just topo, perttraj and reftraj)
	      
	      self.Path['ID'] = '%s%s_P%d_%d' % (self.scriptName[0],self.Hash[:2],sim.steps['ps_ct'],ID)
	      svIDs[sim.steps['co_ct'],sim.steps['in_ct']] = ID
	      
	      run_script(self.Path,sim,CDB=self.ParaSim['CDB'][0],hide=hide)	## run the C++ simulation
	      
	      # end of pd_stp
	      sim.steps['ct'] += sim.steps['pd']
	    
	    elif self.Const['special'] in [2,3,4,5,30]:	# get decorrelated initial conditions from warmup run (natural distribution)
	      
	      for sim.steps['pd_ct'] in range(sim.steps['pd']):
		print_status(sim,time_start,pertSize=self.pert_size[idx_pert])
		
		self.ParaNet['init'] = Data['measureStates'][sim.steps['pd_ct']]
		if self.topo == 'p':
		  self.ParaNet['init'] = np.append(self.ParaNet['init'],initPoisson)
		
		self.writeNet()
		
		HashDataOut = hashlib.sha1(np.array([self.Path['Net'],self.Path['Sim'],self.Path['Topo']])).hexdigest()	# construct and ...
		self.Path['Out'] = '%sP%02d/%03d/D%s.nc' % (self.Path['results'],sim.steps['ps_ct'],ID,HashDataOut)		# ... assign Out-Name
		
		svLinks[sim.steps['co_ct'],sim.steps['in_ct'],sim.steps['pd_ct'],:] = [list(self.Path['Topo']),list(self.Path['Out']),list(self.Path['Ref'])]	# archive file names (save just topo, perttraj and reftraj)
		
		self.Path['ID'] = '%s%s_P%d_%d' % (self.scriptName[0],self.Hash[:2],sim.steps['ps_ct'],ID)
		svIDs[sim.steps['co_ct'],sim.steps['in_ct'],sim.steps['pd_ct']] = ID
		
		run_script(self.Path,sim,CDB=self.ParaSim['CDB'][0],hide=hide)	## run the C++ simulation
		
		trash = self.ParaNet['init']
		
		# end of pd_stp
		sim.steps['ct'] += 1
	      
	      
	      
	  # end of in_stp
	# end of co_stp
      print ']'
      datout_sv.close()
      if sim.mode in ['eft']:
	if sim.cluster['q']: # directly read results to promote fast finding
	  read_max=min(1000,0.1*sim.steps['total'])
	  correct_ct = 0
	  while correct_ct < read_max:
	    
	    handle = subprocess.Popen('qstat -u aschmidt | grep %s%s_P%d_* | wc -l'%(self.scriptName[0],self.Hash[:2],sim.steps['ps_ct']),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
	    out,err = handle.communicate()
	    out = int(out)
	    if not out:
	      time.sleep(30)
	    
	    p, correct_ct = self.inter_read(sim,read_max)	# when on cluster, dont process it
	    if correct_ct < read_max:
	      if not out:
		wait_str = 'No measurements in line anymore - break without interreading'
		p = 0.5
		correct=1
		break
	      else:
		wait_str = 'Not enough measurements were given (%d/%d), wait for more!' % (correct_ct,read_max)
		time.sleep(5)
	      sys.stdout.write('\r{0}'.format(wait_str))
	      sys.stdout.flush()
	      
	      #print '%s_P_%d_* | wc -l' % (self.scriptName[0],sim.steps['ps_ct'])
	      
	      
	  # now start the real reading
	  if not sim.steps['ps_ct']:
	    self.Path['read_cluster'] = self.Path['session'] + 'read/'
	    if not os.path.exists(self.Path['read_cluster']):
	      os.mkdir(self.Path['read_cluster'])
	    runIt(sim.cluster,'python Code/py_code/read_cluster.py %s %d' % (self.Path['info_file'],sim.cluster['scratch']),self.Path['read_cluster'],name='read_%s'%self.Hash)
	
	else:
	  self.prob_read(sim,steps_rd=[sim.steps['ps_ct']])
	  p = self.DataSv['p'][sim.steps['ps_ct']]
	  correct_ct = self.DataSv['correct_ct'][sim.steps['ps_ct']]
	  
	print '\nFraction of restored perturbations is %5.3g (%d measurements)' % (p,correct_ct)
      elif sim.mode in ['eft_tmp']:
	p = 0.5
	correct_ct = 1
      #somewhere here read, whether p has already been sampled (before scaling down) and if so - read in and skip this step!
      if sim.mode in ['eft','eft_tmp']:
	#print idx_pert
	if p < 0.5 and not idx_pert:
	  print 'Scaling down the perturbation sizes by two orders to sample a lower region -> fix saving in such a case!'
	  self.pert_size *= 10**(-2)
	  scaled_down += 1
	  idx_pert -= 1
	  if scaled_down >=3:
	    print 'no stable trajectories can be found in the network, even after going down 6 orders of magnitude! - Interrupt!'
	    break
	else:
	  if ((p > 0.9) or (p < 0.1)):	# fast forward if far from 1-1/e
	    idx_pert+=3
	
	  if ((idx_pert+1 >= sim.steps['ps']) or (p<0.05) or (not correct_ct)):# or (self.pert_size[idx_pert] > 1)):
	    datinit_rd.close()
	    # create this to show read program, that sims are finished
	    ncid = netcdf.netcdf_file(self.Path['results_links'],'w')
	    print self.Path['results_links']
	    ncid.close()
	    break
	if self.pert_size[idx_pert+1] in pert_done:
	  print 'perturbation %g already done, skipping to next' % self.pert_size[idx_pert+1]
	  idx_pert += 1
	
      idx_pert += 1
      
	
      print 'time elapsed: %g sec\n' % (time.time()-time_start)
      
    # dummy comment (add even more)
    if sim.mode == 'eft' and not sim.cluster['q']:
      sim.steps['ps'] = sim.steps['ps_ct'] + 1
      self.fluxSimSv(sim.steps)
      #for dat in os.listdir(self.Path['results']):	# remove all data but the processed one
	#if not 'processed' in dat:
	  #os.remove(self.Path['results'] + dat)
    
    
    #------------------------------ support code -----------------------------
  
  def get_puppetTrain(self,spikes,times,max_time):
    
    puppet_tmp = [times[np.where(spikes==n)[0]] for n in range(self.Const['N'],self.ParaNet['N'][0])]

    # get maximum number of spikes and convert from absolute times to time intervalls
    self.puppetTrainSz = 0
    for n in range(self.Const['N']):
      len_tmp = len(puppet_tmp[n])
      for t in range(len_tmp-1,0,-1):
	puppet_tmp[n][t] -= puppet_tmp[n][t-1]
      if len_tmp > self.puppetTrainSz:
	self.puppetTrainSz = len_tmp
    
    self.puppetTrainSz += 1		# need one length more, to avoid core dumping
    # now write data to array
    self.ParaNet['puppetTrain'] = np.ones((self.Const['N'],self.puppetTrainSz)) + max_time
    for n in range(self.Const['N']):
      self.ParaNet['puppetTrain'][n][:len(puppet_tmp[n])] = puppet_tmp[n]
    self.ParaNet['puppetTrain'] *= 1000
  
  
  
  def create_puppetTrain(self,time_limit):
    #if self.topoConst['drive_type'] == 'poisson':
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
      
  
  
  def initDataSv(self,sim):
    self.DataSv = {}
    self.DataSv['pert_size'] = np.zeros(sim.steps['ps'])
    self.DataSv['p'] = np.zeros(sim.steps['ps'])
    self.DataSv['correct_ct'] = np.zeros(sim.steps['ps'])
    self.DataSv['CDBTimes'] = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    self.DataSv['div_track'] =  np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    self.DataSv['distances'] = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    
    for key in self.DataSv.keys():
      self.DataSv[key][:] = np.nan
  
  
  def get_pert_range(self,sim):
    self.eft_exp()
    
    if sim.steps['ps'] == 1:	# perturbation of expected fluxtube size (should look at others as well?
      if 'special' in self.Const.keys():
	if self.Const['special'] in [0]:
	  print "eft_exp: %g" % self.eft_exp
	  self.pert_size = np.array([self.eft_exp/100.])
	elif self.Const['special'] in [1,6,7,8]:
	  self.pert_size = [self.perturbation*2]
	  print self.perturbation
	elif self.Const['special'] in [2,3,4,5,30]:
	  self.pert_size = np.array([self.eft_exp])
	    
      else:
	self.pert_size = np.array([self.perturbation])
      
    else:
      if sim.mode == 'samples':
	if 'special' in self.Const.keys():
	  if self.Const['special'] in [0]:
	    print "pert: %g" % self.perturbation
	    print "eft_exp: %g" % self.eft_exp
	    self.pert_size = np.array([self.eft_exp/10.])
	    #self.pert_size = self.perturbation*10**(np.linspace(-1,1,sim.steps['ps']))
	  if self.Const['special'] == 2:
	    self.pert_size = np.nan
	  elif self.Const['special'] in [8]:
	    self.pert_size = self.perturbation*10**(np.linspace(-1,1,sim.steps['ps']))
	  
      else:
	if self.Const['N'] <= 1000:
	  self.pert_size = 10**np.linspace(-3,0,sim.steps['ps'])
	elif self.Const['N'] > 1000:
	  self.pert_size = 10**np.linspace(-4,0,sim.steps['ps'])
	#self.pert_size = self.eft_exp*10**np.linspace(-2,3,sim.steps['ps'])	# pertSizes logscaled around eft_exp
  
  
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
      
      ncidPuppet = netcdf.netcdf_file(self.Path['Puppet'],'r',mmap=False)
      puppetTrain = ncidPuppet.variables['puppetTrain'][:]
      self.puppetTrainSz = ncidPuppet.dimensions['puppetTrainSz']
      ncidPuppet.close()
      
      sigma = 0.1
      puppetTrain_tmp = puppetTrain + np.random.normal(scale=sigma,size=(self.Const['N'],self.puppetTrainSz))
      np.copyto(puppetTrain_tmp,puppetTrain,where=puppetTrain_tmp<0)
      
      self.ParaNet['puppetTrain'] = puppetTrain_tmp
      self.ParaNet['init'] = init
      return 0
    
    elif 'samples' in self.scriptName and 'special' in self.Const.keys():
      if self.Const['special'] in [2,3,4,5,30]:
	trash = self.ParaNet.pop('init')
	return 0
      else:
	for i in range(10):	# loop to ensure not perturbing over threshold
	  pert_vect = OrthoVector(np.random.uniform(-1,1,self.Const['N']))	# get orthogonal vector
	  pert_vect *= pert_size/np.linalg.norm(pert_vect)				# ... normalize and scale it ...
	  if self.topo == 'p':	# when poisson neurons present
	    pert_vect = np.append(pert_vect,np.zeros(self.Const['N']))		# (... extend by 0s for poisson neurons ...)
	  init_pert = init + pert_vect						# ...and add it to initial phase point
	  if max(init_pert[:self.Const['N']]) < 1:
	    self.ParaNet['init'] = init_pert
	    return 0
    else:
      for i in range(10):	# loop to ensure not perturbing over threshold
	pert_vect = OrthoVector(np.random.uniform(-1,1,self.Const['N']))	# get orthogonal vector
	#pert_vect = np.random.uniform(-1,1,self.Const['N'])	# get non-orthogonal vector (still orthogonal)
	
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
  
  
  
  
  
  def eft_exp(self):
    self.eft_exp = 1./(math.sqrt(self.Const['N']*self.Const['K'])*self.Const['rateWnt']*self.Const['tauM']/1000.)		# expected tube size from paper as medium guess
  
  
  def create_topology(self,hide):
      
    N = self.Const['N']
    K = self.Const['K']
    if K/float(N) < 0.5:
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
	if (self.Const['K'] >= self.Const['N']/2.) or self.Const['N'] >= 2*10**5:
	  i = 100
	else:
	  runIt(None,run_str)		# get in- and out-degree from here if needed for sims
	  
	  try:
	    # read in the results of topology creation
	    ncid = netcdf.netcdf_file(str_save_tmp,'r',mmap=False)
	    NcItems = np.copy(ncid.variables.items())
	    for item in NcItems:
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
    else:
      graph = random_graph(N,K)
      for item in graph.items():
	self.ParaTopo[item[0]] = item[1]
      found = 1
	
    if (len(np.unique(np.array(self.ParaNet['J0']))) > 1) or (self.topo == 'p'):
      self.ParaTopo['J'] = np.ones(len(self.ParaTopo['postSyn']))*self.ParaNet['J0']/np.sqrt(K)
	
      if self.topo == 'p':
	if self.topoConst['corr'] > 1:
	  n_poisson = N/self.topoConst['corr']
	  assert not (N%self.topoConst['corr']), 'make sure N is a multiple of the number of synpases per poisson neuron!'
	  #print 'Number of poisson neurons added: ', (n_poisson)
	  self.ParaTopo['J'] = np.append(self.ParaTopo['J'],np.ones(N)*self.topoConst['drive_cplg']/np.sqrt(self.Const['K']))
	  preSyn_tmp = np.zeros(len(self.ParaTopo['J']))
	  
	  post_syn_draw = range(N)
	  
	  poisson_wiring = np.zeros((n_poisson,self.topoConst['corr']))
	  
	  for n in range(n_poisson):
	    for i in range(self.topoConst['corr']):
	      poisson_wiring[n,i] = post_syn_draw.pop(np.random.randint(len(post_syn_draw)))
	    
	    self.ParaTopo['postSyn'] = np.append(self.ParaTopo['postSyn'],np.unique(poisson_wiring[n,:]))
	  
	  idx_add = 0
	  idx_deg = 0
	  
	  for inDeg in self.ParaTopo['inDegree']:
	    preSyn_tmp[idx_deg+idx_add:idx_deg+idx_add+inDeg] = self.ParaTopo['preSyn'][idx_deg:idx_deg+inDeg]
	    preSyn_tmp[idx_deg+idx_add+inDeg] = np.where(poisson_wiring == idx_add)[0][0]+self.Const['N']
	    #print preSyn_tmp[idx_deg+inDeg+idx_add]
	    idx_deg += inDeg
	    idx_add += 1
	  self.ParaTopo['preSyn'] = preSyn_tmp
	  
	  self.ParaTopo['inDegree'] = np.append(self.ParaTopo['inDegree']+1,np.zeros(n_poisson))
	  self.ParaTopo['outDegree'] = np.append(self.ParaTopo['outDegree'],np.ones(n_poisson)*self.topoConst['corr'])
	      
	else:
	  self.ParaTopo['J'] = np.append(self.ParaTopo['J'],np.ones(N)*self.topoConst['drive_cplg']/np.sqrt(self.Const['K']))
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

  
  # ---------------------------------------------- data structure stuff ----------------------------
  
  def data_structure(self,sim,hide=1,call=0):
    
    sv_text = 0
    
    #generate output data paths
    self.Path = {}
    
    if sim.cluster['scratch']:
      self.Path['session'] = '/scratch%02d/aschmidt/data/' % sim.cluster['scratch']	#path for storage on scratch
    else:
      self.Path['session'] = 'data/'			#path for local storage
    if (not os.path.exists(self.Path['session']) and not call):
      os.mkdir(self.Path['session'])
    
    self.Path['script'] = self.Path['session'] + self.scriptName + '/'
    if (not os.path.exists(self.Path['script']) and not call):
      print 'creating new directory for this script: %s' % self.Path['script']
      os.mkdir(self.Path['script'])
    
    self.Path['info'] = self.Path['session'] + 'info/' + self.scriptName + '/'
    
    if (not os.path.exists(self.Path['session'] + 'info') and not call):
      os.mkdir(self.Path['session'] + 'info')
    if (not os.path.exists(self.Path['info']) and not call):
      os.mkdir(self.Path['info'])
    
    self.InfoHash(sim,call)
    
    self.Path['inpara'] = self.Path['script'] + 'inparas_' + self.Hash + '/'
    self.Path['results'] = self.Path['script'] + 'results_' + self.Hash + '/'
    if self.break_it:
      return
	
    
    if (not os.path.exists(self.Path['inpara']) and not call):
      os.mkdir(self.Path['inpara'])
    
    self.Path['inTopo_links'] = self.Path['inpara'] + 'datTopo.nc'
    self.Path['inpara_links'] = self.Path['inpara'] + 'datInit.nc'
    
    if (not os.path.exists(self.Path['results']) and not call):
      os.mkdir(self.Path['results'])
    self.Path['results_links'] = self.Path['results'] + 'filenames.nc'
    
    if not hide:
      self.Path['txt'] = self.Path['script'] + 'txt_' + self.Hash + '/'
      if (not os.path.exists(self.Path['txt']) and not call):
	print 'creating %s' % self.Path['txt']	
	os.mkdir(self.Path['txt'])
      	
    else:
      self.Path['txt'] = '/dev/null'


  def InfoHash(self,sim,call=0):
    # creating name for identification of simulation in info file
    self.Hash_array = []
    self.break_it = 0
    #print self.Path
    exists = 0

    if not call:
      para_str_name, para_str_write = self.Info_generate_lines(self.Const,'float','para')
      topo_str_name, topo_str_write = self.Info_generate_lines(self.topoConst,'float','topo')
      step_str_name, step_str_write = self.Info_generate_lines(sim.steps,'int','step')
      saves_str_name, saves_str_write = self.Info_generate_lines(sim.Paras,None,'save')
      self.Path['info_file'] = self.Path['info'] + 'SimInfo_%s%s%s' % (self.scriptName,para_str_name,topo_str_name)
      self.infoName = self.Path['info_file']
      test_str = 'SimInfo_%s%s%s' % (self.scriptName,para_str_name,topo_str_name)
      if self.topo == 'S':
	self.Path['info_file'] += '_'
	test_str += '_'
      self.Path['info_file'] += '.txt'
      for dat in os.listdir(self.Path['info']):
	if test_str in dat:
	  exists = 1
	  break
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
    
    if len(self.DataSv['p'].shape)==1:
      for item in self.DataSv.items():
	self.DataSv[item[0]] = item[1][:steps['ps']]
    #print steps['ps']
    
    processedPath = self.Path['results'] + 'results_processed.nc'
    
    if os.path.exists(processedPath):
      os.remove(processedPath)
    
    ncid = netcdf.netcdf_file(processedPath,'w')
    #print processedPath
    ncid.createDimension('one',1)
    ncid.createDimension('ps_steps',steps['ps'])
    ncid.createDimension('co_steps',steps['co'])
    ncid.createDimension('in_steps',steps['in'])
    ncid.createDimension('pd_steps',steps['pd'])
    
    VarN = ncid.createVariable('N','i', ('one',))
    VarK = ncid.createVariable('K','i', ('one',))
    Varf = ncid.createVariable('f','i', ('one',))
    VartauM = ncid.createVariable('tauM','i', ('one',))
    VarPertSz = ncid.createVariable('pert_size', 'd', ('ps_steps',))
    VarCorr = ncid.createVariable('correct_ct','i',('ps_steps',))

    if len(self.DataSv['p'].shape)==1:
      VarProbs = ncid.createVariable('p','d', ('ps_steps',))
      VarDivTrack = ncid.createVariable('div_track','d',('ps_steps','co_steps','in_steps','pd_steps'))
      VarTimes = ncid.createVariable('CDBTimes','d',('ps_steps','co_steps','in_steps','pd_steps'))
      VarDist = ncid.createVariable('distances','d',('ps_steps','co_steps','in_steps','pd_steps'))
      
      VarTimes[:] = self.DataSv['CDBTimes']
      
    else:
      ncid.createDimension('t_steps',self.DataSv['p'].shape[0])
      
      VarProbs = ncid.createVariable('p','d', ('t_steps','ps_steps'))
      VarDivTrack = ncid.createVariable('div_track','d',('t_steps','ps_steps','co_steps','in_steps','pd_steps'))
      VarDist = ncid.createVariable('distances','d',('t_steps','ps_steps','co_steps','in_steps','pd_steps'))
    
    VarN[:] = self.Const['N']
    VarK[:] = self.Const['K']
    Varf[:] = self.rateWnt
    VartauM[:] = self.Const['tauM']
    
    VarPertSz[:] = self.DataSv['pert_size']
    #print self.DataSv['pert_size']
    VarProbs[:] = self.DataSv['p']
    print self.DataSv['p']
    VarDivTrack[:] = self.DataSv['div_track']
    
    VarDist[:] = self.DataSv['distances']
    VarCorr[:] = self.DataSv['correct_ct']
    
    ncid.close()
  
  
  def clean(self):
    
    if not os.path.exists(self.Path['results']+'/clean.nc'):
      print 'Cleaning up reference trajectories'
      for dat in os.listdir(self.Path['results']):
	if not ('processed' in dat) and not ('filenames' in dat) and not ('test' in dat) and not ('clean.' in dat):
	  try:
	    os.remove(self.Path['results'] + dat)
	  except:
	    shutil.rmtree(self.Path['results'] + dat,ignore_errors=True)	    
      
      print 'Cleaning up input data'
      for dat in os.listdir(self.Path['inpara']):
	if not ('datTopo' in dat) and not ('ParaTop' in dat):
	  try:
	    os.remove(self.Path['inpara'] + dat)
	  except:
	    shutil.rmtree(self.Path['inpara'] + dat,ignore_errors=True)	    

      print 'Cleaning up text files'
      if os.path.exists(self.Path['txt']):
	shutil.rmtree(self.Path['txt'],ignore_errors=True)
      
      # create proxy file to show that directory was cleaned up already
      ncid_proxy = netcdf.netcdf_file(self.Path['results']+'/clean.nc','w')
      ncid_proxy.close()
    else:
      print 'Data is all cleaned up already!'
  
  
  def inter_read(self,sim,read_max=100,errorMsg=0):
    
    ncid = netcdf.netcdf_file(self.Path['results_links'] + '_' + str(sim.steps['ps_ct']),'r',mmap=False)
    correct_ct = 0
    div_track = np.zeros((sim.steps['co'],sim.steps['in'],sim.steps['pd']))
    div_track[:] = np.nan
    
    for co_ct in range(sim.steps['co']):
	
      for in_ct in range(sim.steps['in']):
	link = ''.join(ncid.variables['links'][co_ct,in_ct,-1])
	try:
	  assert os.path.exists(link)
	  ncid_read = netcdf.netcdf_file(link,'r',mmap=False)
	  distances = ncid_read.variables['PTdistances'][:][:sim.steps['pd']/2]
	  ncid_read.close()
	  mask = np.invert(np.isnan(distances)) # mask entries, where no perturbation could be found
	  correct_ct += np.sum(mask)
	  div_track[co_ct][in_ct][mask] = (distances >= 0.1)[mask].astype('int')
	except:
	  1
	  #ncid.close()
	
	if (correct_ct >= read_max) or ((in_ct==sim.steps['in']-1) and (co_ct==sim.steps['co']-1)):	# break, when end or maximum number of allowed readout-files is reached
	  if correct_ct:
	    p = 1 - np.nansum(div_track)/correct_ct
	  else:
	    p = np.nan
	  ncid.close()
	  
	  return p, correct_ct
	    
  
  def prob_read(self,sim,steps_rd=None,errorMsg=0):
    
    sv_all = 0
    
    #ps_steps = len(steps_rd)
    #pert_size = np.zeros(sim.steps['ps'])
    if not steps_rd:			# if not specified, read all
      steps_rd = range(sim.steps['ps'])
    
    if sim.mode in ['eft','eft_tmp']:

      if sim.mode == 'eft':
	div_track = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	dist = np.ones((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	times = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	div_track[:] = np.nan
      
      print self.DataSv.keys()
      if sim.mode == 'eft_tmp':
	self.DataSv['div_track'] = np.zeros((21,sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	self.DataSv['distances'] = np.zeros((21,sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	self.DataSv['p'] = np.zeros((21,sim.steps['ps']))
	trash = self.DataSv.pop('CDBTimes')
	#self.DataSv['correct_ct'] = np.zeros(sim.steps['ps'])
	
	dist = np.ones((21,sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	times = np.zeros((sim.steps['ps'],sim.steps['co'],sim.steps['in'],sim.steps['pd']))
	self.DataSv['div_track'][:] = np.nan
      
      
    
    print 'reading %d files from list %s...' % (sim.steps['ps']*sim.steps['total'],self.Path['results_links'])
    
    idx_ps = 0
    idx_ct = 0
    
    #if sim.mode == 'samples':
    idx_rec_cut = int(sim.TC*self.Const['N']*self.Const['rateWnt'])	# number of spikes in the recurrent neurons after TC secs
    #print 'idx rec cut: ', idx_rec_cut
    if sim.mode == 'samples':
      
      #self.analysis['measureTimes'] = sim.Paras['measureTimes']
      if 'measureTimes' in sim.Paras.keys():
	time_len = len(sim.Paras['measureTimes'])
	self.analysis['measureTimes'] = sim.Paras['measureTimes']
      elif 'measureSpikes' in sim.Paras.keys():
	time_len = len(sim.Paras['measureSpikes'])
	#self.analysis['measureTimes'] = sim.Paras['measureSpikes']
      print self.analysis.keys()
      #time_len /= 10
      
      if self.Const['special'] in [0,1,6,7,8]:		# few sims for plotting
	
	if self.Const['special'] == 0:		# large batch for obtaining conv/div statistics
	  self.analysis['convergeVelocity_1st'] = np.zeros((2,sim.steps['total']*sim.steps['ps']))
	  self.analysis['convergeVelocity_1st'][:] = np.nan
	  self.analysis['convergeVelocity_2nd'] = np.zeros((2,sim.steps['total']*sim.steps['ps']))
	  self.analysis['convergeVelocity_2nd'][:] = np.nan
	    
	  if self.topo == 'p' or self.topoConst['alpha_conv'] > 0:
	    self.analysis['kinkTime'] = np.zeros(sim.steps['total']*sim.steps['ps'])
	    self.analysis['kinkTime'][:] = np.nan
	      
	if self.Const['special'] == 1:
	  
	  self.analysis['distance'] = np.zeros((sim.steps['total']*sim.steps['ps'],time_len))
	  self.analysis['distanceVector'] = np.zeros((sim.steps['total'],time_len,self.Const['N']))
	  self.analysis['distanceVector'][:] = np.nan
	  #self.analysis['distanceOrthVector'] = np.zeros((sim.steps['total'],time_len,self.Const['N']))
	  #self.analysis['distanceOrthVector'][:] = np.nan
	  
	  #self.analysis['train'] = np.zeros((sim.steps['co']*sim.steps['in'],idx_rec_cut,2))
	  #self.analysis['train'][:] = np.nan
	  
	  #self.analysis['trainPert'] = np.zeros((sim.steps['total']*sim.steps['ps'],idx_rec_cut,2))
	  #self.analysis['trainDeltaTime'] = np.zeros((sim.steps['total']*sim.steps['ps'],idx_rec_cut))
	  
	  #self.analysis['trainPert'][:] = np.nan
	  #self.analysis['trainDeltaTime'][:] = np.nan
	  
	  self.analysis['outDegree'] = np.zeros((sim.steps['total']*sim.steps['ps'],self.Const['N']))
	  self.analysis['inDegree'] = np.zeros((sim.steps['total']*sim.steps['ps'],self.Const['N']))

	  #if self.topo == 'p':
	    #idx_drive_cut = int(sim.TC*self.Const['N']*self.topoConst['drive_rate'])

	    #self.analysis['trainDrive'] = np.zeros((sim.steps['total']*sim.steps['ps'],idx_drive_cut,2))
	    #self.analysis['trainDrive'][:] = np.nan
	elif self.Const['special'] == 7:
	  assign = {}
	  assign['train'] = np.zeros((sim.steps['pd']+1,idx_rec_cut,2))
	  assign['train'][:] = np.nan
	  
	elif self.Const['special'] in [8]:
	  assign = {}
	  assign['train'] = np.zeros((sim.steps['pd']+1,idx_rec_cut,2))
	  assign['train'][:] = np.nan
	  print 'shape: ', 	assign['train'].shape
	  self.analysis['spikecross'] = np.zeros((sim.steps['total']*sim.steps['ps'],2))
	  self.analysis['spikecross'][:] = np.nan
	
	distanceOrthVector = np.zeros((time_len,self.Const['N']))
	self.analysis['distanceOrth'] = np.zeros((sim.steps['total']*sim.steps['ps'],time_len))
      
      
      
      elif self.Const['special'] in [2,3,4,5,30]:
	all_to_all = sim.steps['co']*sim.steps['in']*np.sum(np.arange(sim.steps['pd']))
	
	if self.Const['special'] == 2:
	  # compress into one!
	  self.analysis['D_decorr'] = np.zeros((2,all_to_all))
	  self.analysis['D_decorr'][:] = np.nan
	  statesTmp = np.zeros((sim.steps['pd'],time_len,self.Const['N']))
	  distanceOrthVector = np.zeros((time_len,self.Const['N']))
	  #hist_bin = 50
		      
	if self.Const['special'] in [3,4,5,30]:
	  self.analysis['measureStates'] = np.zeros((sim.steps['co'],sim.steps['pd'],time_len,self.Const['N']))
	  
	  thresh_drop = 10**(-3)
	  if self.Const['special'] in [3,30]:
	    self.analysis['T_conv'] = np.zeros((sim.steps['co'],sim.steps['in'],sim.steps['pd'],sim.steps['pd']))
	    self.analysis['T_conv'][:] = np.nan
	    
	    self.analysis['T_conv_measure'] = np.zeros((sim.steps['co'],sim.steps['in'],sim.steps['pd'],sim.steps['pd']))
	    self.analysis['T_conv_measure'][:] = np.nan
	    
	    # get offset due to convergence time
	    
	    if 'D_decorr' in self.analysis.keys():
	      if self.analysis['D_decorr'][0] > 10 or np.isnan(self.analysis['D_decorr'][0]):
		self.analysis['D_decorr'][0] = 1
		print "D decorr is set to 1 (no D found in simulation)"
	    else:
	      self.analysis['D_decorr'] = [1,0,0]
	      print "D decorr is set to 1 (no simulation found for D)"
	    # use second time constant, as network settled to shift already
	    dt = np.log(thresh_drop/self.analysis['D_decorr'][0])*(-0.01)
	    print 'offset due to finite convergence speed: %g' % dt

	    
	  elif self.Const['special'] == 4:
	    #self.analysis['train'] = np.zeros((sim.steps['total'],idx_rec_cut,2))
	    #self.analysis['train'][:] = np.nan
	    
	    distanceOrthVector = np.zeros((time_len,self.Const['N']))
	    self.analysis['distanceOrth'] = np.zeros((all_to_all,time_len))
	    
	    assign = {}
	    assign['train'] = np.zeros((sim.steps['pd'],idx_rec_cut,2))
	    assign['train'][:] = np.nan
	    
	    self.analysis['trainDeltaT'] = np.zeros((all_to_all,idx_rec_cut,3))
	    self.analysis['trainDeltaT'][:] = np.nan
	    
	    self.analysis['phi_corr'] = np.zeros((all_to_all,self.Const['N']))
	    slide_len = 10
	    num_thresh = 10
	    frac_conv = np.zeros((time_len,num_thresh))
	    self.analysis['threshold_peaks'] = np.zeros((all_to_all,num_thresh))
	    
	    
	  elif self.Const['special'] == 5:
	    var = np.zeros((sim.steps['pd'],self.Const['N']))
	## include error measure?!
      # prepare distance measurements
      
      
		      
      
    #if not os.path.exists(self.Path['results_links']):
      #print "Simulations not yet finished"
      #return
    
    for ps_ct in steps_rd:
      correct_ct = 0
      
      try:
	ncid = netcdf.netcdf_file(self.Path['results_links'] + '_%d' % ps_ct,'r',mmap=False)
      except:
	sim.steps['ps'] = ps_ct
	break
      self.DataSv['pert_size'][ps_ct] = ncid.variables['pert_size'].getValue()

      for co_ct in range(sim.steps['co']):
	#try:
	if sim.mode in ['eft','eft_tmp']:
	#if self.Const['special'] in [-1,0,1]:
	  if self.topo == 'p':
	    link_topo = ''.join(ncid.variables['links'][co_ct,0,-3])
	  else:
	    link_topo = ''.join(ncid.variables['links'][co_ct,0,-2])
	elif self.Const['special'] in [2,3,4,30]:
	  link_topo = ''.join(ncid.variables['links'][co_ct,0,0,-3])
	else:
	  link_topo = ''.join(ncid.variables['links'][co_ct,0,-2])
	
	ncidTop = netcdf.netcdf_file(link_topo,'r',mmap=False)
	Topo = {}
	Topo['adjMat'] = postsynapticvector2adjacencymatrix(ncidTop.variables['postSyn'][:],ncidTop.variables['outDegree'][:])
	Topo['inDegree'] = ncidTop.variables['inDegree'][:]
	Topo['outDegree'] = ncidTop.variables['outDegree'][:]
	Topo['postSyn'] = ncidTop.variables['postSyn'][:]
	
	
	#print ncidTop.variables.keys()
	ncidTop.close()
	inDeg_min_idx = np.in1d(Topo['inDegree'],heapq.nsmallest(int(0.1*self.Const['N']),Topo['inDegree']))
	inDeg_max_idx = np.in1d(Topo['inDegree'],heapq.nlargest(int(0.1*self.Const['N']),Topo['inDegree']))
	#except:
	  #if errorMsg:
	    #print linkName + ' does not exist. Something seems to have gone wrong'

	#print heapq.nsmallest(int(0.01*self.Const['N']),Topo['inDegree'])
	#print heapq.nlargest(int(0.01*self.Const['N']),Topo['inDegree'])
	
	#print inDeg_min_idx
	#print inDeg_max_idx
	#print len(Topo['inDegree'][inDeg_min_idx])
	#print len(Topo['inDegree'][inDeg_max_idx])
	#print Topo
	
	
	for in_ct in range(sim.steps['in']):
	  
	  if sim.mode in ['eft','eft_tmp']:
	    try:
	      linkTopo = ''.join(ncid.variables['links'][co_ct,in_ct,-1])
	      assert os.path.exists(linkTopo)
	      ncid_read = netcdf.netcdf_file(linkTopo,'r',mmap=False)
	      
	      if sim.mode == 'eft':
		distances = ncid_read.variables['PTdistances'][:]
		mask = np.invert(np.isnan(distances)) # mask entries, where no perturbation could be found
		self.DataSv['distances'][ps_ct,co_ct,in_ct] = distances
		self.DataSv['div_track'][ps_ct,co_ct,in_ct][mask] = (distances >= 0.1)[mask].astype('int')
		self.DataSv['CDBTimes'][idx_ps][co_ct][in_ct] = ncid_read.variables['PTtimes'][:]
	      else:
		distances = ncid_read.variables['distances'][:]
		mask = np.invert(np.isnan(distances[:,0])) # mask entries, where no perturbation could be found
		self.DataSv['distances'][:,ps_ct,co_ct,in_ct] = distances.transpose()
		
		for p in range(ncid_read.dimensions['PTSz']):
		  if mask[p]:
		    for t in range(ncid_read.dimensions['measureTimesSz']):
		      
		      if t < np.sum(mask):
			self.DataSv['div_track'][t,ps_ct,co_ct,in_ct,p] = (distances[p,t] >= 0.1).astype('int')
		      else:
			self.DataSv['div_track'][t,ps_ct,co_ct,in_ct,p] = self.DataSv['div_track'][t-1,ps_ct,co_ct,in_ct,p]
	      ncid_read.close()
	      
	      correct_ct += np.sum(mask)
	    except:
	      1
	      #print "this batch went wrong..."
	  elif sim.mode == 'samples':
	    
	    if self.Const['special'] in [0,1,6,7,8]:
	      #try:
		str_print = "read simulation %d/%d  " % (in_ct+sim.steps['in']*co_ct,sim.steps['co']*sim.steps['in'])
		sys.stdout.write('\r{0}'.format(str_print))
		sys.stdout.flush()
		try:
		#print ncid.variables['links'][co_ct,in_ct,-1]
		  link = ''.join(ncid.variables['links'][co_ct,in_ct,-1])
		except:
		  link = ''.join(ncid.variables['links'][ps_ct,co_ct,in_ct,-1])
		assert os.path.exists(link), '%s does not exist'
		
		Data = readDataOut(link)
		
		if self.Const['special'] in [7,8]:
		  idx_rec_cut_tmp = min(idx_rec_cut,np.sum(Data['train'][:,0] < self.Const['N']),min([np.sum(x < self.Const['N']) for x in Data['trainPert'][:,:,0]]))
		  assign['train'][0,:idx_rec_cut_tmp] = Data['train'][:idx_rec_cut_tmp]
		  assign['train'][1:,:idx_rec_cut_tmp] = Data['trainPert'][:,:idx_rec_cut_tmp]
		
		#if self.Const['special'] in [1]:
		  
		  #idx_ref = in_ct + sim.steps['in']*co_ct
		  #idx_rec = Data['train'][:,0] < self.Const['N']
		  #idx_rec_cut_tmp = min(idx_rec_cut,np.sum(idx_rec))
		  
		  #self.analysis['train'][idx_ref][:idx_rec_cut_tmp] = Data['train'][idx_rec][:idx_rec_cut_tmp]
		  
		  #if self.topo == 'p':
		    #idx_drive = np.invert(idx_rec)
		    #idx_drive_cut_tmp = min(idx_drive_cut,np.sum(idx_drive))
		    
		    #self.analysis['trainDrive'][idx_ref][:idx_drive_cut_tmp] = Data['train'][idx_drive][:idx_drive_cut_tmp]
		
		for pd_ct in range(sim.steps['pd']):
		  str_print = "\t\t\t, process simulation %d/%d" % (idx_ct,sim.steps['total']*sim.steps['ps'])
		  sys.stdout.write('\r{0}'.format(str_print))
		  sys.stdout.flush()
		  
		  #if self.Const['special'] == 0:
		  try:
		    time_len = np.where(Data['PTtimes'][pd_ct] < Data['measureTimes'])[0][0]
		  except:
		    1
		  #distanceVector[:] = np.nan
		  
		  distanceVector = (Data['measureStates']-Data['measureStatesPert'][pd_ct])[:time_len,:self.Const['N']]
		  #distanceOrthVector[:] = np.nan
		  for t in range(time_len):
		    distanceOrthVector[t] = OrthoVector(distanceVector[t])
		  
		  self.analysis['distanceOrth'][idx_ct][:time_len] = np.sqrt(np.sum(distanceOrthVector[:time_len]**2,axis=1))
		  self.analysis['distanceOrth'][idx_ct][time_len:] = np.nan
		  
		  if self.Const['special'] == 0:

		    mask = np.invert(np.isnan(self.analysis['distanceOrth'][idx_ct]))
		    distance_med = self.analysis['distanceOrth'][idx_ct][mask]
		    #distance_med = sp.signal.medfilt(self.analysis['distanceOrth'][idx_ct][mask],1)
		    #self.analysis['distanceOrth'][idx_ct][mask] = distance_med

		    if distance_med[-1] < 0.01 and not any(distance_med > 1):# only get convergence speed for converged runs
		      
		      try:
			drop_idx = np.argwhere(distance_med<10**(-10))[0][0]
		      except:
			drop_idx = -1
		      
		      slope_guess = 1./self.analysis['measureTimes'][drop_idx]*np.log(distance_med[drop_idx]/distance_med[0])
		      intercept_guess = distance_med[0]
		      
		      if self.topo == 'p':
			popt, pcov = curve_fit(func_exp_intersect,self.analysis['measureTimes'][:drop_idx],distance_med[:drop_idx],sigma=distance_med[:drop_idx],p0=[2*slope_guess,slope_guess,intercept_guess,0.01*intercept_guess])
			perr= np.sqrt(np.diag(pcov))
			
			intersect = np.log(popt[3]/popt[2])/(popt[0]-popt[1])
			#print "init / latter: ", popt[0], popt[1], " intersect @ ", intersect
			
			idx_intersect = np.argmin(abs(self.analysis['measureTimes'][mask]-intersect))
			
			if intersect > 0.005:
			  self.analysis['convergeVelocity_1st'][0,idx_ct] = 1./popt[0]
			  self.analysis['convergeVelocity_1st'][1,idx_ct] = perr[0]/popt[0]**2
			  self.analysis['convergeVelocity_2nd'][0,idx_ct] = 1./popt[1]
			  self.analysis['convergeVelocity_2nd'][1,idx_ct] = perr[1]/popt[1]**2
			  self.analysis['kinkTime'][idx_ct] = intersect
			  
			else:
			  popt, pcov = curve_fit(func_exp,self.analysis['measureTimes'][:drop_idx],distance_med[:drop_idx],sigma=distance_med[:drop_idx],p0=[slope_guess,intercept_guess])
			  
			  perr= np.sqrt(np.diag(pcov))
			  
			  self.analysis['convergeVelocity_1st'][0,idx_ct] = 1./popt[0]
			  self.analysis['convergeVelocity_1st'][1,idx_ct] = perr[0]/popt[0]**2
			
			  self.analysis['convergeVelocity_2nd'][0,idx_ct] = 1./popt[0]
			  self.analysis['convergeVelocity_2nd'][1,idx_ct] = perr[0]/popt[0]**2
			  
		      else:
			popt, pcov = curve_fit(func_exp,self.analysis['measureTimes'][:drop_idx],distance_med[:drop_idx],p0=[slope_guess,intercept_guess],sigma=distance_med[:drop_idx])
			perr= np.sqrt(np.diag(pcov))
			
			self.analysis['convergeVelocity_2nd'][0,idx_ct] = 1./popt[0]
			self.analysis['convergeVelocity_2nd'][1,idx_ct] = perr[0]/popt[0]**2
			
		      if np.isinf(self.analysis['convergeVelocity_1st'][1,idx_ct]):
			self.analysis['convergeVelocity_1st'][:,idx_ct] = np.nan
		      
		      if np.isinf(self.analysis['convergeVelocity_2nd'][1,idx_ct]):
			self.analysis['convergeVelocity_2nd'][:,idx_ct] = np.nan
		      
		      #intersect = np.log(popt[3]/popt[2])/(popt[0]-popt[1])
		      
		    else:
		      self.analysis['convergeVelocity_1st'][:,idx_ct] = np.nan
		      self.analysis['convergeVelocity_2nd'][:,idx_ct] = np.nan
		  
		  elif self.Const['special'] == 1:
		    self.analysis['distance'][idx_ct][:time_len] = np.sqrt(np.sum(distanceVector[:time_len]**2,axis=1))
		    self.analysis['inDegree'][idx_ct] = Topo['inDegree'][:self.Const['N']]
		    self.analysis['outDegree'][idx_ct] = Topo['outDegree'][:self.Const['N']]
		  
		    self.analysis['distanceVector'][idx_ct] = distanceVector
		    #self.analysis['distanceOrthVector'][idx_ct] = distanceOrthVector
		    
		    #idx_rec = Data['train'][pd_ct] < self.Const['N']
		    #idx_rec_cut_tmp = min(idx_rec_cut,np.sum(idx_rec))
		    #self.analysis['trainPert'][idx_ct][:idx_rec_cut_tmp] = Data['trainPert'][pd_ct][idx_rec][:idx_rec_cut_tmp]
		    
		  elif self.Const['special'] == 6:
		    plt.figure()
		    plt.plot(Data['measureTimes'],self.analysis['distanceOrth'][idx_ct])
		    plt.yscale('log')
		    plt.show()
		  
		  elif self.Const['special'] == 7:
		    if self.analysis['distanceOrth'][idx_ct][-1] > 0.01:
		      spike_cross = assignSpikes(assign,0,pd_ct+1,idx_rec_cut_tmp,Topo=Topo,find_decorr=True)
		      print "cross at:"
		      print spike_cross
		      print spike_cross[1,1]
		      #time.sleep(5)
		      idx_cross = np.where((self.analysis['measureTimes']-spike_cross[1,1])>0)[0][0] - 1
		      print idx_cross
		      idx_top = np.where(self.analysis['distanceOrth'][idx_ct] > 0.9*np.max(self.analysis['distanceOrth'][idx_ct]))[0][0]
		      
		      popt, pcov = curve_fit(func_exp,self.analysis['measureTimes'][idx_cross:idx_top],self.analysis['distanceOrth'][pd_ct][idx_cross:idx_top])
		      
		      perr = np.sqrt(np.diag(pcov))
		      
		      plt.figure()
		      plt.plot(Data['measureTimes'],self.analysis['distanceOrth'][idx_ct])
		      plt.plot(spike_cross[1,1],self.analysis['distanceOrth'][idx_ct][idx_cross],'ok')
		      plt.plot(Data['measureTimes'][idx_top],self.analysis['distanceOrth'][idx_ct][idx_top],'ok')
		      plt.plot(Data['measureTimes'],popt[1]*np.exp(popt[0]*Data['measureTimes']),'r')
		      #plt.xscale('log')
		      plt.yscale('log')
		      
		      plt.ylim([10**-3,10**1.5])
		      plt.show()
		  
		  
		  elif self.Const['special'] == 8:
		    #assign['train'][1][:idx_rec_cut_tmp] = Data['trainPert'][pd_ct][:idx_rec_cut_tmp]
			  
		    #print Topo.keys()
		    #print "assign:"
		    spike_cross = assignSpikes(assign,0,pd_ct,idx_rec_cut_tmp,Topo=Topo,find_decorr=True)
		    if len(spike_cross) == 2:
		      #print spike_cross
		      self.analysis['spikecross'][idx_ct] = Topo['outDegree'][spike_cross[:,0].astype('int')]
		      
		      #print spike_cross[:,0], Topo['outDegree'][spike_cross[:,0].astype('int')]
		      
		      #print self.analysis['distanceOrth'][idx_ct]
		      #print self.analysis['distanceOrth'].shape
		      #plt.figure()
		      
		      #plt.plot(self.analysis['measureTimes'],self.analysis['distanceOrth'][idx_ct])
		      #plt.plot([spike_cross[1,1],spike_cross[1,1]],[10**(-6),10**(1)])
		      #plt.yscale('log')
		      #plt.show()
		      
		  idx_ct += 1
	      #except:
		#1
	    elif self.Const['special'] in [2,3,4,5,30]:
	      # this should do all-to-all comparison
	      
	      for pd_ct in range(sim.steps['pd']):	#read in all stuff
		str_print = "read simulation %d/%d  " % (pd_ct + in_ct*sim.steps['pd'] + co_ct*sim.steps['in']*sim.steps['pd'],sim.steps['total'])
		sys.stdout.write('\r{0}'.format(str_print))
		sys.stdout.flush()
		
		link_data = ''.join(ncid.variables['links'][co_ct,in_ct,pd_ct,-2])
		#print link_data
		try:
		  Data = readDataOut(link_data)
		  #print 'reading...'
		  if not pd_ct:
		    self.analysis['measureTimes'] = Data['measureTimes']
		  if self.Const['special'] == 2:
		    statesTmp[pd_ct] = Data['measureStates'][:time_len,:self.Const['N']]
		    
		  if self.Const['special'] in [3,4,5,30]:
		    
		    self.analysis['measureStates'][co_ct][pd_ct] = Data['measureStates'][:time_len,:self.Const['N']]
		    #print Data['measureStates'][:time_len,:self.Const['N']]
		    if self.Const['special'] == 4:
		      idx_rec = Data['train'][:,0] < self.Const['N']
		      idx_rec_cut_tmp = min(idx_rec_cut,np.sum(idx_rec))
		      
		      idx_tmp = pd_ct + sim.steps['pd']*in_ct + sim.steps['pd']*sim.steps['in']*co_ct
		      #self.analysis['train'][idx_tmp][:idx_rec_cut_tmp] = Data['train'][idx_rec][:idx_rec_cut_tmp]
		      assert all(Data['train'][:,1]<sim.TC*1.1), 'nope'
		      assert all(Data['train'][:,1]>0), 'nope'
		      assign['train'][pd_ct,:idx_rec_cut_tmp] = Data['train'][idx_rec][:idx_rec_cut_tmp]

		    if self.Const['special'] == 5:
		      var[pd_ct] = np.var(self.analysis['measureStates'][co_ct][pd_ct],axis=0)
		    #print "var shape:" 
		    #print var[pd_ct].shape
		except:
		  continue
		
		  
	      for pd_ct in range(sim.steps['pd']-1):		#process all stuff
		
		for pd_ct_compare in range(pd_ct+1,sim.steps['pd']):
		  
		  str_print = "\t\t\t, process simulation %d/%d" % (idx_ct,np.sum(np.arange(sim.steps['pd']))*sim.steps['co']*sim.steps['in'])
		  sys.stdout.write('\r{0}'.format(str_print))
		  sys.stdout.flush()
		  
		  if self.Const['special'] == 2:
		    distanceVector = (statesTmp[pd_ct]-statesTmp[pd_ct_compare])[:,:self.Const['N']]
		  else:
		    if (self.analysis['measureStates'][co_ct][pd_ct] == 0).all() or (self.analysis['measureStates'][co_ct][pd_ct_compare] == 0).all():
		      if self.Const['special'] in [3,30]:
			self.analysis['T_conv'][co_ct][in_ct][pd_ct][pd_ct_compare] = np.nan
			self.analysis['T_conv_measure'][co_ct][in_ct][pd_ct][pd_ct_compare] = np.nan
		      else:
			print 'fix reading here!'
			#print self.analysis['measureStates'][co_ct]
			#print self.analysis['measureStates'].shape
			#print link_data
		      idx_ct += 1
		      continue
		      
		    distanceVector = (self.analysis['measureStates'][co_ct][pd_ct]-self.analysis['measureStates'][co_ct][pd_ct_compare])[:,:self.Const['N']]
		  distance = np.sqrt(np.sum(distanceVector**2,axis=1))
		  #for t in range(time_len):
		    #distanceOrthVector[t] = OrthoVector(distanceVector[t])
		
		  #distanceOrth = np.sqrt(np.sum(distanceOrthVector[:time_len]**2,axis=1))
		  
		  if self.Const['special'] == 2:
		    
		    for t in range(time_len):
		      distanceOrthVector[t] = OrthoVector(distanceVector[t])
		    distanceOrth = np.sqrt(np.sum(distanceOrthVector[:time_len]**2,axis=1))
		    
		    distance_med = sp.signal.medfilt(distanceOrth,5)
		    try:
		      idx_end = np.where(distance_med<0.1)[0][0]
		    except:
		      idx_end = -1
		    
		    
		    # run first fit to obtain guess
		    try:
		      popt,pcov = curve_fit(func_const,self.analysis['measureTimes'][:idx_end],distance_med[:idx_end])
		    
		      # cut initial and final behaviour (due to transient behaviour or convergence)
		      idx_start = np.where(distance_med<popt[0])[0][0]
		      idx_end = np.where(distance_med>popt[0])[0][-1]
		      if self.Const['N'] <= 1000:
			thresh_idx = 100
		      else:
			thresh_idx = 10
		      if idx_end - idx_start < thresh_idx:
			#print "no decorrelation could be found, end-start=%d" % (idx_end-idx_start)
			self.analysis['D_decorr'][:,idx_ct] = np.nan
		      else:		      
			# run second fit to obtain real values
			popt2,pcov2 = curve_fit(func_const,self.analysis['measureTimes'][idx_start:idx_end],distance_med[idx_start:idx_end])
			perr2 = np.sqrt(np.diag(pcov))
			#print " D = %g pm %g" % (popt2[0],perr2[0])
			self.analysis['D_decorr'][0,idx_ct] = popt2[0]
			self.analysis['D_decorr'][1,idx_ct] = perr2[0]
			
			#fig = plt.figure(figsize=(4,3))
			#ax0 = plt.axes([0.2,0.6,0.75,0.35])
			#ax1 = plt.axes([0.2,0.2,0.75,0.35])
			#ax0.plot(self.analysis['measureTimes'],distance_med,color='k')
			#ax0.plot(self.analysis['measureTimes'],distance,color='orangered')
			
			#ax0.plot(self.analysis['measureTimes'],popt[0]*np.ones(len(distance_med)),'--k')
			#ax0.plot(self.analysis['measureTimes'],popt2[0]*np.ones(len(distance_med)),'--r')
			##plt.plot(self.analysis['measureTimes'][idx_start],distance_med[idx_start],'or')
			##plt.plot(self.analysis['measureTimes'][idx_end],distance_med[idx_end],'or')
			#ax0.set_yscale('log')
			#ax1.hist(distanceVector[0],bins=101,histtype='step')
			#ax1.set_xlim([-1,1])
			#ax1.set_ylim([0,50])
			#ax1.set_ylabel([])
			#plt.ion()
			#plt.show(block=False)
			#for t in range(0,100,2):
			  #ax1.cla()
			  #ax1.hist(distanceVector[t],bins=101,histtype='step',color='k')
			  #ax1.set_xlim([-1,1])
			  #ax1.set_ylim([0,50])
			  #ax1.set_ylabel([])
			  #fig.canvas.draw()
			  
			#ax1.hist(distanceVector[0],bins=101,histtype='step',color='lightgrey')
			#ax1.hist(distanceVector[10],bins=101,histtype='step',color='grey')
			#ax1.hist(distanceVector[200],bins=101,histtype='step',color='k')
			#ax1.set_xlim([-1,1])
			#plt.show()
		    except:
		      self.analysis['D_decorr'][:,idx_ct] = np.nan
		    
		    
		    
		    #if self.analysis['D_decorr'][0,idx_ct] > 5:
		    #print np.max(distanceOrthVector,axis=0)
		    
		    
		      
		      
		  if self.Const['special'] in [3,4,30]:
		    
		    if self.Const['special'] == 4:
		      
		      for t in range(time_len):
			distanceOrthVector[t] = OrthoVector(distanceVector[t])
		      
		      self.analysis['distanceOrth'][idx_ct][:time_len] = np.sqrt(np.sum(distanceOrthVector[:time_len]**2,axis=1))
		      
		      print 'he'
		      #print assign['train'][pd_ct]
		      if not (np.isnan(assign['train'][pd_ct])).all() and not (np.isnan(assign['train'][pd_ct_compare])).all():
			print 'tada'
			deltaT = assignSpikes(assign,pd_ct,pd_ct_compare,idx_rec_cut_tmp,Topo=Topo,find_decorr=False)
			self.analysis['trainDeltaT'][idx_ct,:idx_rec_cut_tmp,0] = assign['train'][pd_ct,:idx_rec_cut_tmp,0]
			self.analysis['trainDeltaT'][idx_ct,:idx_rec_cut_tmp,1] = assign['train'][pd_ct,:idx_rec_cut_tmp,1]
			self.analysis['trainDeltaT'][idx_ct,:idx_rec_cut_tmp,2] = deltaT[:idx_rec_cut_tmp]
		      
		      for n in range(self.Const['N']):
			output = np.corrcoef(self.analysis['measureStates'][co_ct,pd_ct][:,n],self.analysis['measureStates'][co_ct,pd_ct_compare][:,n])
			self.analysis['phi_corr'][idx_ct,n] = output[0,1]
		      
		    #print 'find how trajectories collapse onto each other -> all at once or first onto some batches, later together on other ones?'
		    if distance[-1] < thresh_drop:
		      
		      idx_drop = np.where(distance < thresh_drop)[0][0]
		      #this is, where trajectories are close in distance
		      t = Data['measureTimes'][idx_drop]

		      
		      if self.Const['special'] in [3,30]:
			
			dt_tmp = np.log(distance[idx_drop]/self.analysis['D_decorr'][0])*(-0.01)
			self.analysis['T_conv'][co_ct][in_ct][pd_ct][pd_ct_compare] = max(0,t - dt_tmp)
			self.analysis['T_conv_measure'][co_ct][in_ct][pd_ct][pd_ct_compare] = t
			
			#print '\tdt: ', dt_tmp

			#idx_RC = np.argmin(abs(Data['measureTimes']-(t-dt_tmp)))
			#plt.figure(figsize=(2,2))
			#ax0 = plt.axes([0.34,0.22,0.65,0.7])
			#ax0.plot(self.analysis['measureTimes'],distance)
			##plt.plot([t-dt_tmp,t-dt_tmp],[0.0001,1],'r')
			#ax0.plot(Data['measureTimes'],self.analysis['D_decorr'][0]*np.ones(len(Data['measureTimes'])),':r')
			
			#b = distance[idx_drop] * np.exp(Data['measureTimes'][idx_drop]/0.01)
			##print 'b: ', b
			#ax0.plot(Data['measureTimes'],b*np.exp(-Data['measureTimes']/0.01),'--r')
			#ax0.plot(Data['measureTimes'][idx_RC],distance[idx_RC],'or')
			#ax0.plot(Data['measureTimes'][idx_drop],distance[idx_drop],'or')
			#ax0.set_xticks(np.linspace(0,0.5,6))
			#ax0.set_xlim([0,0.25])
			#ax0.set_yscale('log')
			#ax0.set_ylim([10**(-6),10**(1)])
			#ax0.set_yticks(np.logspace(-5,1,4))
			#ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
	
			#plt.show()
			
		    elif self.Const['special'] in [3,30]:
		      self.analysis['T_conv'][co_ct][in_ct][pd_ct][pd_ct_compare] = np.inf
		      self.analysis['T_conv_measure'][co_ct][in_ct][pd_ct][pd_ct_compare] = np.inf
		    
		  
		  elif self.Const['special'] == 5:
		    ## can changes in decorrelation distance be explained by altered temporal fluctuations?
		    
		    print var
		    print np.mean(var)
		    for n in range(self.Const['N']):
		      print "mean distance: %g" % np.mean(abs(distanceVector[:,n]))
		      print "variance: %g" % np.mean(var,axis=0)[n]
		  
		 

		    #if not co_ct:
		      #idx_tmp = idx_ct
		    #else:
		      #idx_tmp = idx_ct%co_ct
		    #self.analysis['distanceVector'][co_ct][idx_tmp][:] = distanceVector
		    
		  idx_ct += 1		
		
		
		
		
		#for t in range(time_len):
		  #print np.histogram(Data_pert['measureStates'][t][:self.Const['N']],bins=hist_bin,range=[-1.1,1])
		  #self.analysis['phase_hist'][idx_ct][t] = np.histogram(Data_pert['measureStates'][t][:self.Const['N']],bins=hist_bin,range=[-0.5,1])[0]
		  
	      #if self.Const['special'] == 3 and self.topo == 'p':		# find totally stable nets - only makes sense in poisson networks
		#idx_time = len(Data_pert['measureTimes'])-1
		#dist_vec = np.linalg.norm(OrthoVector((Data_unpert['measureStates'][idx_time] - Data_pert['finalStates'])[:self.Const['N']]))
		
		#if dist_vec < 10**(-10):
		  #div_track[co_ct][pd_ct] = 1
		  #self.analysis['T_conv'][co_ct][pd_ct] = Data_pert['TC']
		#else:
		  #div_track[co_ct][pd_ct] = 0
		  #self.analysis['T_conv'][co_ct][pd_ct] = -1
	      
	      
	  else:
	    link_unperturbed = ''.join(ncid.variables['links'][co_ct,in_ct,0,-1])
	    try:
	      Data_unpert = readDataOut(link_unperturbed,1)
	    except:
	      print 'Data is not here anymore! Reading stopped!'
	      ncid.close()
	      self.break_it = 1
	      return
	    if sim.mode == 'mosaic':
	      ct_fluxtube = 0
	      FluxtubeState = np.zeros((sim.steps['pd'],self.Const['N']))	# initialize array of maximum possible number of fluxtubes
	      FluxtubeState[ct_fluxtube] = Data_unpert['finalStates'][:self.Const['N']]
	      
	      n = int(np.sqrt(sim.steps['pd']))
	      fluxtube = np.zeros((n,n))
	      fluxtube[:] = np.nan
	    
	    for pd_ct in range(sim.steps['pd']):
	      link_perturbed = ''.join(ncid.variables['links'][co_ct,in_ct,pd_ct,-2])

	      try:
		Data_pert = readDataOut(link_perturbed,1)
		
		if sim.mode == 'mosaic':
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
		    
		#elif sim.mode == 'samples':
		  
		  #if not (self.Const['special'] == 3):
		    
		    
		    #distanceVector_min = (Data_unpert['measureStates']-Data_pert['measureStates'])[:time_len].transpose()[inDeg_min_idx].transpose()
		    #distanceVector_max = (Data_unpert['measureStates']-Data_pert['measureStates'])[:time_len].transpose()[inDeg_max_idx].transpose()
		    
		    #for t in range(time_len):
		      #distanceOrthVector[t] = OrthoVector(distanceVector[t])
		      #if self.topoConst['alpha_conv'] > 0:
		      #distanceOrthVector_min[t] = OrthoVector(distanceVector_min[t])
		      #distanceOrthVector_max[t] = OrthoVector(distanceVector_max[t])
		    #distanceOrth = np.sqrt(np.sum(distanceOrthVector**2,axis=1))
		    #self.analysis['distanceOrth'][idx_ct] = np.sqrt(np.sum(distanceVector**2,axis=1))
		    
		    #if self.topoConst['alpha_conv'] > 0:
		    #distanceOrthVector_min[t] = OrthoVector(distanceVector_min[t])
		    #distanceOrthVector_max[t] = OrthoVector(distanceVector_max[t])
		  
		  
		  #if self.topoConst['alpha_conv'] > 0:
		  #self.analysis['distanceOrth_min'][idx_ct] = np.sqrt(np.sum(distanceOrthVector_min**2,axis=1)/len(Topo['inDegree'][inDeg_min_idx])*self.Const['N'])
		  #self.analysis['distanceOrth_max'][idx_ct] = np.sqrt(np.sum(distanceOrthVector_max**2,axis=1)/len(Topo['inDegree'][inDeg_max_idx])*self.Const['N'])
		  
		#if 'special' in self.Const.keys():		  
		  
		  #if plot:
		    #levs = range(160)
		    #rb = mcolors.LinearSegmentedColormap.from_list(name='red_black',colors=[(0,(0,0,0)),(1,(1,0,0))],N=len(levs)-1,)
		  
		    #fig,axes = plt.subplots(nrows=2,ncols=1)
		    
		    #axes[0].scatter(self.analysis['trainTime'][0][idx_ct],self.analysis['trainNeuron'][0][idx_ct],c=np.log(1+self.analysis['trainDeltaTime'][idx_ct]),cmap=rb,vmin=0,vmax=max(np.log(1+self.analysis['trainDeltaTime'][idx_ct]))/100,marker='D')
		    #axes[0].set_xlim([0,1])#self.analysis['measureTimes'][2*fit_idx[1]-1]])
		    #axes[0].set_ylim([0,np.max(self.analysis['trainNeuron'])])
		    
		    ##print max(self.analysis['distanceOrth'][idx_ct] - self.analysis['distance'][idx_ct])
		    #axes[1].plot(self.analysis['measureTimes'],self.analysis['distanceOrth'][idx_ct],'-k')
		    #axes[1].plot(self.analysis['measureTimes'],self.analysis['distance'][idx_ct],'-r')
		    #if self.analysis['distanceOrth'][idx_ct][-1] < 0.01:
		      #axes[1].plot(self.analysis['measureTimes'][fit_idx],self.analysis['distanceOrth'][idx_ct][fit_idx],'or')
		      #axes[1].plot(self.analysis['measureTimes'][fit_idx],self.analysis['distance'][idx_ct][fit_idx],'or')

		    #axes[1].set_xlim([0,1])#self.analysis['measureTimes'][2*fit_idx[1]-1]])
		    ##axes[1].set_xlim([0,self.analysis['measureTimes'][-1]])
		    #axes[1].set_ylim([10**(-16),10**2])
		    #axes[1].set_yscale('log')
		    
		    #plt.show(block=False)
		    

	      except:
		if errorMsg:
		  print linkName + ' does not exist. Something seems to have gone wrong'

		
	      idx_ct += 1
      ncid.close()
      
      if sim.mode == 'eft':
	if correct_ct:
	  self.DataSv['p'][ps_ct] = 1 - np.nansum(self.DataSv['div_track'][ps_ct])/correct_ct
	else:
	  self.DataSv['p'][ps_ct] = np.nan

	self.DataSv['correct_ct'][ps_ct] = correct_ct
	
	if sim.steps['ps'] > 1:
	  print 'Fraction of restored perturbations at perturbation size ps = %5.3g is %5.3g (%d measurements)' % (self.DataSv['pert_size'][ps_ct],self.DataSv['p'][ps_ct],correct_ct)
    
      if sim.mode == 'eft_tmp':
	for idx_time_pert in range(21):
	  mask = np.invert(np.isnan(self.DataSv['div_track'][idx_time_pert,idx_ps]))
	  #print "correct: ", correct_ct
	  #correct_ct = np.sum(mask)
	  #print "correct: ", correct_ct
	  
	  
	  if correct_ct:
	    p = np.sum(self.DataSv['div_track'][idx_time_pert,idx_ps][mask])/correct_ct
	  else:
	    p = np.nan
	  
	  #if sim.steps['ps'] > 1:
	  print 'Time %d: Fraction of restored perturbations at perturbation size ps = %5.3g is %5.3g (%d measurements)' % (idx_time_pert,self.DataSv['pert_size'][ps_ct],p,correct_ct)
	  
	  #self.DataSv['pert_size'][ps_ct] = self.DataSv['pert_size'][ps_ct]
	  #self.DataSv['div_track'][idx_time_pert][ps_ct] = div_track[idx_time_pert][idx_ps]
	  #print self.DataSv['div_track'][idx_time_pert,ps_ct]
	  self.DataSv['p'][idx_time_pert,ps_ct] = p
	  #self.DataSv['CDBTimes'][idx_time_pert][ps_ct] = times[idx_time_pert][idx_ps]
	  #self.DataSv['distances'][idx_time_pert,ps_ct] = dist[idx_time_pert][idx_ps]
	  self.DataSv['correct_ct'][ps_ct] = correct_ct
      print '\n'
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
    
    #elif sim.mode == 'samples' and ('convergeVelocity' in self.analysis.keys()):
      # remove these - too costly!!
      
      
      #for dat in os.listdir(self.Path['results']):
	#if not ('results' in dat) and not ('filenames' in dat):
	  #os.remove(self.Path['results'] + dat)
      #self.analysis.pop('distance2')
      #self.analysis.pop('distanceVector')
      #self.analysis.pop('distanceOrthVector')

    if sim.mode == 'samples':
      
      #print self.analysis
      if self.Const['special'] in [2,3,30]:
	for key in ['distanceOrth','measureStates']:
	  try:
	    trash = self.analysis.pop(key)
	  except:
	    1
	#trash = self.analysis.pop('distanceOrth')
	#trash = self.analysis.pop('measureTimes')
	
	#if self.Const['special'] == 3:
	  #trash = self.analysis.pop('measureStates')
	  #trash = self.analysis.pop('train')
	
	
	#self.analysis['p'] = np.nanmean(div_track)
      
	#self.DataSv = {'p':self.analysis['p'],'T_conv':self.analysis['T_conv']}
	#save_analysis(self.DataSv,self.Path['results']+'results_processed.nc')
    
    #if sim.mode == 'eft_tmp':
      #save_analysis(self.DataSv,self.Path['results']+'results_processed.nc')
    
    
  
  def eft_calc(self,analysis):	# add proper error estimation (from eft calcs for different steps)
    
    mask = np.where(self.DataSv['p'] > 0)[0]
    
    print self.DataSv['p'][mask]
    
    x = np.log(analysis['pert_size'])
    y = np.log(analysis['p'])
    
    err = np.log(analysis['p_err'])
    #print y
    #print err
    
    #analysisFactory = IAnalysisFactory.create()
    #treeFactory = analysisFactory.createTreeFactory()
    #tree = treeFactory.create()

    #### DataPointSet
    #dataPointSetFactory = analysisFactory.createDataPointSetFactory(tree)
    #dataPointSet = dataPointSetFactory.create('dataPointSet', 'Asymmetric', 2)
    
    #print 'fill dataset'
    #### Fill DataPointSet
    #for i in range(len(analysis['p'])):
      #### Asymmetric errors
      #dataPoint = dataPointSet.addPoint()
      #dataPoint.coordinate(0).setValue(x[i])
      #dataPoint.coordinate(1).setValue(y[i])
      
      #print err[i,0]
      #print err[i,1]
      
      #dataPoint.coordinate(1).setErrorMinus(err[i,0])
      #dataPoint.coordinate(1).setErrorPlus(err[i,1])
      
      
    #fitFactory = analysisFactory.createFitFactory()
    #fitter = fitFactory.createFitter()
    #fitResult = fitter.fit(dataPointSet, 'P1')
    #print '### Asymmetric:  p0 + p1 * x'
    #print 'Parameter:', fitResult.fittedParameters()
    #print 'Error    :', fitResult.errors()

      
      
    var_p = None
    #print 'he'
    #print self.DataSv['pert_size'][mask]
    
    if (self.DataSv['p']==0).all():
      eft, var_eft = 0, 0
      var_p = np.zeros(len(self.DataSv['p']))
    elif (self.DataSv['p'][mask]==1).all():
      eft = self.ParaSim['D_decorr']
      var_eft = np.nan
      var_p = np.zeros(len(self.DataSv['p']))
    else:
      if not type(var_p)==np.ndarray:	# rather take variance error
	var_p = np.sqrt(self.DataSv['correct_ct']*self.DataSv['p']*(1-self.DataSv['p']))/self.DataSv['correct_ct']
      #try:
	#idx_stop = np.where(self.DataSv['p'] > 0.05)[0][-1]
      #except:
      idx_stop = len(self.DataSv['p'])
      #try:
	#idx_start = np.where(self.DataSv['p'] < 0.9)[0][0]
      #except:
      idx_start = 0
      
      p_fit = self.DataSv['p'][mask][idx_start:idx_stop]
      #print p_fit
      try:
	m, b, r_value, p_value, std_err = stats.linregress(self.DataSv['pert_size'][mask][idx_start:idx_stop],np.log(p_fit))

	eft = np.abs((np.log(np.exp(-1)) - b)/m)
	
	var_eft = np.abs((np.log(np.exp(-1)) - b)/(m**2)*std_err)		# get proper error for eft! only slope-error taken into account here
	if eft > self.ParaSim['D_decorr']:
	  eft = D_decorr
	  var_eft = np.nan
	
      except:
	eft = None
	var_eft = None
	print 'perturbation sizes are not in range of simulated data'
      #print eft
    return eft, var_eft, var_p
    
  def eft_calc_time_pert(self,var_p=None):	# add proper error estimation (from eft calcs for different steps)
    
    time_pert = len(self.DataSv['p'])
    eft = np.zeros(time_pert)
    var_eft = np.zeros(time_pert)
    
    if (self.DataSv['p']==0).all():	# add here: maximum distance, calculated by c++ solver
      eft[:] = np.nan
      var_eft[:] = np.nan
      var_p = np.zeros((time_pert,self.DataSv['p'].shape[1]))
      
    else:
      for idx_time_pert in range(time_pert):
	#if not type(var_p)==np.ndarray:	# rather take variance error
	  #var_p = np.sqrt(self.DataSv['correct_ct']*self.DataSv['p']*(1-self.DataSv['p']))/self.DataSv['correct_ct']
	try:
	  idx_stop = np.where(self.DataSv['p'][idx_time_pert] < 0.1)[0][0]
	except:
	  idx_stop = len(self.DataSv['p'][idx_time_pert])
	try:
	  idx_start = np.where(self.DataSv['p'][idx_time_pert] < 0.8)[0][0]
	except:
	  idx_start = 0
	mask = np.where(self.DataSv['p'][idx_time_pert] > 0)
	p_fit = self.DataSv['p'][idx_time_pert][mask][idx_start:idx_stop]
	#print p_fit
	if all(p_fit==1):
	  eft[idx_time_pert] = 10	# change to max decorr_dist
	try:
	  m, b, r_value, p_value, std_err = stats.linregress(self.DataSv['pert_size'][mask][idx_start:idx_stop],np.log(p_fit))

	  eft[idx_time_pert] = np.abs((np.log(1-np.exp(-1)) - b)/m)
	  var_eft[idx_time_pert] = np.abs((np.log(1-np.exp(-1)) - b)/(m**2)*std_err)		# get proper error for eft! only slope-error taken into account here
	  #eft[idx_time_pert] = -(1 + b)/m
	  #var_eft[idx_time_pert] = (1 + b)/(m**2)*std_err		# get proper error for eft! only slope-error taken into account here
	  

	  if eft[idx_time_pert] > self.ParaSim['D_decorr']:
	    eft[idx_time_pert] = self.ParaSim['D_decorr']
	
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
  
  #get first estimate by slope from small intervall
  steps_est = np.where(x>(x[fix_idx]+sign*0.05))[0][0] - fix_idx
  
  popt, pcov = curve_fit(func_linear,x[fix_idx:fix_idx+steps_est:sign],y[fix_idx:fix_idx+steps_est:sign])
  perr = np.sqrt(np.diag(pcov))
  m_est = popt[0]
  #m_est, b_est, r_value, p_value, std_err = stats.linregress(x[fix_idx:fix_idx+steps_est:sign],y[fix_idx:fix_idx+steps_est:sign])
  
  fit_err = np.zeros(final_idx)
  min_err = 1
  min_ct = 0
  min_idx = 0
  j = test_idx
  
  if (final_idx - first_idx <= 10):
    m = np.nan
    b = np.nan
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
  
  return m, b, min_idx

  
  
def sv_links(sim,Path,pert,linkNr=2):

  datout_sv = netcdf.netcdf_file(Path['results_links'] + '_' + str(sim.steps['ps_ct']),'w')	# save filenames to this

  datout_sv.createDimension('one',1)
  datout_sv.createDimension('connectivities',sim.steps['co'])
  datout_sv.createDimension('initial conditions',sim.steps['in'])
  
  datout_sv.createDimension('Simulations',sim.steps['co']*sim.steps['in']*sim.steps['pd'])
  datout_sv.createDimension('links',linkNr)
  datout_sv.createDimension('linksizes',len(Path['results'])+52)	#hashnames are 51 long
  
  if linkNr==2:
    svLinks = datout_sv.createVariable('links','S1',('connectivities','initial conditions','links','linksizes'))
    svIDs = datout_sv.createVariable('IDs','i',('connectivities','initial conditions'))
  elif linkNr==3:
    datout_sv.createDimension('perturbation directions',sim.steps['pd'])
    svLinks = datout_sv.createVariable('links','S1',('connectivities','initial conditions','perturbation directions','links','linksizes'))
    svIDs = datout_sv.createVariable('IDs','i',('connectivities','initial conditions','perturbation directions'))
  
  #if sim.cluster['q']:
    #svIDs = datout_sv.createVariable('IDs','i',('connectivities','initial conditions','perturbation directions'))
  
  #else:
    #svIDs = None
  
  svPert = datout_sv.createVariable('pert_size', 'f4', ('one',))
  svPert[0] = pert
  
  return datout_sv, svLinks, svIDs
  
    
def OrthoVector(vector):
  N = len(vector)
  norm_traj = np.ones(N)/np.linalg.norm(np.ones(N))	# norm vector of trajectory
  shift = np.dot(vector,norm_traj)
  #vector -= shift*norm_traj		# ... substract projection on trajectory (changes input vector as well!)
  #print "shift: %g" % shift
  return vector - shift*norm_traj
  #return vector


  

def analyze(sim_mode,net=None,sim=None,keys=None,restrictions=None,plot=0,call=0,save=0,reread=0):

  scratch = 0
  if not net:
    scratch = int(raw_input('You want to access data from local (0) or from scratch (1)? '))
    if scratch:
      sessionPath = '/scratch%02d/aschmidt/data/' % scratch	#path for storage on scratch
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
    net.data_structure(sim,1,call)
    net.initDataSv(sim)
    
    if sim_mode in ['eft','eft_tmp']:
      out = analyze_eft(net,sim,plot,save,reread)
    elif sim_mode == 'LE':
      out = analyze_LE(net,sim,plot,save)
    elif sim_mode == 'statistics':
      out = analyze_statistics(net,sim,plot,save,reread)
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
  results_file = net.Path['results'] + 'results_processed.nc'
  results_file_test = net.Path['results'] + 'test.nc'
  
  if not os.path.exists(results_file) and os.path.exists(results_file_test):
    os.system('mv %s %s' % (results_file_test,results_file))
    print 'renaming %s to %s' % (results_file_test,results_file)
  print "Hash: ", net.Hash
  #results_file = net.Path['results'] + 'test.nc'
  # construct plot strings
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
  #print 'reading  %s' % results_file
  #return
  # read from already processed data
  if not os.path.exists(results_file) or reread == 2:
    print "not yet processed"
    net.prob_read(sim)
    if not net.break_it:
      net.fluxSimSv(sim.steps)
  else:
    #print 'reading  %s' % results_file
    ncid = netcdf.netcdf_file(results_file,'r',mmap=False)
    try:
      sim.steps['ps'] = ncid.dimensions['ps_steps']
      assert sim.steps['ps']
    except:
      sim.steps['ps'] = len(ncid.variables['pert_size'][:])
    
    #try:
      #net.DataSv['pert_size'] = ncid.variables['pert_size'][:]
    #except:
      #net.DataSv['pert_size'] = ncid.variables['perturbation_sizes'][:]
    
    for key in ncid.variables.keys():
      net.DataSv[key] = ncid.variables[key][:]

    ncid.close()
    
  #net.clean()	#clean up afterwards
  
  if sim.mode == 'eft_tmp':
    print net.DataSv['p']
    time_pert = len(net.DataSv['p'])
    print 'time: ', time_pert
  
  sim.steps['ps'] = np.sum(net.DataSv['correct_ct'] > 0)

  net.eft_exp()		# mean fluxtube radii (from analysis & from finding)
  
  #if reread and sim.steps['ps']:
    ##print 'fix and re-enable saving (and cleaning) after re-reading'
  #print "CDB times: ", net.DataSv['CDBTimes']
  
  if sim.mode == 'eft':
    #print net.DataSv
    #print net.DataSv.keys()
    #print "over 9000(!!!): ", np.sum(net.DataSv['CDBTimes'] > 0.5)
    #print "MAX POWER(!!!): ", np.nanmax(net.DataSv['CDBTimes'])
    if not ('eft_bs' in net.DataSv.keys()) or reread == 1:
      #if not 'p_conv' in net.DataSv.keys():
      analysis['pert_size'] = net.DataSv['pert_size'][:sim.steps['ps']]
      
      p_conv = np.zeros((sim.steps['ps'],sim.steps['total']))
      p_conv[:] = np.nan
      
      analysis['p'] = np.zeros((sim.steps['ps'],3))
      analysis['p'][:] = np.nan
      
      p = np.zeros((sim.steps['ps'],sim.steps['co']))
      
      for ps_ct in range(sim.steps['ps']):
	mask = np.invert(np.isnan(net.DataSv['div_track'][ps_ct]))
	p_conv[ps_ct][:np.sum(mask)] = np.cumsum(net.DataSv['div_track'][ps_ct][mask])/(np.arange(np.sum(mask))+1).astype('float')
	
	for co_ct in range(sim.steps['co']):
	  mask = np.invert(np.isnan(net.DataSv['div_track'][ps_ct][co_ct]))
	  p[ps_ct][co_ct] = 1 - np.sum(net.DataSv['div_track'][ps_ct][co_ct][mask])/float(np.sum(mask))
	
	analysis['p'][ps_ct] = bootstrap(p[ps_ct],10000)
      print "p: ", analysis['p']
      
      eft_co = np.zeros((2,sim.steps['co']))
      eft_co[:] = np.nan
      idx_start = 0
      idx_stop = sim.steps['ps']
      
      for co_ct in range(sim.steps['co']):
	try:
	  idx_stop = np.where(p[:,co_ct]==0)[0][0]
	except:
	  idx_stop = sim.steps['ps']
	
	mask = np.invert(np.isnan(p[idx_start:idx_stop,co_ct]))
	#print p[:,co_ct]
	popt, pcov = curve_fit(func_exp,analysis['pert_size'][idx_start:idx_stop][mask],p[idx_start:idx_stop,co_ct][mask])
	perr=np.sqrt(np.diag(pcov))
	eft_co[0,co_ct] = 1./popt[0] * (- np.log(popt[1]) - 1)
	#np.abs((np.log(np.exp(-1)) - popt[1])/popt[0])
	eft_co[1,co_ct] = np.sqrt((perr[0]/popt[0]**2*(-1-np.log(popt[1])))**2 + (perr[1]/(popt[0]*popt[1]))**2)
	#print 'eft: ', eft_co[:,co_ct]
      net.DataSv['eft_bs'] = bootstrap(eft_co,10000)
      
      
      for key in analysis.keys():
	net.DataSv[key] = analysis[key]
      net.DataSv['p_conv'] = p_conv.transpose()
      
      try:
	trash = net.DataSv.pop('p_err')
      except:
	1
      save_analysis(net.DataSv,results_file)
    #else:
      
    analysis['eft'] = net.DataSv['eft_bs']
    analysis['p'] = net.DataSv['p']
    analysis['pert_size'] = net.DataSv['pert_size']
    #print 'mean eft: %g' % net.DataSv['eft_bs'][0]
    #print r'95perc. confidence: ', net.DataSv['eft_bs'][1:]
    
  
  elif sim.mode == 'eft_tmp':
    try:
      trash = net.DataSv.pop('CDBTimes')
    except:
      1
    
    if not ('eft_bs' in net.DataSv.keys()) or reread == 1:
      #if not 'p_conv' in net.DataSv.keys():
      analysis['pert_size'] = net.DataSv['pert_size']
      
      analysis['p_conv'] = np.zeros((time_pert,sim.steps['total'],sim.steps['ps']))
      p_conv = np.zeros((sim.steps['ps'],sim.steps['total']))
      p_conv[:] = np.nan
      
      net.DataSv['eft_bs'] = np.zeros((time_pert,3))
      
      analysis['p'] = np.zeros((time_pert,sim.steps['ps'],3))
      analysis['p'][:] = np.nan
      
      p = np.zeros((sim.steps['ps'],sim.steps['co']))
      print 'steps: ', sim.steps['ps']
      for t in range(time_pert):
	for ps_ct in range(sim.steps['ps']):
	  mask = np.invert(np.isnan(net.DataSv['div_track'][t,ps_ct]))
	  p_conv[ps_ct][:np.sum(mask)] = np.cumsum(net.DataSv['div_track'][t,ps_ct][mask])/(np.arange(np.sum(mask))+1).astype('float')
	  
	  for co_ct in range(sim.steps['co']):
	    mask = np.invert(np.isnan(net.DataSv['div_track'][t,ps_ct,co_ct]))
	    p[ps_ct,co_ct] = 1 - np.sum(net.DataSv['div_track'][t,ps_ct,co_ct][mask])/float(np.sum(mask))
	  
	  analysis['p'][t,ps_ct] = bootstrap(p[ps_ct],10000)
	print "p: ", analysis['p'][t]
	
	eft_co = np.zeros((2,sim.steps['co']))
	eft_co[:] = np.nan
	idx_start = 0
	idx_stop = sim.steps['ps']
	
	for co_ct in range(sim.steps['co']):
	  try:
	    idx_stop = np.where(p[:,co_ct]==0)[0][0]
	  except:
	    idx_stop = sim.steps['ps']
	  
	  mask = np.invert(np.isnan(p[idx_start:idx_stop,co_ct]))
	  #print p[:,co_ct]
	  popt, pcov = curve_fit(func_exp,analysis['pert_size'][idx_start:idx_stop][mask],p[idx_start:idx_stop,co_ct][mask])
	  perr=np.sqrt(np.diag(pcov))
	  eft_co[0,co_ct] = 1./popt[0] * (- np.log(popt[1]) - 1)
	  #np.abs((np.log(np.exp(-1)) - popt[1])/popt[0])
	  eft_co[1,co_ct] = np.sqrt((perr[0]/popt[0]**2*(-1-np.log(popt[1])))**2 + (perr[1]/(popt[0]*popt[1]))**2)
	  #print 'eft: ', eft_co[:,co_ct]
	net.DataSv['eft_bs'][t] = bootstrap(eft_co,10000)
	
	
	for key in analysis.keys():
	  net.DataSv[key] = analysis[key]
	net.DataSv['p_conv'][t] = p_conv.transpose()
      
      try:
	trash = net.DataSv.pop('p_err')
      except:
	1
      save_analysis(net.DataSv,results_file)
      
    analysis['eft'] = net.DataSv['eft_bs']
    analysis['p'] = net.DataSv['p']
    analysis['pert_size'] = net.DataSv['pert_size']
    
    
  if 'special' in net.Const.keys():
    if net.Const['special'] == 1:
      prob_convergence(net.Path['results'] + 'results_processed.nc')
	  
  
  if plot:
    sim.steps['ps'] = np.sum(net.DataSv['correct_ct']>0)
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
    if np.isnan(net.ParaSim['D_decorr']):
      D_thresh = 0.1
    else:
      D_thresh = 0.7*net.ParaSim['D_decorr']
    
    #print net.DataSv['CDBTimes'][np.where(net.DataSv['distances']>D_thresh)]
    if 'CDBTimes' in net.DataSv.keys():
      
      plt.figure(figsize=(3,2))
      
      ax = plt.axes([0.22,0.22,0.74,0.74])
      
      for ps_ct in range(sim.steps['ps']):
	# get histogram of times, after which simulations diverged
	#print net.DataSv['CDBTimes']
	#print net.DataSv.keys()
	times_histo = np.histogram(net.DataSv['CDBTimes'][ps_ct][np.where(net.DataSv['distances'][ps_ct]>D_thresh)],bins=101,range=[0,np.max(net.DataSv['CDBTimes'])])
	
	col = ps_ct/float(sim.steps['ps'])
	ax.plot(times_histo[1][:-1],1-np.cumsum(times_histo[0])/float(np.sum(net.DataSv['distances'][ps_ct]>D_thresh)),color=(col,col,col))
	#ax.plot(times_histo[1][:-1],1-np.cumsum(times_histo[0])/float(sim.steps['total']),color=(col,col,col))
	ax.set_xlabel(r'$\displaystyle T_{div}$',fontsize=14)
	ax.set_ylabel('survival rate')

      ax.set_xlim([0,1])
      #ax.set_yscale('log')
      #ax.set_ylim([10**(-3),1])
      #plt.legend()
      
      #plt.suptitle(title_str,fontsize=14,y=0.97)
      if save:
	save_str = './pics/SONET/pert_recall_N=%d_K=%d_conv=%3.1g.pdf' % (net.Const['N'],net.Const['K'],net.topoConst['alpha_conv'])
	plt.savefig(save_str)
	print "Figure saved to %s." % save_str
      plt.show(block=False)
      
    
    if sim.mode == 'eft':
      
      plt.figure(figsize=(3,2))
      ax1 = plt.axes([0.22,0.22,0.35,0.74])
      ax2 = plt.axes([0.6,0.22,0.35,0.74])
  
      #net.DataSv['distances'] = net.DataSv['distances'].reshape(np.prod(net.DataSv['distances'].shape))
      
      #ax0.hist(np.log10(net.DataSv['distances']),bins=100,range=[-8,1])
      
      #ax0.set_xlim([-8,1])
      
      net.DataSv['p_conv'] = 1-net.DataSv['p_conv'].transpose()
      
      co_conv = sim.steps['co']#5
      in_conv = sim.steps['in']#min(25,sim.steps['in'])
      pd_conv = sim.steps['pd']#100
      
      steps_conv = co_conv*in_conv*pd_conv
      
      for ps_ct in range(sim.steps['ps']):
	col = ps_ct/float(sim.steps['ps'])
	mask = np.invert(np.isnan(net.DataSv['p_conv'][ps_ct]))
	ax1.plot(np.arange(len(net.DataSv['p_conv'][ps_ct]))[mask]+1,net.DataSv['p_conv'][ps_ct][mask],color=(col,col,col),linewidth=0.5)
	p_conv_tmp = 1 - np.cumsum(net.DataSv['div_track'][ps_ct][:co_conv,:in_conv,:pd_conv])/(np.arange(steps_conv)+1).astype('float')
	mask = np.invert(np.isnan(p_conv_tmp))
	
	ax2.plot(analysis['pert_size'][ps_ct],analysis['p'][ps_ct,0],'.',color=(col,col,col),markersize=5)
      ax1.set_xscale('log')
      
      ax1.set_xticks(10**np.linspace(2,4,2))
      ax1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: r'$\displaystyle 10^{%d}$' % np.log10(x)))
      ax1.set_xlim([1,sim.steps['total']])
      ax1.set_xlabel(r'simulation \#',fontsize=12)
      #ax1.set_ylim([0.01,1])
      ax1.set_ylabel(r'$\displaystyle f_{RP}$',fontsize=14)
      
      
      # curves
      appr_plot_range = 10**np.linspace(-4,1,1000)
      appr_exp_curve = np.exp(-appr_plot_range/net.eft_exp)
      
      appr_curve = np.exp(-appr_plot_range/analysis['eft'][0])
      ax2.plot(appr_plot_range,appr_curve,'--',color='k',linewidth=0.5)
      
      eft_error = analysis['eft'][1:] - analysis['eft'][0]
      
      ax2.errorbar(analysis['eft'][0],1./np.exp(1),xerr=[[eft_error[0]],[eft_error[1]]],ecolor='r',elinewidth=0.5)
      
      ax2.plot([analysis['eft'][0],analysis['eft'][0]],[0,1],'--',color='lightgrey')
      
      mask = np.invert(np.isnan(analysis['p'][:,0]))
      plt_bootstrap(analysis['pert_size'],analysis['p'],ax2,'k','-',mask=mask)
      
      ax2.plot(appr_plot_range,appr_curve,'--',color='r',linewidth=0.5)
      ax2.tick_params(axis='y',which='both',left='off',labelleft='off')
      
      plt.annotate(r'$\displaystyle \varepsilon_{FT}$',xy=[analysis['eft'][0],np.exp(-1)],xytext=[0.1,0.5],arrowprops=dict(arrowstyle="->"),fontsize=12)
      
      #ax2.set_yscale('log')
      #ax2.plot([D_decorr,D_decorr],[0,1],'k',linewidth=2)
      #ax2.plot([D_spike_fail,D_spike_fail],[0,1],'r',linewidth=2)
      
      plt.setp(ax2,xscale='log',xticks=10**np.linspace(-3,-1,2),xlim=[10**(-4),1])
      
      ax2.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(x))))
      ax2.set_xlabel(r'$\displaystyle \varepsilon$',fontsize=14)
      ax2.set_ylim([0,1])
      
      
      #plt.suptitle(title_str,fontsize=14,y=0.97)
      if save:
	save_str = './pics/SONET/fluxtubes/convergence_N=%d_K=%d_conv=%3.1g.pdf' % (net.Const['N'],net.Const['K'],net.topoConst['alpha_conv'])
	plt.savefig(save_str)
	print "Figure saved to %s." % save_str
  
      plt.show(block=False)
      
    elif sim.mode == 'eft_tmp':
      
      plt.figure(figsize=(3,2))
      ax1 = plt.axes([0.2,0.2,0.75,0.75])
      
      for t in range(time_pert):
	col = t/float(time_pert)
	print analysis['pert_size']
	print analysis['p'][t]
	ax1.plot(analysis['pert_size'][:sim.steps['ps']],analysis['p'][t],'.',color=(col,col,col),markersize=5)
	
	mask = np.invert(np.isnan(analysis['p'][t,:,0]))
	plt_bootstrap(analysis['pert_size'][:sim.steps['ps']],analysis['p'][t],ax1,'k','-',mask=mask)
      
	#for ps_ct in range(sim.steps['ps']):
	
	
	
      ax1.set_xscale('log')
      
      ax1.set_xticks(10**np.linspace(2,4,2))
      ax1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: r'$\displaystyle 10^{%d}$' % np.log10(x)))
      ax1.set_xlim([1,sim.steps['total']])
      ax1.set_xlabel(r'simulation \#',fontsize=12)
      #ax1.set_ylim([0.01,1])
      ax1.set_ylabel(r'$\displaystyle f_{RP}$',fontsize=14)
      
      
      ##curves
      #appr_plot_range = 10**np.linspace(-4,1,1000)
      #appr_exp_curve = np.exp(-appr_plot_range/net.eft_exp)
      
      #appr_curve = np.exp(-appr_plot_range/analysis['eft'][0])
      #ax1.plot(appr_plot_range,appr_curve,'--',color='k',linewidth=0.5)
      
      eft_error = analysis['eft'][1:] - analysis['eft'][0]
      
      #ax1.errorbar(analysis['eft'][0],1./np.exp(1),xerr=[[eft_error[0]],[eft_error[1]]],ecolor='r',elinewidth=0.5)
      
      #ax1.plot([analysis['eft'][0],analysis['eft'][0]],[0,1],'--',color='lightgrey')
      
      
      
      #ax1.plot(appr_plot_range,appr_curve,'--',color='r',linewidth=0.5)
      ax1.tick_params(axis='y',which='both',left='off',labelleft='off')
      
      #plt.annotate(r'$\displaystyle \varepsilon_{FT}$',xy=[analysis['eft'][0],np.exp(-1)],xytext=[0.1,0.5],arrowprops=dict(arrowstyle="->"),fontsize=12)
      
      #ax1.set_yscale('log')
      #ax1.plot([D_decorr,D_decorr],[0,1],'k',linewidth=2)
      #ax1.plot([D_spike_fail,D_spike_fail],[0,1],'r',linewidth=2)
      
      mask = np.invert(np.isnan(analysis['eft'][1:,0]))
      print analysis['eft'][1:]
      print np.linspace(0.1,2,time_pert-1)
      sp = plt.axes([.25, .25, .3, .3])
      plt_bootstrap(np.linspace(0.1,2,time_pert-1),analysis['eft'][1:],sp,'k','-',mask=mask)
      #sp.plot(range(1,time_pert),analysis['eft'][1:],color='k')
      sp.set_ylim([10**(-2),10**(-1)])
      sp.set_yscale('log')
      
      sp.set_xticks(np.linspace(0,2,5))
      #sp.set_xticklabels(np.linspace(0,2,time_pert),fontsize=10)
      sp.set_xlabel('time in s',fontsize=14)
      sp.set_ylabel(r'$\displaystyle \varepsilon_{FT}$',fontsize=14)
      
      plt.setp(ax1,xscale='log',xticks=10**np.linspace(-3,-1,2),xlim=[10**(-4),1])
      
      ax1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(x))))
      ax1.set_xlabel(r'$\displaystyle \varepsilon$',fontsize=14)
      ax1.set_ylim([0,1])
      
      
      #plt.suptitle(title_str,fontsize=14,y=0.97)
      if save:
	save_str = './pics/SONET/fluxtubes/convergence_N=%d_K=%d_conv=%3.1g.pdf' % (net.Const['N'],net.Const['K'],net.topoConst['alpha_conv'])
	plt.savefig(save_str)
	print "Figure saved to %s." % save_str
  
      plt.show(block=False)
      
    
  try:
    print 'found flux tube radius: %5.3f+/-%5.3f' %(net.DataSv['eft_bs'][0],np.max(net.DataSv['eft_bs'][1:]-net.DataSv['eft_bs'][0]))
  except:
    print "No fluxtube radius found."

  return analysis

def dummy_plot(ccmap,border,steps):
  plt.figure()
  # Using contourf to provide my colorbar info, then clearing the figure
  Z = [[0,0],[0,0]]
  levels = np.linspace(0,border,steps)
  cbar_dummy = plt.contourf(Z, levels, cmap=ccmap)
  plt.clf()
  return cbar_dummy


def analyze_LE(net,sim,plot,save):
  print net.Hash
  try:
    net.analysis
  except:
    net.analysis = {}
  
  results_file = net.Path['results'] + 'results_processed.nc'
  
  #if (not os.path.exists(results_file)) or reread==2:
    #print "read it"
    #ncid = netcdf.netcdf_file(results_file,'r',mmap=False)
    
  #print net.Path['results']
  
  #for dat in os.listdir(net.Path['results'] + 'init/'):
    #ncid = netcdf.netcdf_file(net.Path['results'] + 'init/' + dat)
    
    #print ncid.variables.keys()
    #ncid.close()
    
    #net.prob_read(sim)					# get net.analysis
    #if not net.break_it:
      #save_analysis(net.analysis,results_file)
  #if not len(net.analysis.keys()):			# if stuff should be reread from results file
    #ncid = netcdf.netcdf_file(results_file,'r',mmap=False)
    ## hand over all variables
    #for key in ncid.variables.keys():
      ##if key in ['convergeVelocity','kinkTime','measureTimes','distanceOrth']:
	##['distance','trainTimeDrive','trainNeuronDrive','trainNeuron','trainTime']):
      #net.analysis[key] = ncid.variables[key][:]
  # read data
  Data = {}
  
  for dat in os.listdir(net.Path['inpara']):
    if 'ParaTop' in dat:
      ncid = netcdf.netcdf_file(net.Path['inpara']+dat,'r',mmap=None)
      inDeg = ncid.variables['inDegree'][:]
      ncid.close()
  
  fileDir = net.Path['results'] + 'init/'
  #print fileDir
  for dat in os.listdir(fileDir):
    if 'ini' in dat:
      #print dat
      Data = readDataOut(fileDir+dat)
      #print Data.keys()
      if net.topo == 'p':
	net.analysis['LyapunovExponents'] = Data['subLyapunovExponents']
	#print Data['LyapunovExponents']
	#print Data['subLyapunovExponents']
	
	#plt.figure()
	#plt.plot(np.arange(len(Data['subLyapunovExponents'])),Data['subLyapunovExponents'])
	#plt.show(block=False)
	#ncid = netcdf.netcdf_file(net.Path['results']+'init/'+dat,'r',mmap=None)
	
	#dimLE = ncid.dimensions['LEonsSz']
      dimConv = 0
      if 'LEconvergenceSz' in ncid.dimensions:
	print 'conv'
	dimConv = len(Data['LEconvergence'])
	#Data['LEtimes'] = ncid.variables['LEtimes'][:]
	
      #if net.topo == 'p':
	#Data['LyapunovExponents'] = ncid.variables['LEons'][:][net.Const['N']:]
	#if dimConv:
	  #Data['LEconvergence'] = ncid.variables['LEconvergence'][:].reshape((dimConv/dimLE, dimLE))[:,net.Const['N']:]
      #else:
	#Data['LyapunovExponents'] = ncid.variables['LEons'][:][net.Const['N']:]
	#if dimConv:
	  #Data['LEconvergence'] = ncid.variables['LEconvergence'][:].reshape((dimConv/dimLE, dimLE))
      
      #net.analysis['LyapunovExponents'] = Data['LyapunovExponents']
      #Data['LEconvergence'] = ncid.variables['LEconvergence'][:].reshape((dimConv/dimLE, dimLE))[:,self.Const['N']]
      
      #if 'subLyapunovExponents' in ncid.variables.keys():
	#Data['LyapunovExponents'] = ncid.variables['subLyapunovExponents'][:]
	
	#Data['LEconvergence'] = ncid.variables['subLEconvergence'][:].reshape((dimConv/dimLE, dimLE))
      #else:
	#Data['LyapunovExponents'] = ncid.variables['LEons'][:]
	#
	
      #print ncid.variables['LEons'][:]
      #print Data['LyapunovExponents'].shape
  # now get convergence of mean contribution to LE with in-degree
  #print np.mean(Data['localLE'],axis=0)
  #print 'max', np.sum(Data['LyapunovExponents'] > -1), np.sum(inDeg < 3)
  
  #print 'min', np.sum(Data['LyapunovExponents'] < -99), np.sum(inDeg < 3)
  #print net.analysis
  #print Data
  #print 'max: %g' % Data['LyapunovExponents'][1]
  #if net.Const['special'] == 0:
    #print 'mean: %g' % np.mean(Data['LyapunovExponents'])
    #print 'second: %g' % Data['LyapunovExponents'][2]
  
  

  #fig,axes = plt.subplots(nrows=2,ncols=1)
  #axes[0].plot(np.linspace(1/float(N),1,N), analysis['LyapunovExponents'])
  #axes[1].plot(np.linspace(1/float(N),1,N),np.mean(Data['localLE'],axis=0))
  #plt.show()
  
  if plot:
    N = len(Data['LyapunovExponents'])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
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
    
    #plot CLV, calculate "Chaos index"
    

    if 'localLE' in Data.keys():
      
      plt.figure(figsize=(4,3))
      
      ax0 = plt.axes([0.15,0.15,0.7,0.8])
      axcb = plt.axes([0.9,0.15,0.02,0.8])
      #gs = gridspec.GridSpec(2, 30)
    
      #ax0 = plt.subplot(gs[0:2,0:16])
      #ax1 = plt.subplot(gs[0,20:30])
      #ax2 = plt.subplot(gs[1,20:30])
      #axcb = plt.(gs[0:2,17])
      
      levs = range(160)
      assert len(levs) % 2 == 0, 'N levels must be even.'
      min_val = -1.5#Data['localLE'].min()
      max_val = 0.1#Data['localLE'].max()
      
      zero_val = abs(min_val)/(max_val-min_val)
      rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors=[(0,(0,0,1)),(zero_val,(1,1.,1)),(1,(1,0,0))],N=len(levs)-1,)
      dim0 = Data['localLE'].shape[0]
      dim1 = N#Data['localLE'].shape[1]
      times = np.tile(Data['LEtimes'][0:dim0,np.newaxis], (1,dim1))	#times = np.tile(Data['LEtimes'][0:dim0], (dim1,1))
      index = np.tile(np.linspace(1/float(dim1),1,dim1), (dim0,1))	#index = np.tile(np.linspace(1/float(dim1),1,dim1)[:,np.newaxis], (1,dim0))
      
      if Data['localLE'].shape[1] > N:
	N_idx = N
      else:
	N_idx = 0
      Data['localLE'] = Data['localLE'].transpose()[N_idx:].transpose()
	  
      pc = ax0.pcolor(times,index,Data['localLE'],cmap=rwb,vmin=min_val,vmax=max_val)
      ax0.set_xlabel('time (ms)',fontsize=14)
      ax0.set_ylabel('i / N',fontsize=14)
      ax0.set_xlim([0,np.max(Data['LEtimes'][0:dim0])])
      ax0.set_title('local LE per spike',fontsize=14)
      plt.colorbar(pc, cax=axcb)
    #else:
    #gs = gridspec.GridSpec(2, 1)
    
    if net.Const['special'] == 0:
      fig = plt.figure(figsize=(2.5,2))
      
      ax1 = plt.axes([0.3,0.2,0.65,0.75])
      #ax2 = plt.subplot(gs[1,0])
      ## get covariance of localLE with indegree/outdegre and/or firing rate
      
      if net.topo == 'p':
	ax1.plot(np.linspace(0,1,N), Data['subLyapunovExponents'],'k')
      else:
	ax1.plot(np.linspace(1/float(N),1,N-1), Data['LyapunovExponents'][1:],'k')
	
      ax1.plot([0,1], [-100,-100],'--',color='grey')
      
      #ax1.plot(np.linspace(1/float(N),1,N), Data['LEclv'][N_idx:])
      #ax1.yaxis.set_label_position('right')
      #ax1.yaxis.set_ticks_position('right')
      ax1.set_xticks(np.linspace(0,1,3))
      ax1.set_yticks(np.linspace(-100,0,6))
      
      ax1.set_ylabel(r'LE $\displaystyle \lambda_i / s^{-1}$',fontsize=14)
      ax1.set_xlabel('i/N',fontsize=12)
      ax1.set_ylim([-110,10])
    
    
    #pc = ax1.pcolor(times, index, Data['localLE'],cmap=rwb)
    #ax1.set_xlabel('time (ms)')
    #ax1.set_ylabel('i / N')
    #ax1.set_xlim([0,np.max(Data['LEtimes'][0:dim0])])
    #ax1.set_title('local Lyapunov exponents per spike')
    
    #if dimConv:
      #ax2.plot(Data['LEtimes'][:len(Data['LEconvergence'])],Data['LEconvergence'].transpose()[1])
      #ax2.plot(Data['LEtimes'][:len(Data['LEconvergence'])],Data['LEconvergence'].transpose()[N/2])
      #ax2.plot(Data['LEtimes'][:len(Data['LEconvergence'])],Data['LEconvergence'].transpose()[-1])
      #ax2.yaxis.set_label_position('right')
      #ax2.yaxis.set_ticks_position('right')
      #ax2.set_xlabel('time',fontsize=14)
      #ax2.set_ylabel('LE convergence',fontsize=14)
      #ax2.set_xlim([0, max(Data['LEtimes'])])
      ##ax2.set_title('LE convergence',fontsize=16)
      #ax2.set_ylim([-110,10])
    
    
    #ax3.set_ylim([np.min(analysis['LyapunovExponents']),np.max(analysis['LyapunovExponents'])])
    #plt.suptitle(title_str,fontsize=18)
    #r'Lyapunov spectra for N=%d, K=%d, $\displaystyle\nu=%g %s$'%(net.Const['N'],net.Const['K'],net.Const['rateWnt'],',\; '.join(str_add.split('_'))),fontsize=16)
    
    if save:
      save_str = './pics/%s/LE/single_N=%d_K=%d_f=%3.1f%s.pdf' % (net.scriptName.split('_')[0],net.Const['N'],net.Const['K'],net.Const['rateWnt'],para_str)
      plt.savefig(save_str)
      print "Figure saved to %s." % save_str
    #else:
    plt.show(block=False)
  
  return Data


def analyze_statistics(net,sim,plot,save,reread):
  
  if net.Const['special'] == 1:
    for dat in os.listdir(net.Path['results'] + 'init/'):
      if 'ini' in dat:
	#print dat
	Data = readDataOut(net.Path['results']+'init/'+dat)
	print Data.keys()
	break
    
    plot_neurons = np.where(Data['train'][:,0]<1000)[0]
    
    time_cut = 1#1./(2*net.Const['rateWnt'])
    
    idx_1sec = np.argwhere(Data['train'][:,1]>time_cut)[0][0]
    spikes = 0
    while not (spikes > 3):
      idx_neuron = np.random.randint(1000)
      spikes = np.sum(Data['train'][:,0][:idx_1sec] == idx_neuron)
    idx_spikes = np.where(Data['train'][:,0][:idx_1sec]==idx_neuron)[0]
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
  
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
    plt.figure(figsize=(6,3))
    ax1 = plt.axes([0.1,0.575,0.5,0.375])
    ax2 = plt.axes([0.1,0.125,0.5,0.375],sharex=ax1)
    ax3 = plt.axes([0.675,0.65,0.3,0.3])
    ax4 = plt.axes([0.675,0.175,0.3,0.3])
    
    spike_time = Data['train'][:,1][idx_spikes[0]]
    ax1.plot([0,time_cut],[idx_neuron,idx_neuron],'b',linewidth=6,alpha=0.4)
    ax1.plot([spike_time,spike_time],[0,1000],'y',linewidth=6,alpha=0.6)
    ax1.plot(Data['train'][:,1][plot_neurons],Data['train'][:,0][plot_neurons],'sk',markersize=2)
    ax1.plot(Data['train'][:,1][idx_spikes],Data['train'][:,0][idx_spikes],'sr',markersize=5)
    ax1.set_ylabel(r'Neuron \#',fontsize=12)
    #ax1.set_xticks([])
    ax1.yaxis.set_ticks_position('none')
    ax1.xaxis.set_ticks_position('bottom')
    #ax1.set_xlim([0,time_cut])
    
    #print Data['measureStates'].shape
    
    ax2.plot(Data['measureTimes'],Data['measureStates'][:,idx_neuron],'-k')
    
    for i in idx_spikes:
      ax2.plot([Data['train'][i,1],Data['train'][i,1]],[-1,2],'r',linewidth=1)
    ax2.plot([0,time_cut],[1,1],'--k',linewidth=0.5)
    ax2.set_xlabel('t in s',fontsize=12)
    ax2.set_ylabel(r'$\displaystyle \phi$',fontsize=14)
    ax2.set_xlim([0,time_cut])
    ax2.set_ylim([-0.2,1.2])
    ax2.set_yticks(np.linspace(0,1,3))
    ax2.spines['right'].set_color('none')
    ax2.yaxis.set_ticks_position('left')
    ax2.spines['top'].set_color('none')
    ax2.xaxis.set_ticks_position('bottom')
    
    rate_var = np.var(Data['rateNeurons'])
    rate_array = np.linspace(0,4*net.Const['rateWnt'],10000)
    hist_rate = np.histogram(Data['rateNeurons'],bins=50,range=[0,4*net.Const['rateWnt']])
    hist_vals = hist_rate[0]/float(net.Const['N'])
    mu = hist_rate[1][np.argmax(hist_vals)]
    gauss = np.max(hist_vals)/np.sqrt(2*math.pi*rate_var)*np.exp(-(rate_array-mu)**2/(2*rate_var))
    ax3.plot(hist_rate[1][:-1],hist_vals,'k')
    ax3.plot(rate_array,gauss,'--',color='lightgrey')
    ax3.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '%d'%y))
    ax3.set_xticks(np.linspace(0,4*net.Const['rateWnt'],5))
    ax3.set_xlim([0,4*net.Const['rateWnt']])
    ax3.set_yticks([])
    ax3.set_ylim([0,1.2*np.max(hist_vals)])
    ax3.set_xlabel(r'$\displaystyle \langle ISI \rangle = \nu$ in Hz',fontsize=14)
    #ax3.set_ylabel(r'$\displaystyle \rho(\nu)$',fontsize=14)
    ax3.yaxis.set_label_position('right')
    ax3.spines['left'].set_color('none')
    ax3.spines['right'].set_color('none')
    ax3.yaxis.set_ticks_position('right')
    ax3.spines['top'].set_color('none')
    ax3.xaxis.set_ticks_position('bottom')

    
    hist_CV = np.histogram(Data['cvNeurons'],bins=50,range=[0,1.7])
    idx_CV_max = np.argmax(hist_CV[0])
    
    #print idx_CV_max
    #print hist_CV[1][idx_CV_max]
    ax4.plot([hist_CV[1][idx_CV_max],hist_CV[1][idx_CV_max]],[0,10000],'b',linewidth=6,alpha=0.4)
    ax4.plot(hist_CV[1][:-1],hist_CV[0],'k')
    ax4.set_xticks(np.linspace(0,1.5,4))
    ax4.set_xlim([0,1.7])
    ax4.set_ylim([0,1.2*np.max(hist_CV[0])])
    ax4.set_yticks([])
    ax4.set_xlabel(r'$\displaystyle CV_{ISI}$',fontsize=14)
    #ax4.set_ylabel(r'$\displaystyle \rho(CV)$',fontsize=14)
    ax4.yaxis.set_label_position('right')
    ax4.spines['left'].set_color('none')
    ax4.spines['right'].set_color('none')
    ax4.yaxis.set_ticks_position('right')
    ax4.spines['top'].set_color('none')
    ax4.xaxis.set_ticks_position('bottom')

    
    print "now get synchrony!"
    print Data['chi']
    plt.show()
    
    return
    
  analysis = {}
  Data = {}
  stats = ['topology','moments','subthreshold']		# 
  #stats = ['topology','moments','ISIdecorr','subthreshold','trainCorr'] #, choose which ones to plot
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  numsim = sim.steps['total']
  N = net.Const['N']
  
  str_add = []
  for item in net.topoConst.items():
    if 'alpha' in item[0] and not (net.topo == 'p'):
      str_add.append('%s=%g' % (item[0],item[1]))
    elif not ('alpha' in item[0]):
      str_add.append('%s=%g' % (item[0],item[1]))
  
  if os.path.exists(net.Path['results'] + 'results_processed.nc') and not reread:
    try:
      ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r',mmap=None)
      print ncid.variables.keys()
      for key in ncid.variables.keys():
	analysis[key] = ncid.variables[key][:]
	#Data[key] = ncid.variables[key][:]
	
      ncid.close()
    except:
      print 'results file corrupted - read again'
  
  else:
  
    read_list = ['chi','rateNeurons','finalCurrents']
    read_list_topo = []
    
    # specify values that should be read out
    #read_list = ['rateNeurons', 'cvNeurons','chi','finalCurrents','mean_ISIdecorr', 'cv_ISIdecorr','trainNeuron','trainTime','measure1stStateVarNeurons','measureTimes']
    #read_list_topo = ['p_hat', 'alpha_recip_hat', 'alpha_conv_hat', 'alpha_div_hat', 'alpha_chain_hat', 'inDegree', 'outDegree']
    
    if 'topology' in stats:
      read_list_topo.extend(['p_hat', 'alpha_recip_hat', 'alpha_conv_hat', 'alpha_div_hat', 'alpha_chain_hat', 'inDegree', 'outDegree'])
    if 'moments' in stats:
      read_list.append('cvNeurons')
    if 'ISIdecorr' in stats:
      read_list.extend(['mean_ISIdecorr','cv_ISIdecorr'])
    if 'subthreshold' in stats:
      read_list.extend(['measureStates','measureTimes'])
    
    idx_file = 0
    
    
    print 'read stuff...'
    pathResults = net.Path['results'] + 'init/'
    
    for dat in os.listdir(pathResults):
      if 'ini' in dat:
	
	Data, ncidTopo = readDataStats(pathResults + dat,Data,read_list,numsim,idx_file)	# read results file
	Data = readDataStats(ncidTopo,Data,read_list_topo,numsim,idx_file)				# read topology file
	
	if not idx_file and ('ISI_decorr' in stats):
	  Data['mean_ISIdecorr'] = np.zeros((numsim,N))
	  
	if ('trainCorr' in stats) and (not ('Xcor' in analysis.keys()) or reread):
	  
	  if not idx_file:
	    bins = 1000
	    bin_hist = 101
	    normed_train = np.zeros((N,bins))
	    var_train = np.zeros(N)
	    analysis['Xcor'] = np.zeros(numsim)
	    analysis['xhist_Xcor'] = np.linspace(-1,1,bin_hist)
	    analysis['yhist_Xcor'] = np.zeros(bin_hist)
	    analysis['high_Xcor'] = np.zeros((numsim,2))
	    analysis['max_Xcor'] = np.zeros((numsim,2))
	    dt = sim.TC/float(bins)

	  ncid = netcdf.netcdf_file(pathResults + dat,'r',mmap=None)
	  
	  for n in range(net.Const['N']):
	    trainTimeNeuron = ncid.variables['trainTime'][:][np.where(ncid.variables['trainNeuron'][:] == n)[0]]
	    binned_train = np.histogram(trainTimeNeuron,bins=bins,range=[0,sim.TC])[0]
	    normed_train[n] = binned_train - Data['rateNeurons'][idx_file][n]*dt
	    var_train[n] = np.var(binned_train)
	  ncid.close()
	  
	  xcor = np.zeros((N,N))
	  xcor[:] = np.nan

	  for n in range(net.Const['N']):
	    for m in range(n+1,net.Const['N']):
	      xcor[n][m] = np.mean(normed_train[n]*normed_train[m])/np.sqrt(var_train[n]*var_train[m])
	    #print Data['xcor'][idx_file][n]
	    #time.sleep(10)
	  mask = np.invert(np.isnan(xcor))
	  analysis['Xcor'][idx_file] = np.mean(xcor[mask])
	  analysis['yhist_Xcor'] += np.histogram(xcor[mask],bins=bin_hist,range=[-1,1])[0]
	  analysis['max_Xcor'][idx_file] = [np.min(xcor[mask]),np.max(xcor[mask])]
	  analysis['high_Xcor'][idx_file] = [np.sum(xcor[mask] < analysis['max_Xcor'][idx_file][0]/2.),np.sum(xcor[mask] > analysis['max_Xcor'][idx_file][1]/2.)]
	  
	  #print np.mean(xcor[mask]), np.min(xcor[mask]), np.max(xcor[mask])
	  #print np.sum(xcor[mask] < -0.3), np.sum(xcor[mask] > 0.3)
	idx_file += 1
    
    if 'subthreshold' in stats:
      Data['measureTimes'] = Data['measureTimes'][0]
	  
    #print analysis['Xcor']
    #print analysis['yhist_Xcor']
    #print analysis['max_Xcor']
    #print analysis['high_Xcor']
    #plt.figure()
    #plt.plot(analysis['xhist_Xcor'],analysis['yhist_Xcor'])
    #plt.yscale('log')
    #plt.show(block=False)
    
    #### symmetry propto Xcor?
    #### correlated spike trains can be linked to neurons with high out-degree?
    #### bin size changes something?
    #### what would CCG change? should I try?
    
    np.place(Data['rateNeurons'],np.isnan(Data['rateNeurons']),0)	# nan-entries == silent neurons
    
    Data['rateNeurons'] = Data['rateNeurons'].transpose()[:N].transpose()
    # compute statistics...
    #print Data['finalCurrents']
  #print 'key: ' , Data.keys()
  #print Data['inDegree']
  #print Data['inDegree'].shape
  # ... topology ...
  if (('topology' in stats) and (not ('inDegree' in analysis.keys()) or reread)):# and not (net.topo == 'p'):
    #analysis['inDegree'] = np.zeros((numsim,N))
    #analysis['outDegree'] = np.zeros((numsim,N))
    analysis['inDegree'] = Data['inDegree'][:,:N]
    analysis['outDegree'] = Data['outDegree'][:,:N]
    
    Data['rateNetwork'] = np.zeros(numsim)
    Data['inDegreeNet'] = np.zeros(numsim)
    Data['outDegreeNet'] = np.zeros(numsim)
    Data['corr_KNu'] = np.zeros(numsim)
    Data['corr_KK'] = np.zeros(numsim)
    
    for idx_file in range(numsim):
      Data['rateNetwork'][idx_file] = np.mean(Data['rateNeurons'][idx_file])
      Data['inDegreeNet'][idx_file] = np.mean(analysis['inDegree'][idx_file])
      Data['outDegreeNet'][idx_file] = np.mean(analysis['outDegree'][idx_file])
      
      
      # correlation of in-degree and firing rate
      Data['corr_KNu'][idx_file] = np.nanmean((Data['rateNeurons'][idx_file]-Data['rateNetwork'][idx_file])*(analysis['inDegree'][idx_file]-Data['inDegreeNet'][idx_file]))/np.sqrt(np.var(Data['rateNeurons'][idx_file])*np.var(analysis['inDegree'][idx_file]))

      # correlation of in- and out-degree
      Data['corr_KK'][idx_file] = np.nanmean((analysis['outDegree'][idx_file]-Data['outDegreeNet'][idx_file])*(analysis['inDegree'][idx_file]-Data['inDegreeNet'][idx_file])) /np.sqrt(np.var(analysis['outDegree'][idx_file])*np.var(analysis['inDegree'][idx_file]))
    
    analysis['corr_KNu'] = np.nanmean(Data['corr_KNu'])
    analysis['corr_KK'] = np.nanmean(Data['corr_KK'])
    
    analysis['inDegree'] = np.reshape(analysis['inDegree'],N*numsim)
    analysis['outDegree'] = np.reshape(analysis['outDegree'],N*numsim)
    
  
  if (plot and ('topology' in stats)):
    K = net.Const['K']
    y_border = 0.1
    ax_border = [2*K,2*K]
    if max(analysis['inDegree']) > 2*K:
      ax_border[0] = N
    if max(analysis['outDegree']) > 2*K:
      ax_border[1] = N
    
    hist_topo_in = np.histogram(analysis['inDegree'],bins=min(ax_border[0],100),range=[0,ax_border[0]])
    hist_topo_out = np.histogram(analysis['outDegree'],bins=min(ax_border[1],100),range=[0,ax_border[1]])
    
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
    plt.figure(figsize=(2.5,2))
    ax1 = plt.axes([0.05,0.2,0.8,0.75])
    #ax2 = plt.axes([0.55,0.2,0.4,0.75])
    
    
    ax1.plot(hist_topo_in[1][:-1],hist_topo_in[0],color='lightgrey',label=r'$\displaystyle K^{in}$')
    ax1.plot(hist_topo_out[1][:-1],hist_topo_out[0],'--k',label=r'$\displaystyle K^{out}$')
    #axes[0].hist(analysis['inDegree'],bins=min(ax_border[0],100),range=[0,ax_border[0]],color='k')
    #axes[1].hist(analysis['outDegree'],bins=min(ax_border[1],100),range=[0,ax_border[1]],color='k')
    
    y_border = max(np.max(hist_topo_in[0]),np.max(hist_topo_out[0]))
    y_border *= 1.2
    
    ax1.set_xlabel('In-/Out-Degree',fontsize=12)
    #ax2.set_xlabel('Out-degree',fontsize=12)
    #ax1.set_ylabel('Probability',fontsize=12)
    ax1.set_xlim([0,ax_border[0]])
    #ax2.set_xlim([0,ax_border[1]])
    ax1.set_ylim([0,y_border])
    #ax2.set_ylim([0,y_border])
    
    ax1.set_xticks(np.linspace(0,max(ax_border),5))
    #ax1.set_xticklabels(np.linspace(0,ax_border[0],5).astype('int'))
    #ax2.set_xticks(np.linspace(0,ax_border[1],5))
    #ax2.set_xticklabels(np.linspace(0,ax_border[1],5).astype('int'))
    
    ax1.set_yticks([])
    #ax2.set_yticks([])
    
    ax1.yaxis.set_label_position('right')
    ax1.spines['left'].set_color('none')
    ax1.spines['right'].set_color('none')
    ax1.yaxis.set_ticks_position('right')
    ax1.spines['top'].set_color('none')
    ax1.xaxis.set_ticks_position('bottom')
    
    #ax2.yaxis.set_label_position('right')
    #ax2.spines['left'].set_color('none')
    #ax2.spines['right'].set_color('none')
    #ax2.yaxis.set_ticks_position('right')
    #ax2.spines['top'].set_color('none')
    #ax2.xaxis.set_ticks_position('bottom')
    #axes[0].set_yticks(np.linspace(0,y_border*N*numsim,3))
    #axes[0].set_yticklabels(np.linspace(0,y_border,3),fontsize=14)
    #axes[1].set_yticks(np.linspace(0,y_border*N*numsim,3))
    #axes[1].set_yticklabels(np.linspace(0,y_border,3),fontsize=14)
    
    #axes[0].legend(prop={'size':14})
    ax1.legend(prop={'size':12},loc=1)
    #plt.suptitle(r'Erd\H{o}s-R\'{e}nyi network',fontsize=20)
    if save:
      save_str = './pics/%s/statistics/degrees_%s_N=%d.pdf' % (net.scriptName.split('_')[0],'_'.join(str_add),N)
      plt.savefig(save_str)
      print 'Figure saved to ' + save_str
    plt.show(block=False)
  
  #print Data.keys()
  if reread:
    # ... synchrony ...
    analysis['chi'] = np.mean(Data['chi'])
  print 'chi: ', analysis['chi']
  
  # ... moments ...
  if ('moments' in stats) and (not ('rateNeurons' in analysis.keys()) or reread):
    m = 0
    search_range = 5
    #num_bin = 51

    #moment_bars = [np.max(analysis['rateNeurons']),1.5]
    
    num_bin_multiplikator = [5,10]
    
    for key_moment in ['rateNeurons','cvNeurons']:	#iterate over moments
      
      range_max = np.ceil(np.nanmax(Data[key_moment]))
      
      num_bin = int(min(200,range_max*num_bin_multiplikator[m]))
      
      analysis[key_moment+'_hist'] = np.zeros((2,200))	#store moments (moments,axes,bins)
      analysis[key_moment+'_hist'][:] = np.nan
      
      analysis[key_moment+'_peaks'] = np.zeros((2,4))	#store up to 4 peaks (moments,axes,peaks)
      analysis[key_moment+'_peaks'][:] = np.nan
      
      analysis[key_moment] = Data[key_moment].flatten()
      mask = np.invert(np.isnan(analysis[key_moment]))
      # convert to histograms
      tmp_hist = np.histogram(analysis[key_moment][mask],bins=num_bin,range=[0,range_max])
      #print 'second: ', tmp_hist[1][:-1]
      analysis[key_moment+'_hist'][0,:num_bin] = tmp_hist[0]/float(sim.steps['total']*N)
      analysis[key_moment+'_hist'][1,:num_bin] = tmp_hist[1][:-1]
      
      # search for peaks in histograms
      peak_idx = 0
      
      for idx_dat in range(num_bin):
	peak = 1
	moment_search = analysis[key_moment+'_hist'][0][max(0,idx_dat-search_range):min(num_bin,idx_dat+search_range+1)]	
	
	if np.sum(analysis[key_moment+'_hist'][0][idx_dat] <= moment_search)==1:
	  analysis[key_moment+'_peaks'][0][peak_idx] = analysis[key_moment+'_hist'][1][idx_dat]
	  analysis[key_moment+'_peaks'][1][peak_idx] = analysis[key_moment+'_hist'][0][idx_dat]
	  peak_idx += 1
	
	if peak_idx >= 4:
	  break
      
      m += 1
      
      
    
  if (plot and ('moments' in stats)):	# convert data to histograms and find peaks
    m=0
    
    moment_bars = [5*net.Const['rateWnt'],1.6]
    moment_labels = [r'$\displaystyle \nu\,$ in Hz',r'$\displaystyle CV_{ISI}$']
    
    plt.figure(figsize=(2,3))
    ax0 = plt.axes([0.1,0.6,0.85,0.3])
    ax1 = plt.axes([0.1,0.15,0.85,0.3])
    
    mask_hist = np.invert(np.isnan(analysis['rateNeurons_hist'][0]))
    ax0.plot(analysis['rateNeurons_hist'][1,mask_hist],analysis['rateNeurons_hist'][0,mask_hist],'k')
    ax0.plot(analysis['rateNeurons_peaks'][0],analysis['rateNeurons_peaks'][1],'.',color='grey',markersize=10,label='peak')
    ax0.set_xlim([0,moment_bars[0]])
    ax0.set_ylim([0,1.2*np.max(analysis['rateNeurons_hist'][0,mask_hist])])
    #ax1.set_ylabel('Probability',fontsize=12)
    ax0.set_xlabel(moment_labels[0])
    ax0.set_yticks([])
    
    mask_hist = np.invert(np.isnan(analysis['cvNeurons_hist'][0]))
    ax1.plot(analysis['cvNeurons_hist'][1,mask_hist],analysis['cvNeurons_hist'][0,mask_hist],'k')
    ax1.plot(analysis['cvNeurons_peaks'][0],analysis['cvNeurons_peaks'][1],'.',color='grey',markersize=10,label='peak')
    ax1.set_xlim([0,moment_bars[1]])
    ax1.set_ylim([0,1.2*np.max(analysis['cvNeurons_hist'][0,mask_hist])])
    #ax1.set_ylabel('Probability',fontsize=12)
    ax1.set_xlabel(moment_labels[1])
    ax1.set_yticks([])
    
    ax1.set_xticks(np.linspace(0,1.5,4))
    
    ax1.yaxis.set_label_position('right')
    ax1.spines['left'].set_color('none')
    ax1.spines['right'].set_color('none')
    ax1.yaxis.set_ticks_position('right')
    ax1.spines['top'].set_color('none')
    ax1.xaxis.set_ticks_position('bottom')
    
    ax0.spines['left'].set_color('none')
    ax0.spines['right'].set_color('none')
    ax0.yaxis.set_ticks_position('none')
    ax0.spines['top'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.xaxis.set_ticks_position('none')
    if save:
      save_str = './pics/%s/statistics/both_moments%d_%s_N=%d.pdf' % (net.scriptName.split('_')[0],m,'_'.join(str_add),N)
      plt.savefig(save_str)
      print 'Figure saved to ' + save_str
    plt.show(block=False)
      
    #for key_moment in ['rateNeurons','cvNeurons']:
      
      #plt.figure(figsize=(2.5,2))
      #gs = gridspec.GridSpec(3, 1)

      #ax0 = plt.axes([0.2,0.55,0.75,0.35])#subplot(gs[0,0])
      #ax1 = plt.axes([0.2,0.2,0.75,0.35])#subplot(gs[1:3,0])

      #mask_data = np.invert(np.isnan(analysis[key_moment]))
      #mask_hist = np.invert(np.isnan(analysis[key_moment+'_hist'][0]))
      
      #boxprops = dict(linestyle='-', linewidth=1, color='k')
      #flierprops = dict(marker='+', markerfacecolor='k',color='k', markersize=2,linestyle='none')
      #whiskerprops = dict(linestyle='--', linewidth=1, color='k')
      
      #ax0.boxplot(analysis[key_moment][mask_data],notch=1,vert=False,widths=0.5,boxprops=boxprops,flierprops=flierprops,whiskerprops=whiskerprops)
      #ax0.plot(analysis[key_moment+'_peaks'][0],np.ones(4),'.',color='grey',markersize=10,label=moment_labels[m] + 'peak')
      #ax0.set_xlim([0,moment_bars[m]])
      #ax0.set_xticklabels([])
      #ax0.set_yticklabels([])
      
      ##ax1.hist(,bins=100,range=(0,moment_bars[m]),linewidth=1,color='k',histtype='step')
      
      #ax1.plot(analysis[key_moment+'_hist'][1,mask_hist],analysis[key_moment+'_hist'][0,mask_hist],'k')
      #ax1.plot(analysis[key_moment+'_peaks'][0],analysis[key_moment+'_peaks'][1],'.',color='grey',markersize=10,label='peak')
      #ax1.set_xlim([0,moment_bars[m]])
      #ax1.set_ylim([0,1.2*np.max(analysis[key_moment+'_hist'][0,mask_hist])])
      ##ax1.set_ylabel('Probability',fontsize=12)
      #ax1.set_xlabel(moment_labels[m])
      #ax1.set_yticks([])
      #if m == 1:
	#ax1.set_xticks(np.linspace(0,1.5,4))
      
      
      #ax1.yaxis.set_label_position('right')
      #ax1.spines['left'].set_color('none')
      #ax1.spines['right'].set_color('none')
      #ax1.yaxis.set_ticks_position('right')
      #ax1.spines['top'].set_color('none')
      #ax1.xaxis.set_ticks_position('bottom')
      
      #ax0.spines['left'].set_color('none')
      #ax0.spines['right'].set_color('none')
      #ax0.yaxis.set_ticks_position('none')
      #ax0.spines['top'].set_color('none')
      #ax0.spines['bottom'].set_color('none')
      #ax0.xaxis.set_ticks_position('none')
      
      ##ax1.legend(loc=m+1,prop={'size':12},numpoints=1)
      ##plt.suptitle(r'Erd\H{o}s-R\'{e}nyi network',fontsize=20)
      #if save:
	#save_str = './pics/%s/statistics/moments%d_%s_N=%d.pdf' % (net.scriptName.split('_')[0],m,'_'.join(str_add),N)
	#plt.savefig(save_str)
	#print 'Figure saved to ' + save_str

      #plt.show(block=False)
	
      #m += 1
  
  if reread:
    Data['corr_NuCV'] = np.zeros(numsim)
  
    for idx_file in range(numsim):
      mask = np.invert(np.isnan(Data['cvNeurons'][idx_file]))
      Data['corr_NuCV'][idx_file] = np.nanmean((Data['rateNeurons'][idx_file][mask]-np.nanmean(Data['rateNeurons'][idx_file][mask]))*(Data['cvNeurons'][idx_file][mask]-np.nanmean(Data['cvNeurons'][idx_file][mask])))/np.sqrt(np.nanvar(Data['rateNeurons'][idx_file][mask])*np.nanvar(Data['cvNeurons'][idx_file][mask]))
      
    #print Data['corr_NuCV']
  
    analysis['corr_NuCV'] =  np.nanmean(Data['corr_NuCV'])
  
    analysis['finalCurrents'] = np.nanmean(Data['finalCurrents'],axis=0)[0]
    
  # ISI stats of potential decorrelation events
  if ('ISIdecorr' in stats) and (not ('mean_ISIdecorr' in analysis.keys()) or reread):
    # put them all together
    analysis['mean_ISIdecorr'] = np.reshape(Data['mean_ISIdecorr'].transpose(2,0,1),(2,N*numsim))
    analysis['cv_ISIdecorr'] = np.reshape(Data['cv_ISIdecorr'].transpose(2,0,1),(2,N*numsim))
    
  if plot and ('ISIdecorr' in stats):
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
      save_str = './pics/%s/statistics/ISIdecorr_%s_N=%d.pdf' % (net.scriptName.split('_')[0],'_'.join(str_add),N)
      plt.savefig(save_str)
      print 'Figure saved to ' + save_str
    #plt.suptitle(r'ISI of potential decorr. events (%d of %d neurons active)' % (N_eff/numsim,N),fontsize=20)
    plt.show(block=False)
  
  
  #if ('trainCorr' in stats):
    
    #plt.figure()
    #plt.plot(
    
    
    
  # subthreshold statistics: phase distributions
  
  #single trajectories (randomly drawn from first simulation)
  if ('subthreshold' in stats) and reread:
    analysis['measureTimes'] = Data['measureTimes']
    #print analysis['measureTimes']
    analysis['measureStates'] = Data['measureStates']
    
    analysis['phase_moments'] = np.array([np.mean(analysis['measureStates']),np.var(analysis['measureStates'])])

    print "phase mean: %g" % analysis['phase_moments'][0]
    print "phase variance: %g" % analysis['phase_moments'][1]
  
  #if not 'D_spike_fail' in analysis.keys():
  print 'currents: ', analysis['finalCurrents']
  tauM = 0.01
  
  T_free = -tauM*np.log(1.-1./(1+analysis['finalCurrents']))
  print 'T free: ', T_free
  
  PRC = np.zeros(3)
  phi_mean = 0.5
  print 'phases: ', analysis['phase_moments']
  PRC[0] = - tauM/T_free * np.log(np.exp(-analysis['phase_moments'][0]*T_free/tauM) + net.Const['J0']/(net.Const['K']*(1+analysis['finalCurrents']))) - analysis['phase_moments'][0]
  
  phase_low = analysis['phase_moments'][0] - analysis['phase_moments'][1]
  phase_up = analysis['phase_moments'][0] + analysis['phase_moments'][1]
  PRC[1] = - tauM/T_free * np.log(np.exp(-phase_low*T_free/tauM) + net.Const['J0']/(net.Const['K']*(1+analysis['finalCurrents']))) - phase_low
  PRC[2] = - tauM/T_free * np.log(np.exp(-phase_up*T_free/tauM) + net.Const['J0']/(net.Const['K']*(1+analysis['finalCurrents']))) - phase_up
  
  print 'PRC: ', PRC[0]
  print 'average spike failure distance: ', PRC[0]*net.Const['K']
  analysis['D_spike_fail'] = PRC*net.Const['K']
  print analysis['D_spike_fail']
    
  
  if ('subthreshold' in stats) and plot:
    #print Data['measureTimes'].shape
    #print Data['measureStates'].shape
    #print analysis['measureStates'][0,:,np.random.randint(N)].shape
    
    #plt.plot(analysis['measureTimes'],analysis['measureStates'][0,:,np.random.randint(N)])
    #plt.ylim([-1,1])
    #plt.show(block=False)
    
    plt.figure(figsize=(2,1.5))
    ax = plt.axes([0.1,0.3,0.7,0.65])
    # distribution of phases in network (assuming ergodicity, thus distribution is the same at all times)
    #phase_distr, phase_bins = np.histogram(analysis['measureStates'].reshape(np.prod(analysis['measureStates'].shape)),bins=100,range=[-1,1])

    ax.hist(analysis['measureStates'].reshape(np.prod(analysis['measureStates'].shape)),bins=100,range=[-1.5,1],color='k')
    #print "Neuron fraction prior to spiking: %4.3f" % (phase_distr[-1]/float(np.sum(phase_distr))*100)
    ax.set_yticks([])
    ax.set_xlim([-1.5,1.1])
    ax.set_xticks(np.linspace(-1,1,3))
    ax.set_xlabel(r'$\displaystyle \phi$')
    ax.set_ylabel(r'$\displaystyle \rho(\phi)$')
    if save:
      save_str = './pics/%s/statistics/phase_distr_%s_N=%d.pdf' % (net.scriptName.split('_')[0],'_'.join(str_add),N)
      plt.savefig(save_str)
      print 'Figure saved to ' + save_str
      
    plt.show(block=False)
  
  
  # remove old results file
  results_file = net.Path['results'] + 'results_processed.nc'
  if os.path.exists(results_file):
    os.remove(results_file)
  print 'Saving results to %s' % results_file


  # save results in nc-file
  ncid= netcdf.netcdf_file(results_file,'w')
  
  for key in analysis.keys():
    dim_tup = []
    for dim_len in analysis[key].shape:
      dim_name = 'dim_%s_%d' % (key,dim_len)
      if not dim_name in ncid.dimensions:
	ncid.createDimension(dim_name,dim_len)
      dim_tup.append(dim_name)
    if len(dim_tup)==0:
      dim_name = 'one'
      if not 'one' in ncid.dimensions:
	ncid.createDimension('one', 1)
      dim_tup = ('one',)
    elif len(dim_tup)==1:
      dim_tup = (dim_name,)
    else:
      dim_tup = tuple(dim_tup)
    dt = analysis[key].dtype

    NcVar = ncid.createVariable(key,dt.char,dim_tup)
    NcVar[:] = analysis[key]
  ncid.close()
  
  #for key in ['measureStates']:
    #try:
      #trash = analysis.pop(key)
    #except:
      #1
  #print analysis
  
  return analysis
    

def analyze_mosaic(net,sim,plot,save):
  
  ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r',mmap=False)
  n = np.copy(ncid.dimensions['n'])
  pert_size = ncid.variables['pert_size'].getValue()
  fluxtube = ncid.variables['fluxtubes'][:].copy()
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
  
  try:
    net.analysis
  except:
    net.analysis = {}
  
  try:
    print 'start read'
  except:
    1
  results_file = net.Path['results'] + 'results_processed.nc'
  
  if (not os.path.exists(results_file)) or (reread==2 and not os.path.exists(net.Path['results'] + 'clean.nc')):
    print "read it"
    
    net.prob_read(sim)		# get net.analysis
    if not net.break_it:
      save_analysis(net.analysis,results_file)
  
  else: # not len(net.analysis.keys()):	# if stuff should be reread from results file
    ncid = netcdf.netcdf_file(results_file,'r',mmap=False)
    
    # hand over all variables
    for key in ncid.variables.keys():
      #if not (key in ['phi_corr','measureStates']):
      net.analysis[key] = ncid.variables[key][:]
  
  #net.clean() #clean up afterwards
  
  ## implement reading from topo-file here to get connections of spike-crossing neurons
  
  if 'SONET' in net.scriptName:
    para_str = ''
    title_str = r'$'
    for key in net.topoConst.keys():
      if 'alpha' in key:
	para_str += '_%s=%g' % (key,net.topoConst[key])
	str_tmp = key.split('_')
	title_str += r'\displaystyle \%s_{%s}=%g,\;' % (str_tmp[0],str_tmp[1],net.topoConst[key])
    title_str += r'$'
    
    title_paras = ['N','K']
    title_str = ''
    #para_str = 'drive_rate=%g_drive_cplg=%g' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
    #title_str = r'$\displaystyle\bar{\nu}_p=%g, J_p=%g$, ' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
    for para in title_paras:
      title_str += '$%s=%g$, ' %(para,net.Const[para])
      
    if net.Const['special'] == 2:
      fileName = 'get_D_decorr.nc'
      D_get_str = './Code/other/calc_I %s %d %d %g' % (fileName,net.Const['N'],net.Const['K'],net.Const['rateWnt'])
      
      runIt(None,D_get_str)
      
      try:
	ncid = netcdf.netcdf_file(fileName,'r')
	
	net.analysis['D_decorr_analytic'] = ncid.variables['D_decorr'].getValue()
	net.analysis['D_spike_fail'] = ncid.variables['D_spike_failure'].getValue()
	net.analysis['D_spike_fail_reset'] = ncid.variables['D_spike_failure_reset'].getValue()
	ncid.close()
      except:
	net.analysis['D_decorr_analytic'] = np.nan
	net.analysis['D_spike_fail'] = np.nan
	net.analysis['D_spike_fail_reset'] = np.nan
      
      #print "spikefail:\t %g" % net.analysis['D_spike_fail']
      #print "spikefail (reset):\t %g" % net.analysis['D_spike_fail_reset']
      #print "decorr dist: \t %g" % net.analysis['D_decorr_analytic']
    
  else:
    title_paras = ['N','K','J0','rateWnt']
    para_str = 'drive_rate=%g_drive_cplg=%g' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
    title_str = r'$\displaystyle\bar{\nu}_p=%g, J_p=%g$, ' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
    
    for para in title_paras:
      title_str += '$%s=%g$, ' %(para,net.Const[para])
  
  title_str = title_str.rstrip(' ').rstrip(',')
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  mpl.rcParams['xtick.labelsize'] = 12
  mpl.rcParams['ytick.labelsize'] = 12
  mpl.rcParams['font.size'] = 12
  
  if net.Const['special'] == 0:	# large badge of sample to get div/conv statistics
    
    time_len = len(net.analysis['distanceOrth'][0])
    
    if not 'tau_conv_first' in net.analysis.keys() or (reread == 2):
      # get the average convergence characteristic timeconstant
      #print net.analysis['convergeVelocity']
      try:
	net.analysis['tau_conv_first'] = bootstrap(net.analysis['convergeVelocity_1st'],10000)
	net.analysis['tau_conv_later'] = bootstrap(net.analysis['convergeVelocity_2nd'],10000)
      except:
	net.analysis['tau_conv_first'] = bootstrap(net.analysis['convergeVelocity'][0],10000)
	net.analysis['tau_conv_later'] = bootstrap(net.analysis['convergeVelocity'][1],10000)
	trash = net.analysis.pop('convergeVelocity')
      save_analysis(net.analysis,results_file)
    
    #print "tau_conv1: ", net.analysis['tau_conv_first']
    #print "tau_conv2: ", net.analysis['tau_conv_later']
    
    if plot:
      
      plt.figure(figsize=(4,3))
      
      mpl.rcParams['xtick.labelsize'] = 12
      mpl.rcParams['ytick.labelsize'] = 12
      
      if net.topo == 'p' or (net.topoConst['alpha_conv'] > 0):
	
	net.analysis['distanceOrth'][np.where(net.analysis['distanceOrth']==0)] = np.nan 
      
	#print net.analysis['distanceOrth']
	
	if net.topo == 'p':
	  
	  #gs = gridspec.GridSpec(2, 2)
	  ax0 = plt.axes([0.175,0.6,0.7,0.35])
	  ax1 = plt.axes([0.175,0.15,0.35,0.3])
	  #ax2 = plt.axes([0.6,0.7,0.3,0.25])
	  #ax2 = plt.subplot(gs[1,1])
	  ax2 = plt.axes([0.55,0.15,0.325,0.3])
	  
	  decay1 = net.analysis['distanceOrth'][0][0]*np.exp(net.analysis['measureTimes']*net.analysis['tau_conv_first'][0])
	  decay2 = net.analysis['distanceOrth'][0][0]*np.exp(net.analysis['measureTimes']*net.analysis['tau_conv_later'][0])
	  
	  #print net.analysis['convergeVelocity']
	  if not (np.isnan(net.analysis['convergeVelocity_1st'])).all():
	    ax1.plot(range(len(net.analysis['convergeVelocity_1st'][0])),-net.analysis['convergeVelocity_1st'][0],'ko',marker='.',markersize=2,label=r'$\displaystyle \tau_{C}^{(1)}$')
	  ax1.set_ylabel(r'$\displaystyle \tau_{C}$')
	  ax1.set_xlabel('trial \#')
	  ax1.set_xticks(np.linspace(0,sim.steps['total'],3))
	  #ax1.set_ylim([10,100])
	  
	  
	  #print net.analysis['kinkTime']
	  ax2.plot(range(len(net.analysis['kinkTime'])),net.analysis['kinkTime'],'ok',marker='.',markersize=2)
	
	  ax2.set_ylabel('kink time in s',fontsize=10)
	  #ax2.yaxis.set_label_position('right')
	  ax2.yaxis.set_ticks_position('right')
	  ax2.set_xlabel('trial \#')
	  ax2.set_xticks(np.linspace(0,sim.steps['total'],3))
	  ax2.set_ylim([0,0.6])
	  ax2.set_yticks(np.linspace(0,0.6,3))
	  ax2.yaxis.set_label_position('right')
	  #ax2.set_ylabel('later convergence speed')

	  #ax2.yaxis.set_ticks_position('right')  

	else:
	  ax0 = plt.axes([0.2,0.15,0.75,0.75])
	  ax1 = plt.axes([.29, .25, .3, .3])
	  #ax0.plot(net.analysis['measureTimes'],decay1,'--k',linewidth=2,label='approximate initial convergence')
	  #ax0.plot(net.analysis['measureTimes'],decay2,'--r',linewidth=5,label='approximate later convergence')
	  ax1.set_ylabel(r'$\displaystyle \tau_{C}$',fontsize=14)

	
	
	ax1.plot(range(len(net.analysis['convergeVelocity_2nd'][0])),-net.analysis['convergeVelocity_2nd'][0],'o',marker='.',color='grey',markersize=2,label=r'$\displaystyle \tau_{C}^{(2)}$')
	#ax2.set_xlabel('trial \#')
	#ax2.set_xlim([0,sim.steps['total']])
	ax1.set_yscale('log')
	ax1.set_ylim([10**(-2.5),10**0])
	#ax2.set_yscale('log')
	#ax2.set_xticklabels([])

      else:
	ax0 = plt.axes([0.2,0.15,0.75,0.75])
      
      
      #ax2.set_yscale('log')
      #ax2.set_xticks(np.linspace(0,0.1,3))
      #ax2.set_xticklabels(np.linspace(0,0.1,3),fontsize=12)
      #ax2.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '{:2g}'.format(y)))
      
      
      #ax2.set_xlim([0,0.1])
      #ax2.set_ylim([10**(-4),10**(1.5)])
      #ax2.set_yticks(10**np.linspace(-3,1,3))
      #ax2.set_yticklabels(10**np.linspace(-3,1,3),fontsize=12)
      ##ax0.xaxis.get_major_formatter().set_powerlimits((0,1))
      #ax2.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
      
      
      #x_array = net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])]
      mask = np.where(net.analysis['distanceOrth'].transpose()[-1] < 0.01)[0]
      
      div = np.zeros(sim.steps['total']).astype('bool')#net.analysis['distanceOrth'][:,time_len-1] > 0.1
      
      for n in range(sim.steps['total']):
	mask = np.invert(np.isnan(net.analysis['distanceOrth'][n]))
	if not any(mask):
	  div[n] = False
	else:
	  div[n] = net.analysis['distanceOrth'][n][mask][-1] > 0.1
            
      if any(div):
	ax0.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][div][:100,:time_len].transpose(),'-',color='orangered',linewidth=0.5)
      if not all(div):
	ax0.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][np.invert(div)][:100,:time_len].transpose(),'-',color='dodgerblue',linewidth=0.5)
      
      #if any(div):
	#ax2.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][div][:10,:time_len].transpose(),'-',color='orangered',linewidth=0.5)
      #if not all(div):
	#ax2.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][np.invert(div)][:10,:time_len].transpose(),'-',color='dodgerblue',linewidth=0.5)


      decay1 = net.analysis['distanceOrth'][0][0]*np.exp(net.analysis['measureTimes']/net.analysis['tau_conv_first'][0])
      
      decay2 = net.analysis['distanceOrth'][0][0]*np.exp(net.analysis['measureTimes']/net.analysis['tau_conv_later'][0])
      #print np.nanmin(net.analysis['distanceOrth'][np.invert(div)])
      
      ax0.plot(net.analysis['measureTimes'][:time_len],decay1,'--k',linewidth=2,label=r'$\displaystyle \tau_{C}^{(1)}$')
      ax0.plot(net.analysis['measureTimes'][:time_len],decay2,'--',color='grey',linewidth=2,label=r'$\displaystyle \tau_{C}^{(2)}$')
      
      mpl.rcParams['xtick.labelsize'] = 12
      mpl.rcParams['ytick.labelsize'] = 12
      ax0.set_xlim([0,1])
      #ax0.set_xlim([0,sim.TC])
      ax0.set_yscale('log')
      
      ax0.set_ylim([10**(-11),10**1.5])
      ax0.set_yticks(10**np.linspace(-11,1,4))
      ax0.set_yticklabels(10**np.linspace(-11,1,4),fontsize=12)
      #ax0.xaxis.get_major_formatter().set_powerlimits((0,1))
      #ax0.set_ylim([10**(-8),np.ceil(net.ParaSim['D_decorr'])*2])
      ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
      
      
      ax0.set_xlabel('time in s',fontsize=12)
      
      ax0.set_ylabel(r'$\displaystyle D_{\phi}$',fontsize=14)
      
      # here comes comparison to convergence from lyapunov exponent:
      #lyap_decay = net.analysis['distanceOrth'][0][0]*np.exp(-62.6*x_array)
      #ax0.plot(x_array,lyap_decay,'--k',linewidth=1)#,label='approximation from 1st LE')
      
      #ax0.plot([0,sim.TC],[net.ParaSim['D_decorr']]*2,'--r',linewidth=2)#,label=r'$\displaystyle D_{\phi}^{decorr}$')
      #ax0.plot(x_array,np.ones(len(x_array))*net.ParaSim['D_spike_failure'],'--r',linewidth=1,label=r'$\displaystyle D_{\phi}^{spike failure}$')

      
      #plt.suptitle(title_str,fontsize=14,y=0.97)
      ax0.legend(loc=1,prop={'size':10})
      #ax1.legend(loc=3,prop={'size':10})
      if save:
	save_str = './pics/%s/samples/huge_batch_div_N=%d_K=%d_f=%3.1f%s.pdf' % (net.scriptName.split('_')[0],net.Const['N'],net.Const['K'],net.Const['rateWnt'],para_str)
	plt.savefig(save_str)
	print "Figure saved to %s." % save_str
	
      plt.show(block=False)
      #if len(mask) < 10:
	#for i in range(net.analysis['distanceOrth_max'].shape[0]):
	  #if net.analysis['distanceOrth'][i][-1] < 0.01:
	    #plt.figure()
	    #plt.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])],net.analysis['distanceOrth'][i],'-b',linewidth=0.1)	# orthogonalized distance
	    #if net.topoConst['alpha_conv'] > 0:
	      ##print net.analysis['distanceOrth_min'][i]
	      #plt.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth_min'][0])],net.analysis['distanceOrth_min'][i],'--r',linewidth=2,label=r'10\% lowest $\displaystyle K^{in}$')	# orthogonalized distance
	      #plt.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth_max'][0])],net.analysis['distanceOrth_max'][i],'--g',linewidth=2,label=r'10\% highest $\displaystyle K^{in}$')	# orthogonalized distance
	    #plt.ylim([10**(-12),5])
	    #plt.yscale('log')
	    #plt.legend()
	    
	  #plt.show(block=False)
    
    
  if net.Const['special'] == 1:
    
    #print net.analysis.keys()
    if plot:
      mpl.rcParams['xtick.labelsize'] = 12
      mpl.rcParams['ytick.labelsize'] = 12
      
      time_len = len(net.analysis['distanceOrth'][0])
      
      plt.figure(figsize=(3,1.9))
      ax0 = plt.axes([0.25,0.22,0.7,0.49])
      #ax1 = plt.axes([0.6,0.6,0.3,0.3])
      
      div = net.analysis['distanceOrth'][:,time_len-1] > 0.1
      
      
      
      
      try:
	ax0.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_ER'][0],':k',linewidth=2)#,label=r'$\displaystyle D_{\phi}^{(dc)}$')
      except:
	1
      #ax0.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*np.sqrt(2),'--k',linewidth=2,label=r'$\displaystyle D_{\phi}^{(sf)}$')
      #ax0.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len),'--k',linewidth=2)
      #print np.argwhere(net.analysis['distanceOrth'][np.invert(div)][0]<10**(-4))[0][0]
      
      if any(div):
	ax0.plot(net.analysis['measureTimes'][:time_len],sp.signal.medfilt(net.analysis['distanceOrth'][div,:time_len][:10],1).transpose(),'-',color='orangered',linewidth=0.5)
      if not all(div):
	ax0.plot(net.analysis['measureTimes'][:time_len],sp.signal.medfilt(net.analysis['distance'][np.invert(div),:time_len][:10],3).transpose(),'-',color='dodgerblue',linewidth=0.5)
	try:
	  time_border = net.analysis['measureTimes'][np.argwhere(net.analysis['distanceOrth'][np.invert(div)][0]<10**(-3.5))[0]]
	except:
	  time_border = 0.1
      
      if 'tau_conv_ER' in net.analysis.keys():
	ER_decay = np.exp(net.analysis['measureTimes'][:time_len]/net.analysis['tau_conv_ER'][1])*net.analysis['distanceOrth'][0,0]
	ax0.plot(net.analysis['measureTimes'][:time_len],ER_decay,'--k',linewidth=2,label=r'$\displaystyle \tau_C^{ER}$')
	
      #if '1stLE' in net.analysis.keys():
	#if net.topo == 'p':
	  #LE_decay = np.exp(net.analysis['measureTimes'][:time_len]*net.analysis['1stLE'][1])*net.analysis['distanceOrth'][0,0]
	  #LE_decay0 = np.exp(net.analysis['measureTimes'][:time_len]*net.analysis['1stLE'][0])*net.analysis['distanceOrth'][0,0]*0.5
	  #ax0.plot(net.analysis['measureTimes'][:time_len],LE_decay0,'--',color='grey',linewidth=2,label='0th LE')
	#else:
	  #LE_decay = np.exp(net.analysis['measureTimes'][:time_len]*net.analysis['1stLE'])*net.analysis['distanceOrth'][0,0]
	
	#ax0.plot(net.analysis['measureTimes'][:time_len],LE_decay,'--k',linewidth=2,label='1st LE')
	
      
      #ax1.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_ER'][0],'--r',linewidth=2)
            
      #if any(div):
	#ax1.plot(net.analysis['measureTimes'][:time_len],distance_med[div,:time_len].transpose(),'-',color='orangered',linewidth=0.5)
      #if not all(div):
	#ax1.plot(net.analysis['measureTimes'][:time_len],distance_med[np.invert(div),:time_len].transpose(),'-',color='dodgerblue',linewidth=0.5)

      #for i in range(sim.steps['total']):
	#mask = np.where(net.analysis['distanceOrth'][i] > 0)
	#if net.analysis['distanceOrth'][i][mask][-1] > 0.1:
	  #col = 'orangered'
	#else:
	  #col = 'dodgerblue'
	#print net.analysis['distanceOrth'][i]-net.analysis['distanceOrth'][i]
	#ax0.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])],net.analysis['distanceOrth'][i],'-',color=col)	# orthogonalized distance
	#ax0.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][i],'-',color=col)	# orthogonalized distance

      #ax0.set_xlim([0,0.2])#sim.TC])
      
#<<<<<<< HEAD
      plt.setp(ax0,xlim=[0,0.5],xticks=np.linspace(0,0.4,5))
      plt.setp(ax0, yscale='log',ylim=[10**(-5.5),10**1.5],ylabel=r'$\displaystyle D_{\phi}$')
      ax0.set_yticks(10**np.linspace(-5,1,4))
#=======
      #plt.setp(ax0,xlim=[0,0.4],xticks=np.linspace(0,0.4,5))
      #plt.setp(ax0, yscale='log',ylim=[10**(-7.5),10**1.5],ylabel=r'$\displaystyle D_{\phi}$')
      #ax0.set_yticks(10**np.linspace(-7,1,4))
#>>>>>>> 9d8162143eab98042456584a6e59f456d8dd7819

      ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))

      
      ax0.set_xlabel('t in s',fontsize=12)
      
      if net.topo =='p':
	text_box = r'$\displaystyle \nu_p=%5.1f\,$Hz, $\displaystyle J_p = %d$' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
	#print text_box
	l = ax0.legend(title=text_box,ncol=2,bbox_to_anchor=(0.5,1.05),loc='lower center',borderaxespad=0,prop={'size':10})
	#plt.setp(l.get_title(),fontsize=10)
      else:
	ax0.legend(loc=1,prop={'size':10})
      
      if net.topo == 'S':
	text_box = ''
	for motif in net.topoConst.keys():
	  if net.topoConst[motif]:
	    text_box = motif.split('_')[1]
	
	    ax0.text(0.025, 10**(-6), r'$\displaystyle \alpha_{%s}=%g$'%(text_box,net.topoConst[motif]), bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=14)
	    
	    if 'chain' in motif:
	      break
      #else:
	#text_box = r'$\displaystyle \nu_p=%5.1f\,$Hz, $\displaystyle J_p = %d$' % (net.topoConst['drive_rate'],net.topoConst['drive_cplg'])
	#ax0.text(0.05, 10**(1.4), text_box, bbox={'facecolor':'white','alpha':1,'pad':5},fontsize=12)
      
      if save:
	save_str = './pics/%s/samples/small_batch_N=%d_K=%d_f=%3.1f%s.pdf' % (net.scriptName.split('_')[0],net.Const['N'],net.Const['K'],net.Const['rateWnt'],para_str)
	plt.savefig(save_str)
	print "Figure saved to %s." % save_str
      
      plt.show(block=False)
      
      
      
      if plot == 2:
	#print net.analysis['distanceOrth']
	a = 1
	i = np.where(net.analysis['distanceOrth'][:,-1] < 0.1)[0][a]
	#i = np.where(net.analysis['distanceOrth'][:,-1] < 0.1)[0][np.random.randint(18)]
	
	print i
	#i=6
	smallest = heapq.nsmallest(int(0.005*net.Const['N']),net.analysis['inDegree'][i])
	
	smallest = np.unique(smallest)
	
	
	
	fig = plt.figure(figsize=(3,3))
	
	ax = plt.axes([0.1,0.625,0.6,0.35])
	ax_time = plt.axes([0.2,0.15,0.75,0.3])
	
	ax_time.plot(net.analysis['measureTimes'],sp.signal.medfilt(net.analysis['distanceOrth'][i],3),'-',color='dodgerblue',label='orth. distance')
	ax_time.plot(net.analysis['measureTimes'],sp.signal.medfilt(net.analysis['distance'][i],3),'--',color='dodgerblue',label='distance')
	
	ax_time.set_yscale('log')
	ax_time.set_xlabel('t in s')
	
	plt.ion()
	plt.show(block=False)
	
	bin_pert = 201
	pert_hist = np.zeros((time_len,bin_pert))
	borders = 1*10**(-5)
	#pert_vect = (net.analysis['measureStatesRef']-net.analysis['measureStates'][i])
	
	#distance = np.sqrt(np.sum(pert_vect**2,axis=1))
	#if 'LyapunovExponents' in net.analysis.keys():
	  #ax_time.plot(net.analysis['measureTimes'],decay_max,'--r',label='maximum LE')
	  #ax_time.plot(net.analysis['measureTimes'],decay_mean,'--y',label='mean LE')
	ax_time.legend(bbox_to_anchor=(0.4,0.6),loc='lower left',borderaxespad=0,prop={'size':10})
	#ax_time.set_xlabel('time in s',fontsize=12)
	ax_time.set_ylabel('$\displaystyle D_{\phi}$',fontsize=14)
		
	ax_time.set_xticks(np.linspace(0,sim.TC,11))
	#ax_time.set_xticklabels(np.linspace(0,sim.TC,6))
	
	ax_time.set_yticks(10**np.linspace(-9,1,6))
	#ax_time.set_yticklabels(10**np.linspace(-9,1,6))
	
	ax_time.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
	mpl.rcParams['xtick.labelsize'] = 12
	mpl.rcParams['ytick.labelsize'] = 12
      
	for t in [50,80,100,200]:
	#for t in [300]*10:
	#low_idx = ax.plot(0,0)
	#for t in range(time_len):
	  ax_time.set_ylim([10**(-8),10**(0.5)])
	  ax_time.set_xlim([0,0.25])
	  #ax.cla()
	  #print np.max(t)
	  col = 4*(200-t)/float(5*200)
	  #col = 0
	  #print col
	  #pert_vect[t] /= net.analysis['distanceOrth'][0][t]
	  
	  max_hist = max(abs(min(net.analysis['distanceVector'][i][t])),max(net.analysis['distanceVector'][i][t]))
	  #if max_hist > borders:
	    #borders *= np.sqrt(10)
	  #elif max_hist < borders/5:
	    #borders /= np.sqrt(10)
	  #pert_hist[t] = np.histogram(pert_vect[t],bins=bin_pert,range=[-borders,borders])[0]
	  pert_hist[t] = np.histogram(net.analysis['distanceVector'][i][t],bins=bin_pert,range=[-borders,borders])[0]
	  ax.plot(np.linspace(-borders,borders,bin_pert),pert_hist[t],label=r'%g' % net.analysis['measureTimes'][t],color=(col,col,col))
	  ax.plot([0,0],[0,1000],'--',c='g')
	  time_ax = ax_time.plot(net.analysis['measureTimes'][t],net.analysis['distanceOrth'][i][t],c=(col,col,col),marker='.',markersize=10)
	  
	  #if not (net.topo == 'p'):
	    #print "he"
	    #for u in smallest:
	      #idx_u = np.where(net.analysis['inDegree'][i]==u)[0]
	      #low_idx = ax.plot(net.analysis['distanceVector'][i,t,idx_u],np.zeros(len(idx_u)),'ok')
	    
	  #ax.set_yscale('log')
	  ax.legend(title='t in s',bbox_to_anchor=(1.05,0),loc='lower left',borderaxespad=0,prop={'size':10})
	  ax.set_xlabel(r'$\displaystyle \Delta \phi$',fontsize=14)
	  ax.set_ylabel(r'$\displaystyle \rho(\Delta \phi)$',fontsize=14)
	  
	  ax.set_xticks(np.linspace(-borders,borders,3))
	  #ax.set_xticklabels(np.linspace(-borders,borders,5))
	  
	  ax.set_yticks([])
	  #ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '{:.1g}'.format(y/(2*net.Const['N']))))
	  #ax.set_yticklabels(np.linspace(0,200./(2*net.Const['N']),3))
	  ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '{:.3g}'.format(y)))
	  #time.sleep(0.5)
	  ax.set_xlim(-borders/10.,borders)
	  ax.set_ylim(0,net.Const['N']*0.1)
	  
	  fig.canvas.draw()
	    
	  #time_ax.remove()
	#if save:
	  #save_str = './pics/%s/samples/slowed_conv_expl.pdf' % (net.scriptName.split('_')[0])
	  #plt.savefig(save_str)
	  #print "Figure saved to %s." % save_str
      #print np.sum(net.analysis['distanceOrth'].transpose()[-1] > 0.1)
      #cbar_steps = 5
      #levs = range(160)
      #rb = mcolors.LinearSegmentedColormap.from_list(name='red_black',colors=[(0,(1,1,1)),(1,(1,0,0))],N=len(levs)-1,)
      #cbar_dummy = dummy_plot(rb,cbar_steps,cbar_steps+1)
      #plt.close('all')
      #i = np.where(net.analysis['distanceOrth'].transpose()[-1] < 0.01)[0][0]
      #print net.analysis['distanceOrth']
      
      ##for i in range(len(net.analysis['distanceOrth'])):
	##mask = np.where(net.analysis['trainDeltaTime'][i] < 1)[0]
	
	###print net.analysis['trainDeltaTime'][i]
	###col = np.log10(net.analysis['trainDeltaTime'][i][:len(net.analysis['trainTimeRef'])])
	##col = net.analysis['trainDeltaTime'][i][:len(net.analysis['trainTimeRef'])]
	
	##print net.analysis['distanceOrth'].transpose()[-1][i]
	
	##max_col = np.ceil(max(col[mask])*100)/100. # round to next 10ms
	###print max_col
	
	##down_scale = 10.
	##plt.figure(figsize=(4,2))
	##ax = plt.axes([0.2,0.25,0.7,0.7])
	##ax.scatter(net.analysis['trainTimeRef'],net.analysis['trainNeuronRef'],c=col,cmap=rb,vmin=0,vmax=max_col/down_scale,marker='D')
	###ax1.scatter(net.analysis['trainTime'][i],net.analysis['trainNeuron'][i],color='k')
	##ax.set_xlabel('time in s',fontsize=22)
	##ax.set_ylabel('Neuron \#',fontsize=22)
	##ax.set_xlim([0,0.1])
	##ax.set_ylim([0,300])
	
	##print 1000*np.linspace(0,max_col/down_scale,int(cbar_steps/2)+1)
	###print np.linspace(0,max_col/10.,cbar_steps)
	##cbar0 = plt.colorbar(cbar_dummy)
	##cbar0.ax.set_ylabel(r'$\displaystyle \Delta t$ in ms',fontsize=22)
	##cbar0.set_ticks(np.linspace(1./2,cbar_steps-1./2,int(cbar_steps/2)+1))
	##cbar0.set_ticklabels(['%5.3g'%i for i in 1000*np.linspace(0,max_col/down_scale,int(cbar_steps/2)+1)])
  
	
	##plt.show(block=False)
	  
  if net.Const['special'] == 2:	# get decorr distance
    
    #print net.analysis['D_decorr']
    if not ('D_decorr_bs' in net.analysis.keys()) or (reread == 2):
      net.analysis['D_decorr_bs'] = bootstrap(net.analysis['D_decorr'],10000)
      save_analysis(net.analysis,results_file)
    
    #print net.analysis['D_decorr_bs']
    #print net.analysis['D_decorr'].shape
    #print net.analysis.keys()
    
    print "(Bootstrapping) D_decorr: %g, lower: %g, upper %g" % (net.analysis['D_decorr_bs'][0],net.analysis['D_decorr_bs'][1],net.analysis['D_decorr_bs'][2])
  
    
  
  
  
  if net.Const['special'] in [3,30]:
    
    #print net.analysis['D_decorr']
    #print net.analysis['tau_conv']
    #print net.Hash
    if 'D_decorr' in net.analysis.keys():
      if np.isnan(net.analysis['D_decorr'][0]):
	net.analysis['D_decorr'][0] = np.array([1,np.nan,np.nan])
	print "D decorr set to minimum value 1 (since simulation could not find D)"
    else:
      net.analysis['D_decorr'] = np.array([1,np.nan,np.nan])
      print "D decorr set to minimum value 1 (since no simulation for D found)"
    #if np.isnan(net.analysis['tau_conv'][0]):
      #net.analysis['tau_conv'][0] = -0.01		# fastest speed possible
      #print "convergence speed set to maximum value 0.01"
    
    mask = np.invert(np.isnan(net.analysis['T_conv']))
    
    n_finite = np.sum(np.isfinite(net.analysis['T_conv']))
    print 'n finite: ', n_finite
    
    #T_conv_mask = net.analysis['T_conv'][mask]
    #print net.analysis['T_conv'].shape
    #mask = np.invert(np.isnan(net.analysis['T_conv']))
    
    #print np.nanmean(net.analysis['T_conv'])
    print net.analysis['T_conv'].shape
    #print net.analysis['T_conv']
    
    #print np.sum(mask)
    steps = sim.steps['co']*sim.steps['in']*np.sum(np.arange(sim.steps['pd']))
    print "number of values: %d/%d " % (n_finite,np.sum(mask))
    if np.sum(mask):
      
      T_hist = np.histogram(net.analysis['T_conv'][mask],bins=np.append(0,net.analysis['measureTimes']),range=[0,net.analysis['measureTimes'][-1]])
      T_hist_vals = 1 - np.cumsum(T_hist[0]/float(np.sum(mask)))
      
      net.analysis['T_min'] = np.nanmin(net.analysis['T_conv'])
      print 'minimum value of T: ', np.nanmin(net.analysis['T_conv'])
      
      try:
	idx_start = np.where(net.analysis['measureTimes'] >= net.analysis['T_min'])[0][0]
	
	try:
	  idx_end = np.where(T_hist_vals < 0.1)[0][0]
	except:
	  idx_end = -1
	print T_hist_vals[idx_start:idx_end]
	
	#popt, pcov = curve_fit(func_exp,net.analysis['measureTimes'][idx_start:idx_end],T_hist_vals[idx_start:idx_end])

	popt, pcov = curve_fit(func_linear,net.analysis['measureTimes'][idx_start:idx_end]+np.nanmin(net.analysis['T_conv']),np.log(T_hist_vals[idx_start:idx_end]),sigma=T_hist_vals[idx_start:idx_end])
      
	net.analysis['tau_RC'] = (-1./popt[0])
	print "characteristic time: %g" % (net.analysis['tau_RC'])
	
      except:
	net.analysis['tau_RC'] = np.inf
	print "characteristic time: %g" % (net.analysis['tau_RC'])
      
      if plot==1:
	try:
	  plt.figure(figsize=(2,2))
	  ax0 = plt.axes([0.34,0.25,0.63,0.7])
	  
	  idx_drop=np.where(T_hist_vals<np.exp(-1))[0][0]
	  ax0.plot(net.analysis['measureTimes'],T_hist_vals,'k')
	  ax0.plot(net.analysis['measureTimes'],np.exp(popt[0]*(net.analysis['measureTimes']-np.nanmin(net.analysis['T_conv']))),'--r')
	  ax0.set_xticks(np.linspace(0,1,6))
	  ax0.set_yscale('log')
	  ax0.set_xlim([0,0.5])
	  ax0.set_ylim([10**(-2),1])
	  
	  plt.annotate(r'$\displaystyle \tau_{RC}$',xy=[net.analysis['measureTimes'][idx_drop],np.exp(-1)],xytext=[net.analysis['measureTimes'][idx_drop]*0.5,np.exp(-1)*0.2],arrowprops=dict(arrowstyle="->"),fontsize=12)
	  
	  ax0.set_xlabel('t in s')
	  ax0.set_ylabel('$\displaystyle f_{D_{\phi}>\Theta}$',fontsize=14)
	  if save:
	    save_str = './pics/%s/samples/fraction_div_nup=%d.pdf' % (net.scriptName.split('_')[0],int(net.topoConst['drive_rate']))
	    plt.savefig(save_str)
	    print "Figure saved to %s." % save_str

	  plt.show(block=False)
	except:
	  1
	
	#plt.figure()
	#plt.hist(net.analysis['T_conv'][mask],bins=1000)
	#plt.show()
	
	
      if plot == 2:
	thresh_drop = 10**(-4)
	dt = -np.log(thresh_drop/net.analysis['D_decorr'][0])*0.01
	print net.analysis['D_decorr'][0]
	print 'dt: ', dt
	
	print net.analysis.keys()
	for co_ct in range(sim.steps['co']):
	  for in_ct in range(sim.steps['in']):
	    #for i in range(sim.steps['pd']):
	      #print net.analysis['T_conv'][co_ct,in_ct,i]
	    
	    T = net.analysis['T_conv_measure'][co_ct,in_ct]
	    
	    mask = np.invert(np.isnan(T))
	    event_times = np.append(10**(-2),np.unique(T[mask]))
	    event_times = np.append(event_times,0.5)
	    
	    if not len(T[mask]):
	      print 'skip'
	      continue
	    
	    # initiate matrices
	    parent_idx = 0
	    parent_matrix = np.zeros((len(event_times),sim.steps['pd']))
	    
	    parent_slice = np.zeros(sim.steps['pd'])
	    parent_slice[:] = np.nan
	    num_ft = np.zeros(len(event_times))
	    num_ft[:] = np.nan
	    #parent_matrix = np.zeros((len(times),sim.steps['pd']))
	    
	    # assign the final states:
	    parent_slice[0] = parent_idx
	    parent_slice[np.invert(np.isnan(T[0]))] = parent_idx	# for every entry in T[0], one trajectory collapsed onto it 
	    parent_idx += 1
	    
	    # now assign not yet collapsed trajectories
	    for n in range(sim.steps['pd']):
	      if np.isnan(parent_slice[n]):
		parent_slice[n] = parent_idx
		parent_slice[np.invert(np.isnan(T[n]))] = parent_idx
		parent_idx += 1
	    
	    #print parent_slice
	    
	    for t_fill in range(len(event_times)):
	    #for t_fill in range(len(times)):
	      parent_matrix[t_fill] = parent_slice
	      
	    # now, assign numbers to the converging trajectories:
	    #print 'measures: ', net.analysis['measureTimes']
	    print 'events: ', event_times
	    
	    t_list = np.zeros(len(event_times))
	    for i in range(len(event_times)):
	      t_list[i] = net.analysis['measureTimes'][np.argmin(abs(net.analysis['measureTimes'] - event_times[i]))]
	    
	    #print 'real events: ', np.unique(t_list)
	    for t in reversed(event_times):
	      
	      idx_time = np.argwhere(t==event_times)[0,0]
	      #print idx_time, event_times[idx_time]
	      #if idx_time <= 10:
		#time.sleep(5)
	      # re-initiate parent slice
	      parent_slice = np.zeros(sim.steps['pd'])
	      parent_slice[:] = np.nan
	      parent_idx = 0
	      
	      conv_tuples = np.where(T == t)
	      #print 'current assignment: ', parent_matrix[idx_time]
	      #print conv_tuples
	      
	      #if len(conv_tuples[0]) > 1:
		#print 'whaddup?'
	      
	      #for j in range(1):
		#for i in range(len(conv_tuples[0])):
	      #n = conv_tuples[0][0]
	      #print 'splitting: ', n
	      
	      #p = parent_matrix[idx_time][n]
	      #print 'current population:', p

	      #mask_add = parent_matrix[idx_time] > p
	      
	      #parent_slice[np.invert(mask_add)] = parent_matrix[idx_time][np.invert(mask_add)]
	      #parent_slice[mask_add] = parent_matrix[idx_time][mask_add] + 1
	      #print 'new slice 1: ', parent_slice
	      
	      
	      #parent_slice[n] = p + 1
	      #print 'times: ', T[n], T[:,n]
	      #print 'times2: ', T[conv_tuples[1][0]], T[:,conv_tuples[1][0]]
	      #parent_slice[np.where(T[n]<t)[0]] = p + 1
	      #parent_slice[np.where(T[:,n]<t)[0]] = p + 1
	      
	      #parent_idx += 1
	      
	      
	      #print 'new slice 2: ', parent_slice
	      # iterate through populations in current parent-structure
	      for p in np.unique(parent_matrix[idx_time]):
		
		idx_set = np.where(parent_matrix[idx_time] == p)[0]	# this is the currently evaluated population
		
		for j in range(2):
		  for i in range(len(conv_tuples[0])):	# should be <=1, as no synchronous spikes present
		    n = conv_tuples[j][i]
		    if n in idx_set:
		      if np.isnan(parent_slice[n]):
			parent_slice[n] = p + parent_idx
			parent_idx += 1
			
		      p_new = parent_slice[n]
		      
		      update_parent_slice(T,parent_slice,t,n,p_new)
		      
		      parent_slice[np.where(T[n]<t)[0]] = p_new
		      parent_slice[np.where(T[:,n]<t)[0]] = p_new
		      
		if not all(np.isnan(parent_slice[idx_set])):
		  parent_idx -= 1
			  
		#iterate non-diverging trajectories
		for n in idx_set:
		  if np.isnan(parent_slice[n]):
		    parent_slice[n] = p + parent_idx
		    
	      for n in range(sim.steps['pd']):
		indices_update = np.append(np.where(T[:,n]<t)[0],np.where(T[n]<t)[0])
		if len(indices_update):
		  if not all(parent_slice[indices_update] == parent_slice[n]):
		    parent_slice[indices_update] = parent_slice[n]
	      
	      num_entries = len(np.unique(parent_slice))
	      reorder = np.vstack((np.arange(num_entries),np.unique(parent_slice))).transpose()
	      for n in range(sim.steps['pd']):
		idx_reorder = np.where(parent_slice[n]==reorder[:,1])[0][0]
		parent_slice[n] = reorder[idx_reorder,0]
	      
	      #print parent_slice
	      #print len(np.unique(parent_slice))/float(sim.steps['pd'])
	      for t_fill in range(max(0,idx_time - 2),idx_time):
		num_ft[t_fill] = len(np.unique(parent_slice))/float(sim.steps['pd'])
		if np.isnan(num_ft[t_fill]):
		  print parent_slice
		parent_matrix[t_fill] = parent_slice
		
	      #print parent_matrix[t]
	    
	    for t in range(len(event_times)):
	      #print len(set(parent_slice))
	      
	      #plt_array = np.linspace(0,1,np.max(parent_matrix[t]) + 3)
	      plt_array = np.linspace(0,sim.steps['pd']+1,np.max(parent_matrix[t]) + 3)
	      parent_matrix[t] = plt_array[parent_matrix[t].astype('int')+1]
	    #print np.max(parent_matrix)
	    #parent_matrix[t] = parent_matrix[t].astype('float')/np.max(parent_matrix[t])*(sim.steps['pd']-1)
	  
	    #print parent_matrix
	    
	    plt.figure(figsize=(3,2.5))
	    ax1 = plt.axes([0.2,0.2,0.75,0.75])
	    #ax2 = plt.axes([0.15,0.15,0.8,0.25])
	    
	    #t_max = np.nanmax(times)
	    for n in range(sim.steps['pd']):
	      #print parent_matrix[:,n]
	      #ax1.plot(net.analysis['measureTimes'],parent_matrix[:,n],'k')
	      ax1.plot(event_times-dt,parent_matrix[:,n],'k')
	    ax1.set_xscale('log')
	    ax1.set_xlabel('time in s')
	    ax1.set_xlim([10**(-2),0.4])
	    
	    #ax1.set_xlim([10**(-2),np.nanmax(net.analysis['measureTimes'])+1])
	    ax1.set_ylim([0,sim.steps['pd']+1])
	    ax1.set_ylabel(r'\# simulation',fontsize=12)
	    #print num_ft
	    #print net.analysis['measureTimes']
	    #ax2.plot(event_times-dt,num_ft)
	    #ax2.set_xscale('log')
	    ##ax2.set_xlim([10**(-2),np.nanmax(net.analysis['measureTimes'])+1])
	    #ax2.set_xlim([10**(-2),0.4])
	    #ax2.set_ylim([0,1])
	    #ax2.set_xlabel('time in s',fontsize=12)
	    #ax2.set_ylabel(r'fraction diverged',fontsize=12)
	    
	    plt.show(block=False)
    
    trash = net.analysis.pop('D_decorr')
    
  if net.Const['special'] == 4:
    
    time_len = len(net.analysis['distanceOrth'][0])-1
    div = net.analysis['distanceOrth'][:,time_len-1] > 0.1
    #print np.where(net.analysis['distanceOrth'][:,time_len-1] > 0.1)[0]
    try:
      i = np.where(net.analysis['distanceOrth'][:,time_len-1] > 0.1)[0][np.random.randint(400)]
    except:
      i = 0
    #i = div[np.random.randint(50)]
    #print i
    #print net.analysis['distanceOrth'][i,:time_len-1]
    mask = np.isfinite(net.analysis['trainDeltaT'][div,:,2])
    #if not 'assigned' in net.analysis.keys():
    n_spikes = np.sum(np.isfinite(net.analysis['trainDeltaT'][div,:,0]))
    n_assigned = np.sum(mask)
    net.analysis['assigned'] = np.array([n_spikes,n_assigned])
    
    print 'number of spikes: ', n_spikes
    print 'number of assigned spikes: ', n_assigned
    print 'fraction of assigned spikes: ', (n_assigned/float(n_spikes))
    
    
    histo = np.histogram(abs(net.analysis['trainDeltaT'][div,:,2][mask]),bins=np.logspace(-5,-1,101),range=[10**(-5),10**(-1)])
    
    net.analysis['deltaT_hist'] = histo[0]
    net.analysis['deltaT_hist_range'] = histo[1]
    
    
    #print net.analysis['trainDeltaT'][div,:,2][mask]
    mean_dt = np.nanmean(abs(net.analysis['trainDeltaT'][div,:,2]))
    net.analysis['mean_dt'] = mean_dt
    print 'mean spiketime deviation: ', mean_dt
    
    
    if plot == 1:
      
      #distance1norm = np.sum(abs(distanceOrthVector),axis=1)
    
      #distanceOrth = np.sqrt(np.sum(distanceOrthVector**2,axis=1))
      #self.analysis['distanceOrth'] = distanceOrth
      
      ##plot phi corr
      #plt.figure(figsize=(4,3))
      #ax0 = plt.axes([0.15,0.2,0.35,0.7])
      #ax1 = plt.axes([0.6,0.2,0.35,0.7])
      
      #ax0.hist(net.analysis['phi_corr'].ravel(),bins=100,range=[-1,1])
      #phi_corr_mean = np.mean(net.analysis['phi_corr'],axis=1)
      #print 'mean: ', np.mean(phi_corr_mean)
      #print 'len: ', len(phi_corr_mean)
      #ax1.plot(np.linspace(0,1,len(phi_corr_mean)),phi_corr_mean,'o')
      #ax1.set_ylim([-0.1,0.1])
      ##if save:
	##save_str = './pics/%s/phase_corr.pdf' % (net.scriptName.split('_')[0])
	##print 'figure saved as: %s ' % save_str
	##plt.savefig(save_str)
      #plt.show(block=False)
      
      
      print net.analysis.keys()
      plt.figure(figsize=(2,2))
      #ax0 = plt.axes([0.34,0.22,0.63,0.7])
      ax0 = plt.axes([0.34,0.45,0.63,0.5])
      ax1 = plt.axes([0.34,0.25,0.63,0.18])
      i=30
      idx_drop = np.where(net.analysis['distanceOrth'][i] < 10**(-3))[0][0]
      m = net.analysis['tau_conv'][0]
      b = net.analysis['distanceOrth'][i][idx_drop]*np.exp(-net.analysis['measureTimes'][idx_drop]/m)
      
      decay = b*np.exp(net.analysis['measureTimes']/m)
      
      #print 
      
      mean_D = np.mean(net.analysis['distanceOrth'][i][10:idx_drop-10])
      mean_D_array = np.ones(len(net.analysis['measureTimes']))*mean_D
      
      cross = np.argmin(abs(decay-mean_D_array))
      #ax1 = plt.axes([0.62,0.65,0.3,0.3])
      if np.sum(div):
	ax0.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][div][:10].transpose(),color='orangered')
      if not all(div):
	ax0.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][np.invert(div)][i].transpose(),color='dodgerblue')
      
      ax0.plot(net.analysis['measureTimes'],decay,'--r')
      
      ax0.plot(net.analysis['measureTimes'],mean_D_array,'--r')
      
      ax0.annotate(r'$\displaystyle t_{RC}$',xy=[net.analysis['measureTimes'][cross],mean_D],xytext=[net.analysis['measureTimes'][cross]*1.5,mean_D*5],arrowprops=dict(arrowstyle="->"),fontsize=12)
      
      #try:
	#ax0.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_ER'][0],':k',linewidth=2)
	#ax0.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_poisson'][0],'--k',linewidth=2)
      #except:
	#1
      #ax0.set_ylabel(r'$\displaystyle d_{\bar{\phi}}$',fontsize=14)
      ax0.set_xlabel('time in s',fontsize=12)
      
      ax0.set_ylabel(r'$\displaystyle D_{\phi}$, $\displaystyle |\Delta t|$',fontsize=14)
      
      ax0.set_yscale('log')
      #ax0.legend(loc=3)
      ax0.set_xlim([0,0.4])
      ax0.set_ylim([10**(-8),10**2.5])
      ax0.set_xticks([])
      
      ax0.set_yticks(np.logspace(-7,1,3))
      #ax0.set_yticklabels(10**np.linspace(-1,1,3),fontsize=10)
      #ax0.xaxis.get_major_formatter().set_powerlimits((0,1))
      ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%g}$' % np.log10(y)))
      
      #print net.analysis['trainDeltaT'][i,:,2]
      
      t_max = np.sum(np.isfinite(net.analysis['trainDeltaT'][i,:,1]))
      test_range = 20
      
      frac_assigned = np.zeros(t_max)
      frac_assigned[:] = np.nan
      
      for t in range(t_max):
	min_idx = max(0,t-test_range)
	max_idx = min(t+test_range,t_max)
	frac_assigned[t] = np.sum(np.isfinite(net.analysis['trainDeltaT'][i,:,2][min_idx:max_idx]))/float(max_idx-min_idx)
	
      #print 'frac assigned: ', frac_assigned
      
      
      #print net.analysis['trainDeltaT'][i,:,2]
      
      ax0.plot(net.analysis['trainDeltaT'][i,:,1],abs(net.analysis['trainDeltaT'][i,:,2]),'.k',markersize=2,label=r'$\displaystyle |\Delta t|$')
      
      #ax0.legend(numpoints=1,prop={'size':9})
      ax1.plot(net.analysis['trainDeltaT'][i,:,1][:t_max],frac_assigned,'k')
      
      ax1.set_xlim([0,0.4])
      ax1.set_xticks(np.linspace(0,0.4,3))
      ax1.set_ylim([0,1.1])
      #ax1.set_yscale('log')
      
      #ax1.set_ylim([10**(-10),10**(-0.9)])
      ax1.set_yticks(np.linspace(0,1,2))
      ax1.set_xlabel('t in s')
      ax1.set_ylabel(r'$\displaystyle f_{assigned}$')
      
      
      
      
      #if any(div):
	#ax1.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][div][:100,:time_len].transpose(),'-',color='orangered',linewidth=0.5)
      #if not all(div):
	#ax1.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][np.invert(div)][:10,:time_len].transpose(),'-',color='dodgerblue',linewidth=0.5)
      
      #try:
	#ax1.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_ER'][0],':k',linewidth=2)
	#ax1.plot(net.analysis['measureTimes'][:time_len],np.ones(time_len)*net.analysis['D_decorr_poisson'][0],'--k',linewidth=2)
      #except:
	#1
      #ax1.set_yscale('log')
      #ax1.set_xlim([0,0.03])
      #ax1.set_xticks(np.linspace(0,0.02,2))
      #ax1.set_ylim([10**(0.4),10**0.95])
      
      #for i in range(sim.steps['total']):
	#mask = np.where(net.analysis['distanceOrth'][i] > 0)
	#if net.analysis['distanceOrth'][i][mask][-1] > 0.1:
	  #col = 'orangered'
	#else:
	  #col = 'dodgerblue'
	#print net.analysis['distanceOrth'][i]-net.analysis['distanceOrth'][i]
	#ax0.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])],net.analysis['distanceOrth'][i],'-',color=col)	# orthogonalized distance
	#ax0.plot(net.analysis['measureTimes'][:time_len],net.analysis['distanceOrth'][i],'-',color=col)	# orthogonalized distance
	
      if save:
	save_str = './pics/%s/samples/decorr1_sample_nu=%d.pdf' % (net.scriptName.split('_')[0],net.topoConst['drive_rate'])
	print 'figure saved as: %s ' % save_str
	plt.savefig(save_str)
      plt.show(block=False)
    
    
    
    if plot == 2:
      
      plt.figure(figsize=(3,2))
      ax1 = plt.axes([0.3,0.6,0.65,0.35])
      ax2 = plt.axes([0.3,0.2,0.65,0.35])
      
      ax1_border = [10**(-13),10]
      
      
      ax1.plot(net.analysis['trainDeltaT'][i,:,1],abs(net.analysis['trainDeltaT'][i,:,2]),'.k',markersize=2,label=r'$\displaystyle |\Delta t|$')
      
      ax1.set_xticks([])
      ax1.set_yscale('log')
      
      ax1.set_ylim([10**(-5),10**(-0.9)])
      ax1.set_yticks(np.logspace(-5,-1,3))
      ax1.set_ylabel(r'$\displaystyle |\Delta t|$/s')
      
      #ax1.legend(loc=3,prop={'size':10})
      
      #ax2.hist(abs(net.analysis['trainDeltaT'][div,:,2][mask]),bins=np.logspace(-5,-1,101),range=[10**(-5),10**(-1)],color='k',histtype='step')
      #ax2.set_xscale('log')
      #ax2.set_yticks([])
      #ax2.set_xlabel(r'$\displaystyle \Delta t$ in s')
      
      #if save:
	#save_str = './pics/%s/assigned_spikes_scatter.pdf' % (net.scriptName.split('_')[0])
	#print 'figure saved as: %s ' % save_str
	#plt.savefig(save_str)
      
      #plt.show(block=False)
      ax2.set_xlabel('t in s')
      
      #return analysis
      #print net.analysis['distanceOrth'].shape
      slide_len = 10
      
      idx_ct = 0
      time_len = len(net.analysis['measureTimes'])
      for co_ct in range(sim.steps['co']):
	for pd_ct in range(sim.steps['pd']):
	  
	  for pd_ct_compare in range(pd_ct+1,sim.steps['pd']):
	    print "step: ", pd_ct, pd_ct_compare
	    
	    distanceVector = (net.analysis['measureStates'][co_ct][pd_ct]-net.analysis['measureStates'][co_ct][pd_ct_compare])[:,:net.Const['N']]
	    
	    #distanceOrthVector = np.zeros((time_len,net.Const['N']))
	    
	    #for t in range(time_len):
	      #distanceOrthVector[t] = OrthoVector(distanceVector[t])
	    
	    #distanceOrth = np.sqrt(np.sum(distanceOrthVector**2,axis=1))
	    
	    #if plot == 2:
	      #plt.figure(figsize=(2,1.5))
	      #ax0 = plt.axes([0.2,0.2,0.75,0.75])
	      #ax1 = plt.axes([0.3,0.3,0.65,0.45])
	      #ax0.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][idx_ct],color='dodgerblue')
	    
	    #if net.analysis['distanceOrth'][idx_ct][-1] < 0.01:
	    frac_conv = np.zeros((time_len,5))
	    
	    # obtain correlation of different frac_conv functions with delay determined by convergence time
	    # obtain statistics of maximum value of frac_conv that does not lead to convergence
	    # get distribution of single phase distances as video through time
	    # maybe do most of the stuff in pert_read already?!
	    
	    for i in range(1,3):
	      
	      peak = []
	      peak_idx = []
	      
	      col = (3-i)/3.
	      frac_conv[:,i-1] = np.sum(abs(distanceVector)<10**(-i),axis=1)/float(net.Const['N'])
	      
	      #for t in range(time_len):
		#idx_min = max(0,t-slide_len)
		#idx_max = max(0,t+slide_len+1)
		#if np.sum(frac_conv[t,i-1] <= frac_conv[idx_min:idx_max,i-1]) == 1:
		  #peak.append(frac_conv[t,i-1])
		  #peak_idx.append(t)
	      
	      #print peak
	      ax2.plot(net.analysis['measureTimes'],frac_conv[:,i-1],color=(col,col,col),label=r'$\displaystyle \theta_{%d} = 10^{-%d}$'%(i,i))
	      #ax1.plot(net.analysis['measureTimes'][peak_idx],peak,'o',color=(col,col,col))
	      
	      ax2.annotate(r'$\displaystyle \theta_{%d} = 10^{-%d}$'%(i,i),xy=[5.5,np.mean(frac_conv[:,i-1])],xytext=[6.5,np.mean(frac_conv[:,i-1])+0.15],arrowprops=dict(arrowstyle="->"),fontsize=10) #bbox={'facecolor':'white','alpha':0.9,'pad':5})
	    #if plot == 2:
	    #plt.setp(ax0, xlim=[0,1.2], xticks=[], yscale='log', ylim=[10**(-13),10**1.5], yticks=10**np.linspace(-13,1,7), ylabel=r'$\displaystyle D_{\phi}$')
	    #ax0.xaxis.get_major_formatter().set_powerlimits((0,1))
	    #ax0.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
	    
	    plt.setp(ax2, xlabel='t in s',xticks=np.linspace(0,10,6), yticks=np.linspace(0,1,3))
	    ax2.set_ylabel(r'$\displaystyle f_{\phi_i < \theta_k}$',fontsize=14)
	    #plt.setp(ax1, xticks=np.linspace(0,1.2,7), xlabel='t in s', yticks=np.linspace(0,1,3), ylabel=r'$\displaystyle f_{ct}$')
	    
	    
	    
	    #plt.legend(prop={'size':10},fancybox=True,bbox_to_anchor=(0.12,0.68),loc='lower left', shadow=True)
	    if save:
	      save_str = './pics/%s/single_phase_fraction_nu=%d.pdf' % (net.scriptName.split('_')[0],net.topoConst['drive_rate'])
	      print 'figure saved as: %s ' % save_str
	      plt.savefig(save_str)
	    plt.show(block=True)
	  idx_ct += 1
    
      #print net.analysis['threshold_peaks']
      
      for i in range(10):
	plt.figure(figsize=(4,3))
	col = (11-i)/11.
	mask = np.invert(np.isnan(net.analysis['threshold_peaks'][:,i]))
	plt.hist(net.analysis['threshold_peaks'][:,i][mask],bins=20,range=[0,1],color=(col,col,col))
	plt.show(block=False)
      
      print np.nanmax(net.analysis['threshold_peaks'],axis=0)
      plt.figure()
      plt.plot(np.arange(10),np.nanmax(net.analysis['threshold_peaks'],axis=0))
      plt.ylim([0,1])
      plt.show()
    
    
    
    
    
  
  if net.Const['special'] == 5:
    
    #print net.analysis.keys()
    #print net.analysis['measureStates'].shape
    var = np.zeros((sim.steps['co'],net.Const['N'],sim.steps['pd']))
    for co_ct in range(sim.steps['co']):
      for pd_ct in range(sim.steps['pd']):
	for n in range(net.Const['N']):
	  var[co_ct,n,pd_ct] = np.var(net.analysis['measureStates'][co_ct,pd_ct,:,n])
	
      #print var[co_ct,0,:]
  
  
  if net.Const['special'] == 7:
    #print net.analysis.keys()
    #print net.analysis['distanceOrth'].shape
    idx_decorr = np.where(net.analysis['distanceOrth'][:,-1]<1)[0]
    if plot and len(idx_decorr):
      #print idx_decorr
      #print net.analysis['distanceOrth'].shape
      #print net.analysis['distanceOrth'][idx_decorr].shape
      plt.figure()
      plt.plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][idx_decorr].transpose(),'-k')
      plt.yscale('log')
      plt.show()
  
  
  if net.Const['special'] == 8:
    #print net.analysis.keys()
    #print net.analysis['spikecross']
    
    if plot:
      idx_div = net.analysis['distanceOrth'][:,-1] > 0.1
      #print idx_div
      
      #print net.analysis['spikecross'][idx_div]
      #print net.analysis['spikecross'][np.invert(idx_div)]
      
      mask = np.invert(np.isnan(net.analysis['spikecross'][:,0]))

      print "total # of diverged traj: ", np.sum(idx_div)
      print "total # of spike crossing events: ", np.sum(mask)
      
      print "fraction: ", np.sum(idx_div)/float(np.sum(mask))
      
      #mask2 = np.invert(np.isnan(net.analysis['spikecross'][idx_div,1]))
      #print net.analysis['spikecross'][idx_div,1][mask2]
      plt.figure()
      plt.hist(net.analysis['spikecross'][idx_div&mask,1],bins=np.linspace(50,150,100),color='k',alpha=0.5)
      plt.hist(net.analysis['spikecross'][np.invert(idx_div)&mask,1],color='b',alpha=0.5)
      plt.show(block=False)
      
      plt.figure()
      plt.hist(net.analysis['spikecross'][idx_div&mask,0],bins=np.linspace(50,150,100),color='k',alpha=0.5)
      plt.hist(net.analysis['spikecross'][np.invert(idx_div)&mask,0],bins=np.linspace(50,150,100),color='b',alpha=0.5)
      plt.show(block=False)
      #for i in len(net.analysis['distanceOrth']):
	#np.
      #plt.plot()
    
    
    
  if plot==2 and net.Const['special'] == 0:
    # get title and parastring from net.Const/net.topoConst
    
  
    

    #sim_idx = np.random.randint(sim.steps['co']*sim.steps['ps'])

    #fig,axes = plt.subplots(nrows=2,ncols=1)
    #if net.topo == 'p':
      #spiked = np.zeros(net.analysis['trainNeuronDrive'].shape[-1])
      #for t in range(net.analysis['trainNeuronDrive'].shape[-1]):
	#spiked[t] = len(np.unique(net.analysis['trainNeuronDrive'][sim_idx][:t]))/float(net.Const['N'])
      
      #axes[0].plot(net.analysis['trainTimeDrive'][sim_idx],spiked,'-k')
      #axes[0].set_xlim([0,net.analysis['measureTimes'].max()])
      #axes[0].set_ylim([0,1])
      
    #axes[1].plot(net.analysis['measureTimes'],net.analysis['distanceOrth'][sim_idx])
    #axes[1].set_yscale('log')4
    #axes[1].set_ylim([10**-16,10**1])
    
    #for idx_ct in range(sim.steps['total']*sim.steps['ps']):
      #end_state = distance[idx_ct] > 0.01
      #distance[idx_ct] > 0.01      
    
    #print net.analysis['decorrStats'].transpose()[1]
    #for i in range(len(net.analysis['decorrStats'])):
      #if net.analysis['distanceOrth'][i][-1] > 0.01:
	#plt.figure()
	#plt.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])],net.analysis['distanceOrth'][i],'-b')	# orthogonalized distance
	#if not np.isnan(net.analysis['decorrStats'][i][0]):
	  #plt.plot(net.analysis['decorrStats'][i][0],net.analysis['distanceOrth'][i][net.analysis['decorrStats'][i][2]],'or')
	  #plt.plot(net.analysis['measureTimes'][net.analysis['decorrStats'][i][3]],net.analysis['distanceOrth'][i][net.analysis['decorrStats'][i][3]],'or')
	  
	  #plt.plot(net.analysis['measureTimes'],np.exp(net.analysis['decorrStats'][i][1]*(net.analysis['measureTimes']-net.analysis['decorrStats'][i][0]))*net.analysis['distanceOrth'][i][net.analysis['decorrStats'][i][0]],'--r')
	#plt.yscale('log')
	#plt.xlim([0,0.1])
	#plt.ylim([10**(-4),10])
	#plt.show()
  
    plt.figure()
    plt.show(block=False)
    if plot == 2:	# interactive mode on if video/stream
      plt.ion()
    
    
    if net.topo == 'p':
      gs = gridspec.GridSpec(2, 2)
      ax0 = plt.subplot(gs[0,0])
      ax1 = plt.subplot(gs[1,0])
      ax2 = plt.subplot(gs[1,1])
      ax3 = plt.subplot(gs[0,1])
      
      mask = np.invert(np.isnan(net.analysis['convergeVelocity'][0]))
      mean_convVel1 = np.mean(net.analysis['convergeVelocity'][0][mask])
      
      mask = np.invert(np.isnan(net.analysis['convergeVelocity'][1]))
      mean_convVel2 = np.mean(net.analysis['convergeVelocity'][1][mask])
      
      decay1 = np.exp(net.analysis['measureTimes']*mean_convVel1)
      decay2 = np.exp(net.analysis['measureTimes']*mean_convVel2)

      ax0.plot(net.analysis['measureTimes'],decay1,'--k',linewidth=5,label='approximate initial convergence')
      ax0.plot(net.analysis['measureTimes'],decay2,'--r',linewidth=5,label='approximate later convergence')
      
      ax2.plot(range(len(net.analysis['convergeVelocity'][1])),net.analysis['convergeVelocity'][1],'ro')
      ax2.set_ylabel('later convergence speed')
      ax2.set_xlabel('trial \#')

    else:
      gs = gridspec.GridSpec(1, 1)
      ax0 = plt.subplot(gs[0,0])
      #ax2 = plt.subplot(gs[1,0])
      
      conv_traj = net.analysis['decorrStats'].transpose()[1][np.where(net.analysis['decorrStats'].transpose()[1] < 0)[0]]
      decay = np.mean(conv_traj)
      #decay = -67
      decay_arr = 0.002*np.exp(net.analysis['measureTimes']*decay)
      #print decay
      ax0.plot(net.analysis['measureTimes'],decay_arr,'--k',linewidth=5,label='mean convergence rate')
      
      #ax2.plot(range(len(conv_traj)),conv_traj,'ro')
      #ax2.set_ylabel('convergence speed',fontsize=14)
      #ax2.set_xlabel('trial \#',fontsize=14)
      
      
    # try estimated decay    

    #ax0.plot(net.analysis['measureTimes'],np.abs(net.analysis['distance1']).transpose(),'-r')	# orthogonalized distance
    
    ax0.plot(net.analysis['measureTimes'][:len(net.analysis['distanceOrth'][0])],net.analysis['distanceOrth'].transpose(),'-b')	# orthogonalized distance
    #for i in range(len(net.analysis['decorrStats'])):
      #if not np.isnan(net.analysis['decorrStats'][i][0]):
	#ax0.plot(net.analysis['decorrStats'][i][0]],net.analysis['distanceOrth'][i][net.analysis['decorrStats'][i][2]],'or')
	
    ax0.set_xlim([0,0.2])
    #ax0.set_xlim([0,sim.TC])
    
    #ax0.set_ylim([10**(-16),10**1])
    ax0.set_ylim([10**(-8),10**1])
    
    ax0.set_yscale('log')
    ax0.set_xlabel('time in s',fontsize=14)
    ax0.set_ylabel('distance from reference',fontsize=14)
    ax0.legend()
    
    
    if net.topo == 'p':
      ax1.plot(range(len(net.analysis['convergeVelocity'][0])),net.analysis['convergeVelocity'][0],'ko')
      ax1.set_ylabel('initial convergence speed')
      ax1.set_xlabel('trial \#')
      
      ax3.plot(range(len(net.analysis['kinkTime'])),net.analysis['kinkTime'],'ok')

    
    plt.suptitle(title_str,fontsize=18)
    plt.show(block=False)
    
    
    
    plt.figure()
    #plt.plot(drop_time_hist[1][1:],drop_time_hist[0])
    plt.plot(net.analysis['measureTimes'],correct_traj_sum)
    plt.plot(net.analysis['measureTimes'],np.exp(m*net.analysis['measureTimes'])*np.exp(b),'--r')
    plt.xlim([0,0.5])
    plt.ylim([0,1])
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('time in s',fontsize=16)
    plt.ylabel('ratio diverged trajectories',fontsize=16)
    plt.title(title_str)
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
	    #plt.savefig('./pics/vid/%s_%d.pdf'%(png_name,t))
	  #t += 1
	  #if t >= len(net.analysis['measureTimes']) or (net.analysis['distanceOrth'][idx_ct][t-1] < 10**(-8)):
	    #break
      #if save == 2:
	## convert to video
	#video_string = "ffmpeg -f image2 -r 25 -i './pics/vid/" + png_name + "_%d.pdf' ./pics/vid/" + png_name + ".mp4"
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
  #print net.analysis.keys()
  #print net.analysis['D_decorr']
  #print net.analysis['D_decorr_bs']
  
  
  for key in ['measureTimes','train','trainPert','trainDeltaTime','trainDrive','distanceOrth','decorrStats','convergeVelocity','kinkTime','T_conv','distanceVector','distanceOrthVector','trainDeltaT','measureStates','D_decorr']:
    try:
      trash = net.analysis.pop(key)
    except:
      1
  #print net.analysis
  return net.analysis


def update_parent_slice(T,parent_slice,t,n,break_iter=0):
  
  if break_iter > len(parent_slice)**2:
    #print "waaah"
    return 0
  p_new = parent_slice[n]
  
  indices_update = np.append(np.where(T[:,n]<t)[0],np.where(T[n]<t)[0])
  #print indices_update
  p_slice_update = parent_slice[indices_update]
  #np.append(parent_slice[np.where(T[n]<t)[0]],parent_slice[np.where(T[:,n]<t)[0]])
  
  mask = np.invert(np.isnan(p_slice_update))
  if len(np.unique(p_slice_update[mask])) > 1:
    #print "new: ", p_new
    #print p_slice_update
    for m in indices_update:
      parent_slice[m] = p_new
      parent_slice[np.where(T[m]<t)[0]] = p_new
      parent_slice[np.where(T[:,m]<t)[0]] = p_new
      break_iter += update_parent_slice(T,parent_slice,t,m,break_iter)
      
  return 1


def save_analysis(analysis,results_file):
  
  if os.path.exists(results_file):
    os.remove(results_file)
  
  print 'Save results in %s' % results_file
  
  # save results in nc-file
  ncid= netcdf.netcdf_file(results_file,'w')
  #print analysis
  for key in analysis.keys():
    try:
      dim_tup = []
      for dim_len in analysis[key].shape:
	dim_name = 'dim_%s_%d' % (key,dim_len)
	if not dim_name in ncid.dimensions:
	  ncid.createDimension(dim_name,dim_len)
	dim_tup.append(dim_name)
      if len(dim_tup)==0:
	dim_name = 'one'
	if not 'one' in ncid.dimensions:
	  ncid.createDimension('one', 1)
	dim_tup = ('one',)
      elif len(dim_tup)==1:
	dim_tup = (dim_name,)
      else:
	dim_tup = tuple(dim_tup)
      dt = analysis[key].dtype

      NcVar = ncid.createVariable(key,dt.char,dim_tup)
      NcVar[:] = analysis[key]
      #print key
      #print analysis[key]
    except:
      print 'key: %s does not want to be stored. Values:' % key
      print analysis[key]
  ncid.close()

  
  
def erfi(z):
  #print z<0
  return (complex(0.,-1.)*spec.erf(complex(0.,1.)*z)).real

def analyze_drive_cur(net,sim,plot,save):
  
  K = net.Const['K']
  rateWnt = net.Const['rateWnt']
  tauM = net.Const['tauM']/1000.
  V_T = 1
  V_R = 0
  
  analyze = {}
  
  I_sim_tmp = np.zeros(sim.steps['co'])
  i=0
  for dat in os.listdir(net.Path['results'] + '/init'):
    ncid = netcdf.netcdf_file(net.Path['results'] + 'init/' + dat,'r',mmap=False)
    I_sim_tmp[i] = ncid.variables['finalCurrents'][:][0]
    ncid.close()
    i+=1
  
  mask = np.where(I_sim_tmp)
  print I_sim_tmp
  analyze['I_sim'] = np.mean(I_sim_tmp[mask])
  
  # driving current from balance equation
  analyze['I_bal'] = -rateWnt*net.Const['J0']*tauM# - net.topoConst['drive_rate']*net.topoConst['drive_cplg']*tauM
  
  # driving current from selfconsistency equation (Brunel,Hakim)
  
  
  selfcon_sv = 'data/selfcon_K%d_rateWnt%g.nc' % (K,rateWnt)
  runStr = './Code/other/calc_I %s %d %d %g' % (selfcon_sv, net.Const['N'],K,rateWnt)
  runIt(sim.cluster,runStr,net.Path['txt'])
  try:
    ncid = netcdf.netcdf_file(selfcon_sv,'r',mmap=False)
    analyze['I_selfcon'] = np.sqrt(net.Const['K'])*ncid.variables['I_0'].getValue()
    ncid.close()
  
    os.remove(selfcon_sv)
  except:
    1
  
  #print analyze
  #print analyze['I_selfcon']
  plot_V = 1
  plot_phi = 1
  
  ana_keys = ['I_sim','I_selfcon']
  label_text = ['sim.','ana.']
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  mpl.rcParams['xtick.labelsize'] = 14
  mpl.rcParams['ytick.labelsize'] = 12
  mpl.rcParams['font.size'] = 12
  
  plt.figure(figsize=(2.5,2))
  ax = plt.axes([0.15,0.15,0.75,0.7])
  #print analyze.values()
  col = [(x,x,x) for x in np.linspace(0,0.6,2)]
  for i in range(len(ana_keys)):
    if plot:
      
  
      I_plot = analyze[ana_keys[i]]
      #I_plot = analyze['I_sim']/np.sqrt(K)
      #print I_plot
      
      sigma2 = net.Const['J0']**2*rateWnt*tauM
      sigma = np.sqrt(sigma2)

      mu = V_T + I_plot + np.sqrt(K)*net.Const['J0']*rateWnt*tauM
      
      steps = 1001
      erfi_tmp = np.zeros(steps)
      
      
      if plot_V:
	V = np.linspace(-0.5,1,steps)
	
	for j in range(steps):
	  if V[j] < 0:
	    erfi_tmp[j] = erfi((V_R-mu)/sigma)
	  else:
	    erfi_tmp[j] = erfi((V[j]-mu)/sigma)
	
	P_V = np.sqrt(math.pi)*rateWnt*tauM/sigma*np.exp(-(V-mu)**2/sigma2)*(erfi((V_T-mu)/sigma)-erfi_tmp)
	#print P_V
	if not i:
	  ax.plot(V,P_V/np.sum(P_V),'--',color=col[i],label='P(V)')# from %s'%label_text[i])
    
      if plot_phi:
	phi = np.linspace(-0.5,1,steps)
	
	# determine factors of transformation
	
	T_free = tauM*np.log(1+1./I_plot)
	I_c = V_T + I_plot
	V = I_c*(1-np.exp(-phi*T_free/tauM))
	trafo_factor = I_c*T_free/tauM*np.exp(-phi*T_free/tauM)      
	
	#print T_free
	#print I_c
	
	for j in range(steps):
	  if phi[j] < 0:
	    erfi_tmp[j] = erfi((V_R-mu)/sigma)
	  else:
	    erfi_tmp[j] = erfi((V[j]-mu)/sigma)
	
	P_phi = trafo_factor*np.sqrt(math.pi)*rateWnt*tauM/sigma*np.exp(-(V-mu)**2/sigma2)*(erfi((V_T-mu)/sigma)-erfi_tmp)
	
	ax.plot(phi,P_phi/np.sum(P_phi),'-',linewidth=2,color=col[i],label=r'P($\displaystyle \phi^{%s}$)' % label_text[i])

  plt.setp(ax,xlim=[-0.5,1.1],xticks=[0,1],xticklabels=[r'$\displaystyle V_R$',r'$\displaystyle V_T$'])
  plt.setp(ax,ylim=[0,max(P_V/np.sum(P_phi))*1.2],ylabel=r'P($\displaystyle \cdot$)')
  
  ax.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',
    labelleft='off')         # ticks along the top edge are off)
  
  plt.legend(loc=[-0.15,0.7],prop={'size':11},shadow=True)
  if save:
    save_str = './pics/basics/P(V)_N=%d_K=%d_nu=%g.pdf' % (net.Const['N'],net.Const['K'],net.Const['rateWnt'])
    plt.savefig(save_str)
    print "Figure saved to %s." % save_str
  plt.show(block=False)
  
  #except:
    #print "Parameters K=%d and rateWnt=%g gave no solution!" % (K,rateWnt)
    #analyze['I_selfcon'] = np.nan
  
  return analyze


def assignSpikes(Data,idx_ref,idx_pert,idx_cut,Topo=None,find_decorr=True):
  
  assigned = []
  #print "assigning..."
  deltaTime = np.ones(idx_cut)
  deltaTime[:] = np.nan
  
  idx_cut = min(idx_cut,np.sum(np.isfinite(Data['train'][idx_ref][:,0])),np.sum(np.isfinite(Data['train'][idx_pert][:,0])))
  for i in range(idx_cut):
    
    spikeRef = Data['train'][idx_ref][i,0]
    spike = Data['train'][idx_pert][i,0]
    
    spikeTimeRef = Data['train'][idx_ref][i,1]
    spikeTime = Data['train'][idx_pert][i,1]
    
    if np.isnan(spikeRef) or np.isnan(spike):
      deltaTime[i:] = np.nan
      break
    
    if spikeRef==spike and not (i in assigned):
      #print "index: %d/%d" % (i,i)
      deltaTime[i] = abs(spikeTime-spikeTimeRef)
      assigned.append(i)
    else:
      if (Topo['adjMat'][spikeRef][spike] or Topo['adjMat'][spike][spikeRef]) and find_decorr:
	return Data['train'][np.array([idx_ref,idx_pert]),i,:]
      if not find_decorr:
	# search for the belonging spike
	idx_test_list = np.where(Data['train'][idx_pert][:,0] == spikeRef)[0]	# get pos. of neurons AP in perturbed trajectory
	test_Data_spiketimes = Data['train'][idx_pert][idx_test_list,1]
	test_Data_spiketimes_idx = abs(spikeTimeRef-test_Data_spiketimes) < 1*10**(-1)	# < 100ms difference
	idx_test_list = list(idx_test_list[test_Data_spiketimes_idx])
	
	found = 0
	while (found == 0):
	  deltaTime_test = abs(spikeTimeRef-Data['train'][idx_pert][idx_test_list,1])
	  if not len(idx_test_list):
	    #deltaTime[i] = np.nan
	    break
	  else:
	    idx_test = idx_test_list[deltaTime_test.argmin()]
	    
	    if idx_test in assigned:
	      idx_test_list.remove(idx_test)
	    else:
	      deltaTime[i] = Data['train'][idx_pert][idx_test,1]-spikeTimeRef
	      assigned.append(idx_test)
	      #print "index: %d/%d" % (idx_test,i)
	      found = 1
  return deltaTime



def T_free(I_ext,tau_M=0.01):
  
  return tau_M*np.log(I_ext/(I_ext-1.))

def get_Z(phi,I_ext,J,tau_M=0.01):
  return -tau_M/T_free(I_ext)*np.log(np.exp(-phi*T_free(I_ext)/tau_M) - J/I_ext)

def contraction(I_ext=0.1,J=None,tau_M=0.01,plt_range=None):
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  #deltaPhi_array = [1,0.5,0.2]
  deltaPhi_array = [0.5]
  
  I_ext += 1
  
  phi1 = 1
  
  try:
    if J[0]:
      assert len(J)==len(I_ext), "what?"
  except:
    J = np.linspace(-100,-0.1,10001)
  
  #if not plt_range:
    #if all(J[0]==J):
      #plt_range = I_ext
    #else:
      #plt_range = -J
  #mask = np.invert(np.isnan(I_ext))
  #J = J[mask]
  #I_ext = 1+I_ext[mask]
  
  print J
  print I_ext
  
  plt.figure()
  gs = gridspec.GridSpec(3, 1)
  
  ax0 = plt.subplot(gs[0,0])
  ax1 = plt.subplot(gs[1,0])
  ax2 = plt.subplot(gs[2,0])
  
  
  label_text = ['poisson','recurrent']
  #print get_Z(phi1,I_ext,J)
  
  for i in range(2):
    if i:
      J = -np.ones(len(I_ext))*0.1
    deltaPhi = deltaPhi_array[0]
    phi2 = phi1-deltaPhi
    Z = (get_Z(phi1,I_ext,J) - get_Z(phi2,I_ext,J))/deltaPhi
    #print Z
    timeToSpike1 = (1-get_Z(phi1,I_ext,J))*T_free(I_ext)
    timeToSpike2 = (1-get_Z(phi2,I_ext,J))*T_free(I_ext)
    
    deltaTTS = -(timeToSpike1 - timeToSpike2)
    #print -J
    #print get_Z(phi2,I_ext,J)
    ax0.plot(plt_range,get_Z(phi2,I_ext,J),label=label_text[i])
    ax1.plot(plt_range,deltaTTS)
    ax2.plot(plt_range,Z)
    
  ax0.set_ylim([-1,1])
  ax0.set_ylabel(r'$\displaystyle \phi^+$',fontsize=16)
  #ax0.set_yscale('log')
  ax0.legend(prop={'size':16})
  ax0.set_xscale('log')
  
  ax1.set_ylim([10**(-6),1])
  ax1.set_yscale('log')
  ax1.set_ylabel(r'$\displaystyle \Delta $Time to spike / $\displaystyle T_{free}$',fontsize=16)
  ax1.set_xscale('log')
  
  ax2.set_ylim([10**(-2),1])
  ax2.set_yscale('log')
  #ax2.set_yticks([0.1,0.3,1])
  #ax2.set_yticklabels([0.1,0.3,1])
  ax2.set_xlabel(r'$\displaystyle J_p$',fontsize=18)
  #ax2.set_xlabel(r'$\displaystyle \bar{\nu}_p',fontsize=16)
  
  ax2.set_ylabel(r'$\displaystyle \Delta\phi^+/\Delta\phi^-$',fontsize=16)
  ax2.set_xscale('log')
  try:
    plt.suptitle(r'External current $\displaystyle I_{ext}=%g$'%I_ext,fontsize=18)
    plt.savefig('./pics/poisson/phase_contraction_I_ext=%g.pdf'%I_ext)
  except:
    1
  
  plt.show(block=False)
  
  

def readDataStats(fileName,Data,read_list,numsim,idx_file):
  
  ncid = netcdf.netcdf_file(fileName,'r',mmap=False)
  #print ncid.variables.keys()
  for key in read_list:
    if not (key in ncid.variables.keys()):
      continue
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
    ncid = netcdf.netcdf_file(net.Path['results'] + 'results_processed.nc','r',mmap=False)
    
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
    ncid = netcdf.netcdf_file(fileName,'r',mmap=False)  
    
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
	  ncidTop = netcdf.netcdf_file(link_topology,'r',mmap=False)
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
	  DataRef['deltaPhase'] = abs(DataRef['measureStates'][:idx_phase] - Data['measureStates'][:idx_phase])
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

  #save_str = './pics/pert_spread/phase_%3.1f_%3.1f_%3.1f_%3.1f_N=%d.pdf' % (net.const['alpha_recip'],net.const['alpha_conv'],net.const['alpha_div'],net.const['alpha_chain'],net.const['N'])
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




  
  
################## file access and reading ##################

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


def pick_file(scriptName, infoPath, file_excl='', restrictions=None,pick=None):
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
	if not pick:
	  print string
	fileName.append(dat)
	idx += 1
  
  if pick:
    ch_idx = pick
  else:
    ch_idx = int(raw_input('Pick a file: '))
  return fileName[ch_idx]


def clean_erase():
  scriptPath = 'data/SONET/'
  
  for dat in os.listdir(scriptPath):
    #print dat
    if 'results' in dat:
      for fileName in os.listdir(scriptPath + dat):
	if not ('results_processed' in fileName):
	  os.remove(scriptPath + dat + '/' + fileName)
      #break


def erase(scriptPath=None,infoName=None,Hash=None,restrictions=None):
  if infoName and scriptPath:
    tmp_para, tmp_topo, tmp_stp, Hash = read_from_info(infoName)
    print 'removing simulation: %s' % scriptPath
    for dat in os.listdir(scriptPath):
      if Hash in dat:
	shutil.rmtree(scriptPath + '/'+ dat,ignore_errors=True)
    os.remove(infoName)
  else:
    scratch = int(raw_input('You want to access data from local (0) or from scratch (1,2,3,4)? '))
    if scratch:
      sessionPath = '/scratch%02d/aschmidt/data/' % scratch	#path for storage on scratch
    else:
      sessionPath = 'data/'			#path for local storage
    
    scriptName = pick_script(sessionPath)
    
    info_dir = '%sinfo/%s/' % (sessionPath,scriptName)
    #print scriptName
    not_satisfied = 'y'
    while not_satisfied in ['y','all']:
      
      if not_satisfied == 'all':
	fileName = pick_file(scriptName,info_dir,restrictions=restrictions,pick=1)
      else:
	fileName = pick_file(scriptName,info_dir,restrictions=restrictions)
      infoName = info_dir + fileName
      net = network('r',infoName)
      sim = simulation_data('r',infoName,net)
      #net.setup(sim,check=0,call=1)
      net.scriptName = scriptName
      
      print 'removing data with Hash: %s' % net.Hash
      
      for dat in os.listdir(sessionPath + scriptName):
	if net.Hash in dat:
	  shutil.rmtree(sessionPath + scriptName + '/' + dat,ignore_errors=True)
      os.remove(infoName)
      
      
      if not_satisfied == 'all':
	remove_all = 1
      else:
	not_satisfied = raw_input('Another one? ')
      
