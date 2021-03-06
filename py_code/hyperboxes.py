import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time, os, sys
from scipy.io import netcdf
import matplotlib.colors as mcolors
import matplotlib as mpl

#def plot_setup_3D(koords):
  #fig = plt.figure()
  #ax = fig.add_subplot(111)#, projection='3d')
  #plt.ion()
  #ax.scatter(koords[0],koords[1],koords[2],'ok')
  #ax.set_xlabel('x')
  #ax.set_ylabel('y')
  #ax.set_zlabel('z')
  #ax.set_xlim(-0.1,1)
  #ax.set_ylim(-0.1,1)
  #ax.set_zlim(-0.1,1)
  #plt.show()
  #return fig, ax

class network:
  
  def __init__(self,N,K,J,rateWnt,poisson,t_max):

    self.phi_T, self.V_T = 1., 1.
    self.phi_R, self.V_R = 0., 0.
    
    self.tauM = 0.01

    self.N = N
    self.K = K
      
    self.rateWnt = rateWnt
    self.poisson = poisson
    #print poisson
    
    if poisson['rate']:
      self.poisson = {'rate':poisson['rate']}
      if 'train' in poisson.keys():
	self.poisson['train_def'] = True
	self.poisson['train'] = poisson['train']
      else:
	self.poisson['train_def'] = False
    else:
      self.poisson['train_def'] = False
      
    self.TR = 2
    self.TW = 0.1
    self.TC = t_max
    
    # construct network
    self.A = np.random.random((N,N)) < K/float(N) - np.identity(N)
    self.J = self.A * J/np.sqrt(K)
    
    if poisson['rate']:
      self.A = np.append(self.A,np.identity(N),axis=0)
      self.J = np.append(self.J,np.identity(N)*poisson['cplg']/np.sqrt(K),axis=0)
    
    self.A = (self.A).astype('bool')
    
    
  def reset(self,rec=True,poi=True):
    self.s = np.zeros(len(self.A)).astype('int')		# idx of spike
    self.t = 0				# idx of measurement
    self.T = 0				# current time of the simulation
    
    if rec:
      self.phi = np.random.random(self.N)
	
    if self.poisson['rate'] and poi:
      if rec:
	self.phi = np.append(self.phi,np.zeros(self.N))
      #else:
      for n in range(self.N):
	self.phi[self.N+n] = self.phi_T - self.poisson['train'][n][0]/self.T_free
  
  
  def generate_poisson_train(self,t_max):
    
    self.poisson['train'] = []
    
    for n in range(self.N):
      self.poisson['train'].append([])
      T = 0
      while T < t_max:
	spikeTime = -np.log(1-np.random.random())/self.poisson['rate']
	self.poisson['train'][n].append(T + spikeTime)

	T += spikeTime
    #print self.poisson['train']

  def ratefinding(self):
    
    # get first guess of driving current from analytics
    fileName = 'get_analytic_guess.nc'
    
    ana_get_str = './Code/other/calc_I %s %d %d %g' % (fileName,self.N,self.K,self.rateWnt)
    os.system(ana_get_str)
    try:
      ncid = netcdf.netcdf_file(fileName,'r',mmap=False)
      #print ncid.variables['I_0'].getValue()
      I_ext = np.sqrt(self.K)*ncid.variables['I_0'].getValue()
      self.D_decorr = ncid.variables['D_decorr'].getValue()
      self.D_spike_failure = ncid.variables['D_spike_failure'].getValue()
      ncid.close()
      
      os.remove(fileName)
    except:
      I_ext = 0.1
    
    #print self.poisson
    #print np.invert(self.poisson['train_def'])
    if self.poisson['rate'] and np.invert(self.poisson['train_def']):
      self.generate_poisson_train(self.TR)
    #elif self.poisson['rate']:
      #self.phi[self.N:] = self.TR*2
    
    # rate finding factors from LEquipe!
    rateTmp, rateUp, rateDown = 0., 0., 0.
    factorTmp, factorUp, factorDown = 1., 1., 0.
    factorTmpTmp, rateTmpTmp = 0., 0.
    
    iter_rate = 0
    
    while abs(1 - rateTmp/self.rateWnt) > 0.01 and (iter_rate < 30):
      
      self.I_ext = 1 + factorTmp * I_ext
      self.T_free = -self.tauM*np.log((self.V_T - self.I_ext)/(self.V_R - self.I_ext))
      
      self.reset(rec=True,poi=True)
      if self.poisson['train_def']:
	self.phi[self.N:] = -10000

      self.measures = [self.TR]
      
      while self.T < self.TR:
	
	next_spike, dt, dt_bool = self.simple_iteration()
	
	if dt_bool:		#if next spike is reached before measurement
	  self.update_spike(next_spike)
	
      rateTmp = np.sum(self.s[:self.N])/float(self.TR*self.N)
      
      # rate finding from LEquipe!
      if (rateTmp > self.rateWnt):
	factorUp = factorTmp
	rateUp = rateTmp
	
	if (iter_rate == 0):
	  factorTmpTmp=factorUp
	  rateTmpTmp=rateUp

      else:
	factorDown = factorTmp
	rateDown = rateTmp
      
      if ((iter_rate==1) and not (factorTmpTmp==0) and (rateTmp > self.rateWnt)):
	factorDown = -rateUp*(factorUp - factorTmpTmp)/(rateUp - rateTmpTmp) + factorUp   # adjusts lower bound of current to 0-rate value (otherwise it is assumed to be 0)
	      
      # guess the next factor for the next current (assuming a linear slope)
      factorTmp = factorDown + (self.rateWnt - rateDown)/(rateUp - rateDown)*(factorUp - factorDown)

      # At the beginning the current rate can be below the wanted rate, then the following is necessary to find the initial upper bound for the bisection.
      if (rateUp < self.rateWnt):
	factorUp *= 2 # to put it above desired current (factorUp here can be negative)
	factorTmp = factorUp
      
      print 'With driving current I = %g average firing rate of network: %g' % (self.I_ext,rateTmp)
      iter_rate += 1
    
    
  def warmup(self):
    
    self.measures = [self.TW]
    
    if self.poisson['rate'] and not self.poisson['train_def']:
      self.generate_poisson_train(self.TW)
    #elif self.poisson['rate']:
      #self.phi[self.N:] = self.TR*2
      
    self.reset(rec=True,poi=True)
    
    if self.poisson['train_def']:
      self.phi[self.N:] = -100
    
    while self.T < self.TW:
      next_spike, dt, dt_bool = self.simple_iteration()
      
      if dt_bool:		#if next spike is reached before measurement
	self.update_spike(next_spike)
  
    
  def simple_iteration(self):
    
    next_spike = np.argmax(self.phi)
    
    dt_spike = (self.phi_T - self.phi[next_spike])*self.T_free
    dt_measure = self.measures[self.t]-self.T
    
    dt = min(dt_spike,dt_measure)
    self.phi += dt / self.T_free
    
    self.T += dt			# update current time
    
    dt_bool = dt_spike < dt_measure
    
    if not dt_bool:
      self.t += 1
    
    return next_spike, dt, dt_bool
    
  
  def update_spike(self,next_spike):
    
    self.s[next_spike] += 1
    
    #print "next spike: ", next_spike
    if next_spike < self.N:
      self.phi[next_spike] = 0
    else:
      next_poisson = next_spike-self.N
      #print self.poisson['train'][next_poisson]
      #print "next poisson spike (%d): %g" % (self.s[next_spike],next_poisson)
      #print "time: ", self.poisson['train'][next_poisson][self.s[next_spike]]
      self.phi[next_spike] = self.phi_T - (self.poisson['train'][next_poisson][self.s[next_spike]] - self.T)/self.T_free
      
      #print next_poisson
      
    self.phi[self.A[next_spike]] = self.PRC(next_spike)	# implement J

  
  def PRC(self,next_spike):
    return -self.tauM/self.T_free * np.log(np.exp(-self.phi[self.A[next_spike]]*self.T_free/self.tauM) - self.J[next_spike][self.A[next_spike]]/self.I_ext)

    
    
    
def simulation(N,K,J,rateWnt,poisson={'rate':0,'cplg':0,'train':None},n_steps=1,t_steps=101):
  
  #n_steps = 10
  #t_steps = 101
  
  t_max = 1
  
  Data = {}
  
  Data['measureTimes'] = np.linspace(0,t_max,t_steps)
  
  trainNeuron = []
  trainTimes = []	
  train_len = 0
  
  Data['phases'] = np.zeros((n_steps,t_steps,N))
  
  net = network(N,K,J,rateWnt,poisson,t_max)	# separate generating initial conditions from generating topology
  
  net.ratefinding()
  
  if not net.poisson['train_def']:
    
    idx_drive_cut = poisson['rate']*t_max*2
    net.poisson['train_calc'] = np.zeros((n_steps,N,idx_drive_cut))
    net.poisson['train_calc'][:] = t_max*2
    
    for n in range(n_steps):
      net.generate_poisson_train(net.TC)
      
      #print net.poisson['train_calc'][n]
      #print net.poisson['train']
      for j in range(N):
	print len(net.poisson['train'][j]), idx_drive_cut
	idx_cut = min(len(net.poisson['train'][j]),idx_drive_cut)
	net.poisson['train_calc'][n,j,:idx_cut] = net.poisson['train'][j][:idx_cut]
  
  
  
  #print net.J
  for n in range(n_steps):
    sys.stdout.write('\r{0}'.format(n))
    sys.stdout.flush()
    
    net.warmup()
    
    if not net.poisson['train_def']:
      net.poisson['train'] = np.copy(net.poisson['train_calc'][n])
      
      for j in range(N):
	#try:
	  idx_drive_on = np.where(net.poisson['train'][j] > 0.5)[0][0]
	  num_spikes_pert = len(net.poisson['train'][j]) - idx_drive_on
	  
	  
	  
	  idx_drive_on_ref = np.where(net.poisson['train_calc'][0,j] > 0.5)[0][0]
	  num_spikes_ref = len(net.poisson['train_calc'][0,j]) - idx_drive_on_ref
	  #print net
	  num_spikes = min(num_spikes_ref,num_spikes_pert)
	  
	  net.poisson['train'][j,idx_drive_on:idx_drive_on+num_spikes] = net.poisson['train_calc'][0,j,idx_drive_on_ref:idx_drive_on_ref+num_spikes]
	#except:
	  #1
	
      print net.poisson['train']
      #for k in range(len(net.poisson['train'])):
	#while net.poisson['train'][k][0] < 1:
	  #trash = net.poisson['train'][k].pop(0)
    
    net.reset(rec=False,poi=True)
    
    net.measures = Data['measureTimes']
    
    trainNeuron.append([])
    trainTimes.append([])
    last_spike = 0
    
    while net.T < net.TC:
      
      #if net.T < 1:
	#net.J = 0
      #else:
	#net.J = -1
      next_spike, dt, dt_bool = net.simple_iteration()
      
      if dt_bool:		#if next spike is reached before measurement
        net.update_spike(next_spike)
        
        # save spike train
        trainNeuron[n].append(next_spike)
        trainTimes[n].append(net.T)
        
      else:
        # save phases at measure times
        Data['phases'][n][net.t-1] = net.phi[:net.N]

        if n_steps == 1:
          print 'total number of spikes: %d, resulting firing rate: %g' % (np.sum(net.s[:net.N]),np.sum(net.s[:net.N])/(t_max*net.N))
        
        if len(trainNeuron[n]) > train_len:
          train_len = len(trainNeuron[n])
  
  #Data['phases'] = Data['phases'].transpose(1,2,0)
  Data['trainNeuron'] = np.zeros((n_steps,train_len))
  Data['trainTimes'] = np.zeros((n_steps,train_len))
  
  Data['trainNeuron'][:] = np.nan
  Data['trainTimes'][:] = np.nan
  
  for n in range(n_steps):
    Data['trainNeuron'][n][:len(trainNeuron[n])] = np.array(trainNeuron[n])
    Data['trainTimes'][n][:len(trainTimes[n])] = np.array(trainTimes[n])
  
  #print Data['trainTimes']
  
  #print Data['trainTimes'].shape
  print "schteps: ", t_steps
  fileName = 'data/hyper/%s.nc' % construct_filename(N,K,J,rateWnt,poisson,n_steps,t_steps)
  
  saveData(fileName,Data)
  

def plotData(N,K,J,rateWnt,poisson={'rate':0},n_steps=1,t_steps=101,plot=1,save=0):
  
  fileName = 'data/hyper/%s.nc' % construct_filename(N,K,J,rateWnt,poisson,n_steps,t_steps)
  
  print 'reading %s' % fileName
  bins = 56
  
  Data = readData(fileName)
  
  Data['phase_hist'] = np.zeros((t_steps,bins))
  Data['phase_hist1'] = np.zeros((t_steps,bins))
  Data['phase_hist2'] = np.zeros((t_steps,bins))
  
  for t in range(t_steps):
    if n_steps == 1:
      Data['phase_hist'][t] = np.histogram(Data['phases'].transpose(0,1,2)[0][t],bins=bins,range=[-0.1,1])[0]/float(N)
      Data['phase_hist_total'] = np.histogram(Data['phases'],bins=bins,range=[-0.1,1])[0]/float(N*t_steps)
    else:
      Data['phase_hist'][t] = np.histogram(Data['phases'].transpose(2,1,0)[0][t],bins=bins)[0]/float(n_steps)
      Data['phase_hist1'][t] = np.histogram(Data['phases'].transpose(2,1,0)[1][t],bins=bins)[0]/float(n_steps)
      Data['phase_hist2'][t] = np.histogram(Data['phases'].transpose(2,1,0)[2][t],bins=bins)[0]/float(n_steps)
      
      Data['phase_hist_total'] = np.histogram(Data['phases'][0],bins=bins,range=[-0.1,1])[0]/float(N*t_steps)
  
  
  #print Data['trainNeuron']
  #print np.where(Data['trainNeuron']==2)
  
  
  print 'firing rate of neuron 1: %g' % (len(np.where(Data['trainNeuron']==0)[0])/float(max(Data['measureTimes'])*n_steps))
  print 'firing rate of neuron 2: %g' % (len(np.where(Data['trainNeuron']==1)[0])/float(max(Data['measureTimes'])*n_steps))
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='Arial')
  
  mpl.rcParams['font.size'] = 30
  mpl.rcParams['xtick.labelsize'] = 25
  mpl.rcParams['ytick.labelsize'] = 25
  
  #fig = plt.figure(figsize=(6,3))
  if plot == 1:
    print Data['trainTimes']
    print Data['trainNeuron']
    n = 1
    #for n in range(len(Data['trainTimes'])):
    mask = np.invert(np.isnan(Data['trainTimes'][n]))
    mask_ref = np.invert(np.isnan(Data['trainTimes'][0]))
    
    print Data['phases'].shape
    
    
    #num_dat = len(distance)
    #mean_distance = np.mean(distance[int(0.55*num_dat):int(0.75*num_dat)])
    #drop = np.where(distance<0.4)[0][0]
    
    plt.figure(figsize=(2,2))
    ax_dist = plt.axes([0.2,0.2,0.75,0.75])
    
    #ax_dist.plot([0,2],[mean_distance,mean_distance],'-r',linewidth=5,alpha=0.6)
    #ax_dist.plot([Data['measureTimes'][drop],Data['measureTimes'][drop]],[10**(-5),10],'-b',linewidth=5,alpha=0.6)
    
    N_max = 200
    
    #mask_plot = Data['trainNeuron'][n]<N_max
    #ax_train.plot(Data['trainTimes'][n][mask&mask_plot],Data['trainNeuron'][n][mask&mask_plot],'ok',markersize=4,markeredgecolor='none')
    
    #mask_plot = Data['trainNeuron'][0]<N_max
    #ax_train.plot(Data['trainTimes'][0][mask_ref&mask_plot],Data['trainNeuron'][0][mask_ref&mask_plot],'or',markersize=4,markeredgecolor='none')
    #max_idx = np.where(distance < 10**(-5))[0][0]
    for i in range(len(Data['phases'])):
      distanceVector = (Data['phases'][0]-Data['phases'][i])
      distance = np.sqrt(np.sum(distanceVector**2,axis=1))
      ax_dist.plot(Data['measureTimes'],distance,'k')
    ax_dist.set_yscale('log')
    #ptrain = ax_train.scatter(Data['trainTimes'][n][mask],Data['trainNeuron'][n][mask],'k')
    
    ax_dist.plot([1,1],[10**(-5),10**1],'--k')
    #ax_dist.annotate(r'$\displaystyle D_{\phi}^{dc}$', xy=[1.7,mean_distance], xytext=[1.8,mean_distance*0.03], arrowprops=dict(arrowstyle="->"),fontsize=20)
    #ax_dist.annotate(r'$\displaystyle \tau_{RC}$', xy=[Data['measureTimes'][drop],10**(-3)],xytext=[Data['measureTimes'][drop]*0.8,10**(-4)],arrowprops=dict(arrowstyle="->"),fontsize=25)
    
    ax_dist.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
    ax_dist.set_yticks(np.logspace(-4,0,3))
    #ax_dist.set_xticks([0,1,Data['measureTimes'][drop],2])
    #ax_dist.set_xticklabels([-1,0,r'$\displaystyle \tau_{RC}$',1])
    
    ax_dist.set_ylim([10**(-5),10**1])
    ax_dist.set_ylabel(r'$\displaystyle D_{\phi}$')
    #ax_dist.set_ylabel(u'U+03F5')
    ax_dist.set_xlabel('t in s')
    #ax_dist.set_xlim([0,max(Data['measureTimes'])])
    
    ax_dist.spines['top'].set_color('none')
    ax_dist.spines['right'].set_color('none')
    #ax_corr.legend(bbox_to_anchor=(0.8,0.95),loc='upper left',borderaxespad=0,prop={'size':10},ncol=2)
    ax_dist.xaxis.set_ticks_position('bottom')
    ax_dist.yaxis.set_ticks_position('left')
    
    #ax_train.set_xlim([0,max(Data['measureTimes'])])
    #ax_train.set_xticks([])
    #ax_train.set_yticks([0,N_max])
    #ax_train.plot([1,1],[0,N_max],'--k')
    #ax_train.set_ylim([0,N_max])
    #ax_train.set_ylabel(r'\# neuron')
    
    #
    
    plt.show()
    
    
    
    
  if plot == 3:
    
    t_max = 1
    N_max = 500
    n=6
    #mask_plot = Data['trainNeuron'][n]<N_max
    #ax_train.plot(Data['trainTimes'][n][mask&mask_plot],Data['trainNeuron'][n][mask&mask_plot],'ok',markersize=4,markeredgecolor='none')
    #for i in range(len(Data['phases'])):
    distanceVector = (Data['phases'][0]-Data['phases'][n])
    distance = np.sqrt(np.sum(distanceVector**2,axis=1))
    
    num_dat = len(distance)
    mean_distance = np.mean(distance[int(0.55*num_dat):int(0.75*num_dat)])
    drop = np.where(distance<0.4)[0][0]
    
    plt.figure(figsize=(5,4))
    ax_train = plt.axes([0.25,0.44,0.7,0.52])
    ax_dist = plt.axes([0.25,0.17,0.7,0.23])
    
    ax_dist.plot([0,2],[mean_distance,mean_distance],'-r',linewidth=5,alpha=0.6)
    ax_dist.plot([Data['measureTimes'][drop],Data['measureTimes'][drop]],[10**(-5),10],'-b',linewidth=5,alpha=0.6)
    mask_plot = Data['trainNeuron'][0]<N_max
    mask_ref = np.invert(np.isnan(Data['trainTimes'][0]))
    ax_train.plot(Data['trainTimes'][0][mask_plot],Data['trainNeuron'][0][mask_plot],'or',markersize=4,markeredgecolor='none')
    mask_plot = Data['trainNeuron'][n]<N_max
    ax_train.plot(Data['trainTimes'][n][mask_plot],Data['trainNeuron'][n][mask_plot],'ok',markersize=4,markeredgecolor='none')
    ax_train.plot([0.5,0.5],[0,N_max],'--k')
    #print Data['trainTimes'][1][mask_plot]
    #print Data['trainNeuron'][1][mask_ref&mask_plot]
    #print 
    
    max_idx = np.where(distance < 10**(-5))[0][0]
    ax_dist.plot(Data['measureTimes'][:max_idx],distance[:max_idx],'k')
    ax_dist.set_yscale('log')
    #ptrain = ax_train.scatter(Data['trainTimes'][n][mask],Data['trainNeuron'][n][mask],'k')
    
    ax_dist.plot([0.5,0.5],[10**(-5),10**1],'--k')
    ax_dist.annotate(r'$\displaystyle D_{\phi}^{dc}$', xy=[0.8,mean_distance], xytext=[0.85,mean_distance*0.005], arrowprops=dict(arrowstyle="->"),fontsize=20)
    #ax_dist.annotate(r'$\displaystyle \tau_{RC}$', xy=[Data['measureTimes'][drop],10**(-3)],xytext=[Data['measureTimes'][drop]*0.8,10**(-4)],arrowprops=dict(arrowstyle="->"),fontsize=25)
    
    ax_dist.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: r'$\displaystyle 10^{%d}$' % int(np.log10(y))))
    ax_dist.set_yticks(np.logspace(-4,0,2))
    ax_dist.set_xticks([0,0.5,Data['measureTimes'][drop],1])
    ax_dist.set_xticklabels([r'$\displaystyle -0.5$',r'$\displaystyle 0$',r'$\displaystyle \tau_{RC}$',r'$\displaystyle 0.5$'])
    ax_dist.set_xlim([0,t_max])
    
    ax_dist.set_ylim([10**(-5),10**1.4])
    ax_dist.set_ylabel(r'$\displaystyle D_{\phi}$')
    #ax_dist.set_ylabel(u'U+03F5')
    ax_dist.set_xlabel('t in s')
    
    ax_dist.spines['top'].set_color('none')
    ax_dist.spines['right'].set_color('none')
    #ax_corr.legend(bbox_to_anchor=(0.8,0.95),loc='upper left',borderaxespad=0,prop={'size':10},ncol=2)
    ax_dist.xaxis.set_ticks_position('bottom')
    ax_dist.yaxis.set_ticks_position('left')
    
    ax_train.set_xlim([0,t_max])
    ax_train.set_xticks([])
    ax_train.set_yticks([0,N_max])
    ax_train.plot([1,1],[0,N_max],'--k')
    ax_train.set_ylim([0,N_max])
    ax_train.set_ylabel(r'\# neuron')
    
    if save:
      save_str = './pics/hyperbox/drive_on.svg'
      plt.savefig(save_str)
      print "Figure saved to %s." % save_str
    plt.show()

  if plot == 2:
    fig,axes = plt.subplots(nrows=4,ncols=1,figsize=(6,3))
    
    axes[0].set_position([0.1,0.6,0.18,0.33])
    axes[1].set_position([0.3,0.6,0.18,0.33])
    axes[2].set_position([0.57,0.6,0.18,0.33])
    axes[3].set_position([0.77,0.6,0.18,0.33])

    plt.setp(axes[1].get_yticklabels(), visible=False)
    plt.setp(axes[3].get_yticklabels(), visible=False)
    #plt.setp(ax21.get_xticklabels(), visible=False)

    ax21 = plt.axes([0.1,0.35,0.85,0.175])
    ax22 = plt.axes([0.1,0.15,0.85,0.175])
    #ax = plt.axes([0.05,0.05,0.45,0.9])#,projection='3d')
    #ax_hist = plt.axes([0.55,0.55,0.4,0.4])
    #ax_hist2 = plt.axes([0.55,0.05,0.4,0.4])
    
    axes[0].set_xlabel(r'$\displaystyle \phi_1$')
    axes[1].set_xlabel(r'$\displaystyle \phi_1$')
    axes[2].set_xlabel(r'$\displaystyle \phi_1$')
    axes[3].set_xlabel(r'$\displaystyle \phi_1$')
    
    axes[0].set_ylabel(r'$\displaystyle \phi_2$')
    
    axes[0].xaxis.set_label_position('top')
    axes[1].xaxis.set_label_position('top')
    axes[2].xaxis.set_label_position('top')
    axes[3].xaxis.set_label_position('top')
    
    hist_range = np.linspace(-0.1,1,bins)
    
    plt.ion()
    
    #ax_hist2.plot(hist_range,np.histogram(Data['phases'].transpose(2,1,0)[0],bins=bins)[0]/float(t_steps*n_steps),'b')
    #ax_hist2.plot(hist_range,np.histogram(Data['phases'].transpose(2,1,0)[1],bins=bins)[0]/float(t_steps*n_steps),'g')
    #ax_hist2.plot(hist_range,np.histogram(Data['phases'].transpose(2,1,0)[2],bins=bins)[0]/float(t_steps*n_steps),'r')
    
    plt.show()
    
    
    #rb = mcolors.LinearSegmentedColormap.from_list(name='red_black',colors=[(0,(0,0,0)),(1,(1,0,0))],N=2)
    
    Data['phases'] = Data['phases'].transpose(1,2,0)
    
    #print Data['phases'].shape
    ax21.plot(Data['measureTimes'],Data['phases'][:,0],linewidth=0.2)
    ax22.plot(Data['measureTimes'],Data['phases'][:,1],linewidth=0.2)
    
    ax21.set_xlim([0.,0.1])
    ax22.set_xlim([0.,0.1])
    
    ax21.set_ylim([-1,1])
    ax22.set_ylim([-1,1])
    
    ax21.set_yticks(np.linspace(-1,1,3))
    ax22.set_yticks(np.linspace(-1,1,3))
    
    ax22.set_xticks(np.linspace(0,0.1,6))
    ax22.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y,pos: '%d'%(1000*y)))
    #ax22.set_xticklabels(np.linspace(0,100,6))
    
    ax22.set_xlabel('t in ms')
    ax21.set_ylabel(r'$\displaystyle \phi_1$')
    ax22.set_ylabel(r'$\displaystyle \phi_2$')
    
    plt.setp(ax21.get_xticklabels(), visible=False)
    
    
    a = 0
    for t in range(t_steps):
      plt_point = axes[a].scatter(Data['phases'][t][0],Data['phases'][t][1],marker='.',color='k')#,label='t = %g ms' % (Data['measureTimes'][t]*1000))#,phi_sv.transpose(1,2,0)[i][2])
      
      axes[a].set_xlim([-1,1])
      axes[a].set_ylim([-1,1])
      
      ## plot statistics: Histogramm, current firing rate (over last 10ms or so), ...?
      #phist = ax_hist.plot(hist_range,Data['phase_hist'][t],'b')
      #phist = ax_hist.plot(hist_range,Data['phase_hist1'][t],'g')
      ##phist = ax_hist.plot(hist_range,Data['phase_hist2'p][t],'r')
      #phist = ax_hist.plot(hist_range,Data['phase_hist_total'],'--k')
      
      #ax_hist.set_xlim([-0.1,1])
      #ax_hist.set_ylim([0,0.05])
      
      #ax_hist2.set_xlim([-0.1,1])
      #ax_hist2.set_ylim([0,0.05])
      axes[a].set_xticks(np.linspace(-1,1,3))
      axes[a].set_yticks(np.linspace(-1,1,3))
      #axes[a].legend(loc=4,prop={'size':10},numpoints=0)#bbox_to_anchor=(0.5,1.05),loc='lower center',borderaxespad=0,prop={'size':12})
      text_handle = axes[a].text(-0.9, -0.9, r'$\displaystyle t = %g\,$ms' % (Data['measureTimes'][t]*1000), bbox={'facecolor':'white','alpha':1,'pad':5},fontsize=12)
      
      fig.canvas.draw()
      #print Data['measureTimes'][t]
      if (Data['measureTimes'][t] in [0.02,0.021,0.08,0.081]):
	a += 1
	
	if a >= 4:
	  if save:
	    save_str = './pics/hyperbox/snapshot_N=%d_K=%d_nu=%d.pdf' % (N,K,rateWnt)
	    plt.savefig(save_str)
	    print 'Figure saved to ' + save_str	
	  break
      else:
	plt_point.remove()
	text_handle.remove()
      #ptrain.remove()
      #ax_hist.cla()
    plt.ioff()
  

    
def construct_filename(N,K,J,rateWnt,poisson,n_steps=1,t_steps=101):
  
  if poisson['rate'] and not 'train' in poisson.keys():
    return 'N=%d_K=%d_J=%g_rateWnt=%g_poisson_rate=%g_cplg=%g_steps_n=%d_t=%d' % (N,K,J,rateWnt,poisson['rate'],poisson['cplg'],n_steps,t_steps)
  elif poisson['rate']:
    return 'N=%d_K=%d_J=%g_rateWnt=%g_poisson_def_cplg=%g_steps_n=%d_t=%d' % (N,K,J,rateWnt,poisson['cplg'],n_steps,t_steps)
  else:
    return 'N=%d_K=%d_J=%g_rateWnt=%g_steps_n=%d_t=%d' % (N,K,J,rateWnt,n_steps,t_steps)
  
def saveData(fileName,Data):
  
  print "saving file in %s" % fileName
  
  # save results in nc-file
  ncid = netcdf.netcdf_file(fileName,'w')
  print Data.keys()
  for key in Data.keys():
    dim_tup = []
    for dim_len in Data[key].shape:		# add dimensions
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
    dt = Data[key].dtype

    NcVar = ncid.createVariable(key,dt.char,dim_tup)
    NcVar[:] = Data[key]
  ncid.close()
  
def readData(fileName):
  
  Data = {}
  
  ncid = netcdf.netcdf_file(fileName,'r',mmap=False)
  
  for key in ncid.variables.keys():
    Data[key] = ncid.variables[key][:]
  
  ncid.close()
  
  return Data
    