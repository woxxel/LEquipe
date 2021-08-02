import numpy as py
import os

import scipy.special as spec
from numpy import complex

def erfi(z):
  #print z<0
  return (complex(0.,-1.)*spec.erf(complex(0.,1.)*z)).real
  
  
    
def driving_current(rateWnt,mode='ER',N=10000,NeuronType=2):
  
  # get 3 currents: 
  # balanced one (I = -rate*J0*tau_M)
  # selfconsistent one (c++ program, from firing rate)
  # actual one (read from simulations, e.g. average from statistics)

  #erfi = lambda z: complex(0.,-1.)*spec.erf(complex(0.,1.)*z)	#define imaginary error function

  # define and process network stats
  V_T = I_T = 1.
  V_R = 0.
  J0 = -1
  tauM = 10/1000.

  K = 1000
  #sim_array = [2,5,10,20]
  sim_array = [10,20,50,100,200,500]#,1000,2000,5000]
  
  
  # prepare iteration and processing
  #steps = min(20,len(K_array))
  steps = min(20,len(sim_array))
  
  I_bal = np.zeros(steps)
  I_selfcon = np.zeros(steps)
  I_sim = np.zeros(steps)
  
  V_steps = 1001
  V = np.linspace(-0.5,1,V_steps)
  P_V = np.zeros((steps,V_steps))
  
  fig,axes = plt.subplots(nrows=1,ncols=2)
  
  for i in range(len(sim_array)):
    #rateWnt = sim_array[i]
    K = sim_array[i]
    I_bal[i] = -rateWnt*J0*tauM
    sigma2 = J0**2*rateWnt*tauM
    sigma = np.sqrt(sigma2)
  
    # get selfconsistent I from solution of Fokker-Planck eq. (Brunel, Hakim)
    selfcon_sv = 'data/selfcon_K%d_rateWnt%g.nc' % (K,rateWnt)
    runStr = './calc_I %d %g %s' % (K,rateWnt,selfcon_sv)
    runIt(runStr)
    try:
      ncid = netcdf.netcdf_file(selfcon_sv,'r')
      I_selfcon[i] = ncid.variables['I0'].getValue()
      ncid.close()
      
      os.remove(selfcon_sv)
      
      # if found, plot the resulting distribution
      mu = V_T + np.sqrt(K)*(I_selfcon[i] + J0*rateWnt*tauM)
      print mu
      
      erfi_tmp = np.zeros(V_steps)
      
      for j in range(V_steps):
	if V[j] < 0:
	  erfi_tmp[j] = erfi((V_R-mu)/sigma)
	else:
	  erfi_tmp[j] = erfi((V[j]-mu)/sigma)
	#if V[j] > 0.9:
	  #print erfi_tmp[j]
      P_V[i] = np.sqrt(math.pi)*rateWnt*tauM/sigma*np.exp(-(V-mu)**2/sigma2)*(erfi((V_T-mu)/sigma)-erfi_tmp)
      print "success"
      #plt.show()
    except:
      print "Parameters K=%d and rateWnt=%g gave no solution!" % (K,rateWnt)
      I_selfcon[i] = np.nan
    
    
    # Now, add reading of simulation results (get statistics or smth for these, or just make ratefinding)
    
    try:
      if mode == 'ER':
	analyze = scriptFrame('r','drive_cur',1,alpha=[[0],[0],[0],[0]],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[J0],'NeuronType':[NeuronType]},mode_return=1)
      elif mode == 'drive':
	analyze = scriptFrame('r','drive_cur',1,hPconst=[[0],[0],[0],[0]],netConst={'N':[N],'K':[K],'rateWnt':[rateWnt],'J0':[J0],'NeuronType':[NeuronType]},mode_return=1)
      #print analyze
      mask = np.where(analyze['drive_cur'])
      I_sim[i] = np.mean(analyze['drive_cur'][mask])/np.sqrt(K)
      print I_sim
    except:
      print "Simulation is not yet available"
      I_sim[i] = np.nan
  #print I_sim
  #plt.figure()
  print P_V
  print P_V.shape
  print I_selfcon
  axes[1].plot(V,P_V.transpose())

  axes[0].plot(sim_array,I_bal,'-k',label='balance equation')
  axes[0].plot(sim_array,I_selfcon,'or',label='selfconsistency')
  axes[0].plot(sim_array,I_sim,'ob',label='simulations')
  axes[0].set_xlim([0,max(sim_array)])
  axes[0].legend(loc=2)
  axes[1].set_xlim([-0.5,1.5])
  
  plt.show(block=False)
  
  
  
  #def calc_max_dist(self,I_0=None):
    
    # define and process inputs

    
    #nu_array = [1,2,5,10,20]
    
    #plt.figure()
    #for j in range(5):
    #self.Const['rateWnt'] = nu_array[j]
    
    #if not I_0:		# can be taken from simulation after warmup
      #I_0 = I_T - self.Const['J0']*np.sqrt(self.Const['K'])*self.Const['rateWnt']*self.Const['tauM']/1000.
    #print I_0
      #if self.topo == 'p':
	#I_0 += - self.topoConst['drive_cplg']*self.topoConst['drive_rate']*self.Const['tauM']/1000.	
      #print I_0
    
    
    
    
    #if self.topo == 'p':
      #sigma2 += self.topoConst['drive_cplg']**2*self.topoConst['drive_rate']*self.Const['tauM']/1000.
      #mu += self.topoConst['drive_cplg']*self.topoConst['drive_rate']*self.Const['tauM']/1000.
    
    #sigma = np.sqrt(sigma2)
    
    
    
    
    