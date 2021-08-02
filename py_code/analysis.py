import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

def prob_convergence(Path,precision=0.01,conv_intervall=5,steps=None):
  
  ncid = netcdf.netcdf_file(Path,'r',mmap=False)
  div_track = ncid.variables['div_track'][:]
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  # this is testing the optimal number of steps, based on the achieved convergence value with all simulations
  if not steps:
    ps_steps = min(ncid.dimensions['ps_steps'],np.sum(ncid.variables['p'][:]>0))
    co_steps = ncid.dimensions['co_steps']
    in_steps = ncid.dimensions['in_steps']
    pd_steps = ncid.dimensions['pd_steps']
  
  print ps_steps
  x_array = np.arange(pd_steps)+1
  conv_idx_pd = np.zeros((ps_steps,co_steps,in_steps))
  conv_value_pd = np.zeros((ps_steps,co_steps,in_steps))
  trace = np.zeros((co_steps,in_steps,pd_steps))
  trace_mean = np.zeros((ps_steps,pd_steps))
  # always get last entry and check, when simulation stays within precision of 0.01 for 5+ simulations
  non_conv = 0
  
  
  plt.figure(figsize=(4,3),dpi=80)
  
  ax = plt.subplot(111)
  ax.set_position([0.2,0.17,0.75,0.75])
  
  for ps_ct in range(ps_steps):
    col = 2*float(ps_ct+1)/(3*ps_steps)
    for co_ct in range(co_steps):
      for in_ct in range(in_steps):
	conv_trace = np.cumsum(div_track[ps_ct][co_ct][in_ct][:pd_steps])/x_array.astype('float')
	conv_value_pd[ps_ct][co_ct][in_ct] = conv_trace[-1]
	trace[co_ct][in_ct][:] = (conv_trace - conv_value_pd[ps_ct][co_ct][in_ct])**2
	conv_ct = 0
	if conv_value_pd[ps_ct][co_ct][in_ct] in [0,1]:
	  conv_ct = conv_intervall + 1
	  conv_idx_pd[ps_ct][co_ct][in_ct] = -1
	else:
	  
	  for pd_ct in range(pd_steps):
	    if abs(conv_trace[pd_ct] - conv_value_pd[ps_ct][co_ct][in_ct]) < precision:
	      conv_ct += 1
	    else:
	      conv_ct = 0
	    if conv_ct > conv_intervall:
	      conv_idx_pd[ps_ct][co_ct][in_ct] = pd_ct - conv_intervall + 1
	      break
	if conv_ct <= conv_intervall:
	  conv_idx_pd[ps_ct][co_ct][in_ct] = np.nan
	  non_conv += 1
    trace_mean[ps_ct][:] = np.sqrt(np.nansum(np.nansum(trace,axis=1),axis=0)/(co_steps*in_steps))
    plt.plot(x_array,trace_mean[ps_ct],color=(col,col,col),label=r'$\displaystyle \varepsilon = %5.3g$' % ncid.variables['pert_size'][:][ps_ct])
    #time.sleep(20)
  
  plt.plot([100,100],[10**(-3),10**0],'r')
  plt.xscale('log')
  plt.yscale('log')
  
  plt.xlim([1,pd_steps])
  plt.ylim([10**(-3),10**0])
  
  plt.xlabel('\# perturbation directions',fontsize=12)
  plt.ylabel(r'STD($\displaystyle p(\varepsilon)$)',fontsize=12)
  #plt.legend(loc=3,prop={'size':12})
  plt.tick_params(labelsize=10)
  save_str = './pics/SONET/fluxtubes/convergence_pd.pdf'
  plt.savefig(save_str)
  print "Figure saved to %s." % save_str
  
  plt.show(block=False)
  
  print 'non-converged initial conditions: %d' % non_conv
  mask = conv_idx_pd > 0
  print "Fraction of converged simulations within %d perturbation directions: %5.3f" % (pd_steps,1-np.sum(np.isnan(conv_idx_pd))/float(ps_steps*co_steps*in_steps))# this should be > 0.9
  mean = np.nanmean(conv_idx_pd[mask])
  var = np.nanvar(conv_idx_pd[mask])
  print "mean steps to convergence: %3.1f" % mean
  print "coefficient of variation of steps: %5.3f\n" % np.sqrt(var/mean**2)
  
  
  
  
  conv_idx_in = np.zeros((ps_steps,co_steps))
  conv_value_in = np.zeros((ps_steps,co_steps))
  
  trace = np.zeros((co_steps,in_steps))
  trace_mean = np.zeros((ps_steps,in_steps))
  
  x_array = np.arange(np.sum(in_steps)) + 1
  
  non_conv = 0
  
  
  
  plt.figure(figsize=(4,3),dpi=80)
  
  ax = plt.subplot(111)
  ax.set_position([0.2,0.17,0.75,0.75])
  
  for ps_ct in range(ps_steps):
    col = 2*float(ps_ct+1)/(3*ps_steps)
    for co_ct in range(co_steps):
      mask = np.invert(np.isnan(conv_value_pd[ps_ct][co_ct][:in_steps]))
     
      conv_trace = np.cumsum(conv_value_pd[ps_ct][co_ct][:in_steps])/x_array.astype('float')
      
      conv_value_in[ps_ct][co_ct] = conv_trace[-1]
      trace[co_ct][:] = np.nan
      trace[co_ct][:len(conv_trace)] = (conv_trace - conv_value_in[ps_ct][co_ct])**2
      x_array1 = np.arange(np.sum(mask)) + 1
      conv_trace = np.cumsum(conv_value_pd[ps_ct][co_ct][:in_steps][mask])/x_array1.astype('float')
      conv_value_in[ps_ct][co_ct] = conv_trace[-1]
      if conv_value_in[ps_ct][co_ct] in [0,1]:
	conv_idx_in[ps_ct][co_ct] = -1
      else:
	conv_ct = 0
	for in_ct in range(in_steps):
	  if abs(conv_trace[in_ct] - conv_value_in[ps_ct][co_ct]) < precision:
	    conv_ct += 1
	  else:
	    conv_ct = 0
	  if conv_ct > conv_intervall:
	    conv_idx_in[ps_ct][co_ct] = in_ct - conv_intervall + 1
	    break
      if conv_ct <= conv_intervall:
	conv_idx_in[ps_ct][co_ct] = np.nan
	non_conv += 1
    trace_mean[ps_ct][:] = np.sqrt(np.nansum(trace,axis=0)/co_steps)
    plt.plot(x_array,trace_mean[ps_ct],color=(col,col,col))
  
  plt.plot([25,25],[10**(-3),10**0],'r')
  
  plt.xscale('log')
  plt.yscale('log')
  
  plt.xlim([1,in_steps])
  plt.ylim([10**(-3),10**0])
  
  plt.xlabel('\# initial conditions',fontsize=12)
  plt.ylabel(r'STD($\displaystyle p(\varepsilon)$)',fontsize=12)
  #plt.legend(loc=3,prop={'size':12})
  plt.tick_params(labelsize=10)
  save_str = './pics/SONET/fluxtubes/convergence_in_N=.pdf'
  plt.savefig(save_str)
  print "Figure saved to %s." % save_str
  
  plt.show(block=False)
  
  print 'non-converged initial conditions: %d' % non_conv
  mask = conv_idx_in > 0
  print "Fraction of converged simulations within %d initial conditions: %5.3f" % (in_steps,1-np.sum(np.isnan(conv_idx_in))/float(ps_steps*co_steps))# this should be > 0.9
  mean = np.nanmean(conv_idx_in[mask])
  var = np.nanvar(conv_idx_in[mask])
  print "mean steps to convergence: %3.1f" % mean
  print "coefficient of variation of steps: %5.3f\n" % np.sqrt(var/mean**2)










  conv_idx_co = np.zeros((ps_steps))
  conv_value_co = np.zeros((ps_steps))
  
  trace = np.zeros((co_steps))
  trace_mean = np.zeros((ps_steps,co_steps))
  
  x_array = np.arange(co_steps) +1
  
  plt.figure(figsize=(4,3),dpi=80)
  
  ax = plt.subplot(111)
  ax.set_position([0.2,0.17,0.75,0.75])
  
  for ps_ct in range(ps_steps):
    col = 2*float(ps_ct+1)/(3*ps_steps)
    conv_trace = np.cumsum(conv_value_in[ps_ct][:co_steps])/x_array.astype('float')
    conv_value_co[ps_ct] = conv_trace[-1]
    trace[:len(conv_trace)] = (conv_trace - conv_value_co[ps_ct])**2
    
    if conv_value_co[ps_ct] in [0,1]:	#should never, ever happen for normal perturbations
      conv_idx_co[ps_ct] = -1
    else:
      conv_ct = 0
      for co_ct in range(co_steps):
	if abs(conv_trace[co_ct] - conv_value_co[ps_ct]) < precision:
	  conv_ct += 1
	else:
	  conv_ct = 0
	if conv_ct > conv_intervall:
	  conv_idx_co[ps_ct] = co_ct - conv_intervall + 1
	  break
    if conv_ct <= conv_intervall:
      conv_idx_co[ps_ct] = np.nan
    plt.plot(x_array,trace,color=(col,col,col))
  
  plt.plot([5,5],[10**(-5),10**(-2)],'r')
  
  #plt.xscale('log')
  plt.yscale('log')
  plt.tick_params(labelsize=10)
  plt.xlim([1,co_steps])
  plt.ylim([10**(-5),10**(-2)])
  
  plt.xlabel('\# connectivities',fontsize=12)
  plt.ylabel(r'STD($\displaystyle p(\varepsilon)$)',fontsize=12)
  
  save_str = './pics/SONET/fluxtubes/convergence_co.pdf'
  plt.savefig(save_str)
  print "Figure saved to %s." % save_str

  plt.show(block=False)
  mask = conv_idx_co > 0
  print "Fraction of converged simulations within %d connectivities: %5.3f" % (co_steps,1-np.sum(np.isnan(conv_idx_co))/float(ps_steps))# this should be 1
  mean = np.nanmean(conv_idx_co[mask])
  var = np.nanvar(conv_idx_co[mask])
  print "mean steps to convergence: %3.1f" % mean
  print "coefficient of variation of steps: %5.3f\n" % np.sqrt(var/mean**2)

  #co_conv = 5
  #in_conv = 25
  #pd_conv = 100
  
  #plt.figure()
  #steps_conv = co_conv*in_conv*pd_conv
  #for ps_ct in range(ps_steps):
    ##print div_track[ps_ct][:3,:25,:100].shape
    ##mask = np.invert(np.isnan(div_track[ps_ct][:co_conv,:in_conv,:pd_conv]))
    #p_conv = 1 - np.cumsum(div_track[ps_ct][:co_conv,:in_conv,:pd_conv])/(np.arange(steps_conv)+1).astype('float')
    
    #plt.plot(np.arange(steps_conv)+1,p_conv)
    #mask = np.invert(np.isnan(p_conv))
    ##print 1-p_conv
    #print 1-p_conv[mask][-1]
  #plt.show(block=False)
  #print '\n'