import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def bootstrap(X,N,lower=0.025,upper=0.975,n=None):
  
  assert len(X.shape) <= 2, "Please only use 1-dim arrays (or 1 dim with error) for bootstrapping"
  
  if len(X.shape) == 1:
    mask = np.invert(np.isnan(X))
    if not np.sum(mask):
      return [np.nan]*3
    X_masked = X[mask]
    
    bootstrap_distr = np.zeros(N)
    
    for i in range(N):
      bootstrap_distr[i] = np.mean(bootstrap_resample(X_masked))
  
  elif len(X.shape) == 2:
    mask_nan = np.invert(np.isnan(X[0]))
    mask_inf = np.invert(np.isinf(X[1]))
    if not np.sum(mask_nan&mask_inf):
      return np.array([np.nan,np.nan,np.nan])
    
    X_masked = X[:,mask_nan&mask_inf]
    
    bootstrap_distr = np.zeros(N)
    
    for i in range(N):
      resample_X = bootstrap_resample(X_masked)
      bootstrap_distr[i] = np.average(resample_X[0],weights=1./resample_X[1])
  
  bootstrap_distr = np.sort(bootstrap_distr)
  
  #plt.figure()
  #plt.hist(bootstrap_distr,bins=100)
  #plt.show(block=False)
  
  lower_bound = bootstrap_distr[int(lower*N)]
  upper_bound = bootstrap_distr[int(upper*N)]
  mean = np.mean(bootstrap_distr)
  
  return np.array([mean, lower_bound, upper_bound])


def bootstrap_resample(X,n=None):
  
  if n == None:
    n = X.shape[-1]
  
  resample_i = np.floor(np.random.randint(n,size=n)).astype(int)
  
  if len(X.shape) == 1:  
    X_resample = X[resample_i]
  
  elif len(X.shape) == 2:
    X_resample = X[:,resample_i]
  
  return X_resample


def plt_bootstrap(X,Y,ax,col,ls='-',ms='o',label=None,fit_func=None,p0=None,mask=None):
  
  assert len(X.shape) == 1, 'Please specify a one dimensional range array!'
  assert Y.shape[0] == X.shape[0], 'X and Y arrays should have the same length!'
  assert Y.shape[1] == 3, 'Y array should include mean value and upper and lower bound of confidence intervall!'

  if mask == None:
    mask = np.ones(len(X)).astype('bool')

  ax.plot(X[mask],Y[mask,0],'-',color=col,linestyle=ls,marker=ms,label=label,linewidth=2)
  
  ax.fill_between(X[mask],Y[mask,1],Y[mask,2],color=col,linestyle=ls,alpha=0.2,edgecolor=None,linewidth=0)
  
  if fit_func:
    Y_sigma = np.max(Y[:,1:] - Y[:,0].reshape(len(Y),1),axis=1)
    
    popt,pcov = curve_fit(fit_func,X[mask],Y[mask,0],sigma=Y_sigma[mask],p0=p0)
    
    perr = np.sqrt(np.diag(pcov))
    
    print 'fit results: ', popt
    print 'fit errors: ', perr
    ax.plot(X,fit_func(X,popt[0],popt[1]),'--',color='r',linewidth=0.5)
    
    return popt,perr