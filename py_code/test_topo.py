import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os, sys, imp
imp.load_source('support_code', 'py_code/support.py')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from support_code import *

def topo_generate(N=1000,save=0):
  topo('conv',[0,2],1000,5,save)
  topo('div',[0,4],1000,5,save)
  topo('recip',[0,8],1000,5,save)
  topo('chain',[-0.5,0.5],1000,5,save)
  




def topo(alpha,a_range,N=1000,steps=5,save=0):
  
  K = 0.1*N
  J0 = -1
  const = {}
  const['alpha_recip']=0
  const['alpha_conv']=0
  const['alpha_div']=0
  const['alpha_chain']=0
  
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  
  ccmap = mcolors.LinearSegmentedColormap.from_list(name='custom',colors=[(0,(0.8,0.8,0.8)),(1,(0.05,0.05,0.05))],N=steps)
  cbar_dummy = dummy_plot(ccmap,steps,steps+1)
  
  if steps == 1:
    ax0_border = 2*K
    ax1_border = 2*K
  else:
    if alpha == 'conv':
      ax0_border = 2*K
      ax1_border = N
    elif alpha == 'div':
      ax0_border = N
      ax1_border = 2*K
    elif alpha == 'recip':
      ax0_border = 2*K
      ax1_border = 2*K
    elif alpha == 'chain':
      ax0_border = N
      ax1_border = N
      
  if alpha == 'chain':
    const['alpha_conv']=0.5
    const['alpha_div']=1
  alpha_array = np.linspace(a_range[0],a_range[1],steps)
  fig, axes = plt.subplots(nrows=1,ncols=2)
  for a in range(steps):
    print alpha_array[a]
    const['alpha_' + alpha] = alpha_array[a]
    run_str = './SONETs/run_secorder '
    alpha_str = '%5.3f %5.3f %5.3f %5.3g' % (const['alpha_recip'],const['alpha_conv'],const['alpha_div'],const['alpha_chain'])
    run_str += str(N) + ' ' + str(K/float(N)) + ' ' + alpha_str
    #run_str += ' ' + seed
    
    #try:
    devnull = open('/dev/null', 'w')
    oldstdout_fno = os.dup(sys.stdout.fileno())
    os.dup2(devnull.fileno(), 1)
    os.system(run_str)			# run SONet creation with suppressed output
    devnull.close()
    os.dup2(oldstdout_fno, 1)
    os.close(oldstdout_fno)
    
    fileParas = '%d_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f' % (N,K/float(N),const['alpha_recip'],const['alpha_conv'],const['alpha_div'],const['alpha_chain'])
    fileName = 'data/w_' + fileParas + '.dat'
    A = np.zeros((N,N))
    idx = 0
    
    topo_rd = open(fileName,'r')
    for line in topo_rd:
      A[idx,:] = line.split(' ')[:N]
      idx += 1
    post, row_len = adjacencymatrix2postsynapticvector(A.transpose())
    if len(np.unique(np.array(J0))) > 2:
      J = np.ones(len(post))*J0/math.sqrt(K)
    else:
      J = J0/np.sqrt(K)
    topo_rd.close()
    
    # clean up
    for dat in os.listdir('data/'):
      if fileParas in dat:
	os.remove('data/' + dat)
    
    if steps == 1:
      col_plot = 0
    else:
      norm = 4./3*(steps-1)
      col_plot = 0.8-a/norm
    

    
    in_degree = [post.count(i) for i in range(N)]
    out_degree = row_len
    
    in_out_corr = np.mean([(in_degree[i]-K)*(out_degree[i]-K) for i in range(N)])/np.sqrt(np.var(in_degree)*np.var(out_degree))
    #print "correlation: %5.3f" % in_out_corr
    
    #return in_degree, out_degree
    #return post, row_len
    ### plot out degree
    #data = np.histogram(row_len,bins=2*K,range=(0,2*K))
    #ax0_border = max(np.max(row_len),2*K)
    axes[0].hist(out_degree,bins=50,range=(0,ax0_border),color=(col_plot,col_plot,col_plot),alpha=0.8,edgecolor="none")
    #ax1_border = max(np.max(in_degree),2*K)
    axes[1].hist(in_degree,bins=50,range=(0,ax1_border),color=(col_plot,col_plot,col_plot),alpha=0.8,edgecolor="none")
    #plt.plot(data[1][:-1],data[0],color=(a/(2.*steps),a/(2.*steps),0))
    #plt.show()
    #except:
      #print 'no topology found'
  axes[0].set_xlabel('Out-degree',fontsize=18)
  axes[0].set_ylabel('Probability',fontsize=18)
  axes[1].set_xlabel('In-degree',fontsize=18)
  axes[0].set_xlim([0,ax0_border])
  axes[1].set_xlim([0,ax1_border])
  axes[0].set_ylim([0,200])
  axes[1].set_ylim([0,200])
  
  axes[0].set_xticks(np.linspace(0,ax0_border,5))
  axes[0].set_xticklabels(np.linspace(0,ax0_border,5),fontsize=14)
  axes[1].set_xticks(np.linspace(0,ax1_border,5))
  axes[1].set_xticklabels(np.linspace(0,ax1_border,5),fontsize=14)
  
  axes[0].set_yticks(np.linspace(0,200,5))
  axes[0].set_yticklabels(np.linspace(0,0.2,5),fontsize=14)
  axes[1].set_yticks(np.linspace(0,200,5))
  axes[1].set_yticklabels(np.linspace(0,0.2,5),fontsize=14)
  
  
  if steps > 1:
    divider0 = make_axes_locatable(axes[1])
    cax0 = divider0.append_axes("right", size="5%", pad=0.2)
    cbar0 = plt.colorbar(cbar_dummy, cax=cax0)
    cbar0.ax.set_ylabel(r'$\displaystyle\alpha_{%s}$'%alpha,fontsize=22)
    
    cbar0.set_ticks(np.linspace(1./2,steps-1./2,steps))
    cbar0.set_ticklabels(np.linspace(min(alpha_array),max(alpha_array),steps))
    save_str = './pics/SONET/in_out_degree_%s.png' % alpha
  else:
    save_str = './pics/SONET/in_out_degree_%s=%3.1f.png' % (alpha,alpha_array[0])
  if save:
    plt.savefig(save_str)
  #fig.suptitle(r'$\displaystyle\alpha_{%s}$'%alpha,fontsize=20)

  plt.show(block=False)
  
  #return in_degree
  

def dummy_plot(ccmap,border,steps):
  plt.figure()
  # Using contourf to provide my colorbar info, then clearing the figure
  Z = [[0,0],[0,0]]
  levels = np.linspace(0,border,steps)
  cbar_dummy = plt.contourf(Z, levels, cmap=ccmap)
  plt.clf()
  return cbar_dummy