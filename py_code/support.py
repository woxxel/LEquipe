import numpy as np
import hashlib, math, os
  

def random_graph(N,K):
  
  topoConstruct = {}
  
  A = (np.random.random((N,N)) < K/(N-1.))*np.invert(np.eye(N).astype('bool')).astype('int')
  
  #topoConstruct['matrix'] = A
  
  topoConstruct['inDegree'] = np.sum(A,axis=0)
  topoConstruct['outDegree'] = np.sum(A,axis=1)
  
  topoConstruct['preSyn'] = np.where(A.transpose())[1]
  topoConstruct['postSyn'] = np.where(A)[1]
  
  #edges = np.sum(A)
  
  #c1 = N*(N-1.)
  #c2 = N*(N-1)*(N-2)/2.
  
  #topoConstruct['p_hat'] = edges/c1
  
  #p2 = topoConstruct['p_hat']**2

  #A_square = np.dot(A,A)
  
  #N_recip = np.trace(A_square)/2.
  #N_conv = (np.sum(np.dot(A.transpose(),A)) - edges)/2.
  #N_div = (np.sum(np.dot(A,A.transpose())) - edges)/2.
  #N_chain = np.sum(A_square) - 2*N_recip
    
  #p2_recip = N_recip/(c1/2.)
  #p2_conv = N_conv/c2
  #p2_div = N_div/c2
  #p2_chain = N_chain/(c2*2.)
  
  #topoConstruct['alpha_recip_hat'] = p2_recip/p2 - 1
  #topoConstruct['alpha_conv_hat'] = p2_conv/p2 - 1
  #topoConstruct['alpha_div_hat'] = p2_div/p2 - 1
  #topoConstruct['alpha_chain_hat'] = p2_chain/p2 - 1
  
  return topoConstruct


def postsynapticvector2adjacencymatrix(A,B):
  size = int(len(B))
  conn = np.zeros((size,size))
  start_idx = 0
  for i in range(len(B)):
    idx = [int(numb) for numb in A[start_idx:start_idx+B[i]]]
    conn[i,idx] = 1
    start_idx += B[i]
  return conn
  
  
def adjacencymatrix2postsynapticvector(A):
  N = len(A)
  HomogSynapse = 1
  
  row_length = []
  post = []

  if len(np.unique(A)) > 2:
    HomogSynapse = 0
  
  neuron_idx = np.arange(N)
  for n in range(N):
    connected = (A[n,:] != 0)
    connected[n] = 0
    
    row_length.append(sum(connected))
    post.extend(neuron_idx[connected])
  
  return post, row_length
  
  
  
def set_default(Para,check_val,def_val,suppressMessages=1):
  if not (check_val in Para.keys()):
    if (type(def_val) == list):
      Para[check_val] = def_val
    else:
      Para[check_val] = [def_val]
    if not suppressMessages:
      if (def_val == []):
	print '%s is set to default value %s = []' % (check_val,check_val)
      else:
	print '%s is set to default value %s = %g' % (check_val,check_val,def_val)


def check_hetero(check_val,HomogNetwork,suppressMessages=1):
  if (len(check_val[1]) > 1):
    if not suppressMessages:
      print '%s is heterogeneous' % check_val[0]
    return 0
  else:
    return HomogNetwork
    
    
def network_check(Para,N,HomogNetwork=1,suppressMessages=1):
  for i in range(len(Para.keys())):
    if type(Para.values()[i])!=list:	#transform to arrays
      if type(Para.values()[i]) in [int,float,np.float64,np.int64]:
	Para[Para.keys()[i]] = [Para.values()[i]]
      elif type(Para.values()[i]) == np.ndarray:
	try:
	  if len(Para.values()[i]) > 1:
	    Para[Para.keys()[i]] = [item for item in Para.values()[i]]
	  else:
	    Para[Para.keys()[i]] = [Para.values()[i][0]]
	except:
	  Para[Para.keys()[i]] = [Para.values()[i]]
	length = len(Para.values()[i])
	if Para.keys()[i] != 'init':
	  HomogNetwork = check_hetero(Para.items()[i],HomogNetwork)
      elif type(Para.values()[i]) == str:
	Para[Para.keys()[i]] = list(Para.values()[i])
      elif type(Para.values()[i]) == dict:
	Para.values()[i], HomogNetwork = network_check(Para.values()[i],N,HomogNetwork)
      elif type(Para.values()[i]) != dict:
	assert 0, "%s has bad type %s" % (Para.keys()[i],type(Para.values()[i]))
    else:
      if Para.keys()[i] != 'init':
	HomogNetwork = check_hetero(Para.items()[i],HomogNetwork)
  return Para, HomogNetwork


def hashing(Para,values):
  hash_list = []
  for i in range(len(values)):
    if not values[i] in Para.keys():
      #print 'Warning: %s not specified' % values[i]
      hash_list.extend([])
    elif type(Para[values[i]]) == dict:
      for j in range(len(Para[values[i]])):
	hash_list.extend(Para[values[i]].values()[j])
    else:
      try:
	if type(Para[values[i]][0]) == list:
	  for item in Para[values[i]]:
	    hash_list.extend(item)
	else:
	  hash_list.extend(Para[values[i]])
      except:
	hash_list.extend(Para[values[i]])
  Hash = hashlib.sha1(np.array(hash_list)).hexdigest()
  return Hash