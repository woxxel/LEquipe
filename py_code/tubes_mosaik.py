import os,imp

assert os.path.exists('LEquipe'), 'Please run the program in the directory of LEquipe to properly set up the paths.'

directory = os.getcwd()

imp.load_source('read_write_code','py_code/read_write.py')
imp.load_source('file_management_code','py_code/file_management.py')
imp.load_source('support_code', 'py_code/support.py')
imp.load_source('analysis_code', 'py_code/analysis.py')
from read_write_code import *
from file_management_code import *
from support_code import *
from analysis_code import *


############### define the parameters of the network here ###########

def tubes_mosaik(directory = directory, suppressMessages = 1, sv_file = [0,0,0,1]):
  
  ParaNet = {}
  ParaSim = {}
  ParaTopo = {}
  
  steps = 21			#stepnumber for application of different perturbations (has to be uneven, so that (0,0) is included

  ParaNet['NeuronType'] = 2	#LIF neuron
  
  ParaNet['N'] = 1000		#number of neurons N
  ParaNet['K'] = 100 		#number of synapses per neuron K
  ParaNet['J0'] = -1
  ParaNet['tauM'] = 10
  ParaSim['rateWnt'] = 10.	#network-averaged firing rate in Hz
  ParaSim['SR'] = 10
  
  ParaTopo['post'], ParaTopo['row_length'] = random_graph(ParaNet['K'],ParaNet['N'])
  
  ParaSim['saveFinalState'] = 1
  
  directory += '/data/py/'
  if not os.path.exists(directory):
    print 'creating new directory: %s' % directory
    os.mkdir(directory)
    
  HashNet, FileNet, ParaNet = writeNet(directory, ParaNet)
  HashTopo, FileTopo = writeTopo(directory, ParaNet, ParaTopo)
  HashSim, FileSim = writeSim(directory, ParaNet, ParaSim)
  
  HashDataOut = hashlib.sha1(np.array([HashNet, HashTopo, HashSim])).hexdigest()
  
  FileOut = directory + 'DataOut-' + HashDataOut + '.nc'

  ## run the C++ simulation
  os.system('./LEquipe ' + FileNet + ' ' + FileTopo + ' ' + FileSim + ' ' + FileOut)	#get the unperturbed trajectory
  
  Data = readDataOut(FileOut)
  base_state = Data['finalStates'][:,0]
  ParaNet['Iext'] = Data['finalCurrents'][0]
    
  ## find two perpendicular vectors
  pert_vec = np.zeros((2,ParaNet['N'][0]))
  norm = 1
  while 1:
    for i in range(2):
      vec = np.random.uniform(-1,1,ParaNet['N'])
      pert_vec[i] = vec/np.linalg.norm(vec)
    norm = abs(np.sum(pert_vec[0]*pert_vec[1]))
    if norm < 0.1:
      break  
  
  # set limits of perturbations
  bounds = 0.1
  x_range = np.linspace(-bounds,bounds,steps)
  y_range = np.linspace(-bounds,bounds,steps)
  
  state_vec = np.zeros((steps,steps,ParaNet['N'][0]))
  ParaSim['SW'] = 0		# disable ...
  ParaSim['TW'] = 0		# ... warmup
  
  ParaSim['SR'] = 0		# disable ...
  ParaSim['TR'] = 0		# ... rate finding
  
  HashSim, FileSim = writeSim(directory, ParaNet, ParaSim)
  
  for i in range(steps):
    for j in range(steps):
      print 'stepnumber (%d,%d)', (i,j)
      ParaNet['init'] = base_state + x_range[i]*pert_vec[0]+y_range[j]*pert_vec[1]	# get initial state for perturbed trajectory
     
      ## write all parameters to netcdf files to directory data/ and get the hashes
      HashNet, FileNet, ParaNet = writeNet(directory, ParaNet)
      HashDataOut = hashlib.sha1(np.array([HashNet, HashTopo, HashSim])).hexdigest()
      FileOut = directory + 'DataOut-' + HashDataOut + '.nc'

      os.system('./LEquipe ' + FileNet + ' ' + FileTopo + ' ' + FileSim + ' ' + FileOut)

      ## read the output file and plot the results
      Data = readDataOut(FileOut)

      state_vec[i,j] = get_finalState(FileOut)

  
  tubes, dist = find_tubes(state_vec,steps)
  save_tubes(tubes)
  analyse_tubes(tubes,bounds,ParaNet['N'][0],ParaNet['K'][0],ParaSim['rateWnt'][0])
  
  plt.pcolor(x_range,y_range,tubes)
  plt.xlim([-bounds,bounds])
  plt.ylim([-bounds,bounds])
  plt.show()

  if not sv_file[0]:
    os.remove(FileNet)
  if not sv_file[1]:
    os.remove(FileTopo)
  if not sv_file[2]:
    os.remove(FileSim)
  if not sv_file[3]:
    os.remove(FileOut)

  
  #return state_vec,tubes,dist

tubes_mosaik()