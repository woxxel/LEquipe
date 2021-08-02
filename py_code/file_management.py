import os, sys, time, imp, subprocess

imp.load_source('read_write_code','Code/py_code/read_write.py')
from read_write_code import *


def wait_for_godot(fileName):
  # check if initial Data is provided. if not, wait for it
  wait = 0
  waited = 0
  while (wait == 0):      
    try:
      Data = readDataOut(fileName)
      wait = 1
    except:
      if waited == 0:
	print '\n waiting for input file %s ...' % fileName
	waited = 1
      time.sleep(10)
  return Data


def runIt(cluster,run_str,path=None,name=None,hide=1):
  
  clusterstr = ''
  if cluster:
    if cluster['q']:
      if 'name' in cluster.keys():
	clustername = cluster['name']
      else:
	clustername = 'frigg'
      clusterstr = 'qsub -N %s -e %s -o %s -q %s.q -cwd -b y -V ' % (name,path,path,clustername)
  
  if hide:
    devnull = open('/dev/null', 'w')
    oldstdout_fno = os.dup(sys.stdout.fileno())
    os.dup2(devnull.fileno(), 1)
  os.system(clusterstr + run_str)		# run command with suppressed output
  if hide:
    devnull.close()
    os.dup2(oldstdout_fno, 1)
    os.close(oldstdout_fno)
    

def print_status(sim,time_start,max_sim=1100,wait_time=10,pertSize=None):
  # give output with status & time remaining
  time_now = time.time()-time_start
  
  if 'ps_ct' in sim.steps.keys():
    steps = sim.steps['ct'] - sim.steps['ps_ct']*sim.steps['total']
  else:
    steps = sim.steps['ct']
  
  if 'pd_ct' in sim.steps.keys():
    norm = sim.steps['total']
    print_str_1 = '%d. perturbation step (ps=%7.5g): \t(co,in,pd): (%d,%d,%d) of (%d,%d,%d) '% (sim.steps['ps_ct']+1,pertSize,sim.steps['co_ct']+1,sim.steps['in_ct']+1,sim.steps['pd_ct']+1,sim.steps['co'],sim.steps['in'],sim.steps['pd'])
  elif pertSize:
    norm = sim.steps['total']
    print_str_1 = '%d. perturbation step (ps=%7.5g): \t(co,in,pd): (%d,%d,:) of (%d,%d,%d) '% (sim.steps['ps_ct']+1,pertSize,sim.steps['co_ct']+1,sim.steps['in_ct']+1,sim.steps['co'],sim.steps['in'],sim.steps['pd'])
  elif 'in_ct' in sim.steps.keys():
    norm = sim.steps['co']*sim.steps['in']
    print_str_1 = 'Warming up the network (co,in): (%d,%d) of (%d,%d) '% (sim.steps['co_ct']+1,sim.steps['in_ct']+1,sim.steps['co'],sim.steps['in'])
  else:
    norm = sim.steps['co']
    print_str_1 = 'Ratefinding for topologies: (%d) of (%d) '% (sim.steps['co_ct']+1,sim.steps['co'])
    
  ratio = steps/float(norm)
  print_str_2 = '\t%d%% done, about %d seconds left [' % (int(100*ratio),int(time_now*(1./ratio-1))) + '='*int(20*ratio)
  
  print_str_3 = ''
  
  if sim.cluster['q']:
    p = subprocess.Popen('qstat -u aschmidt | wc -l',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
    out,err = p.communicate()
    workload = int(out)
  
    if workload > max_sim:
      print_str_3 += '\t (queue too full - waiting...)'
      while (workload > max_sim): #check workload
	p = subprocess.Popen('qstat -u aschmidt | wc -l',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
	out,err = p.communicate()
	workload = int(out)
	time.sleep(wait_time)
    
  sys.stdout.write('\r{0}'.format(print_str_1 + print_str_2 + print_str_3))
  sys.stdout.flush()


def run_script(Path,sim,CDB,hide):
  
  q = sim.cluster['q']
  
  runstr= ' ./Code/lequipe/LEquipe '
  
  if sim.state:# and sim.mode in ['eft','eft_tmp']:
    jobName = 'run_%02d_%02d' % (sim.steps['ps_ct'],sim.steps['co_ct'])
    if 'pd_ct' in sim.steps.keys():
      if (sim.steps['in_ct'] == 0) and (sim.steps['pd_ct'] == 0):
	# remove old file
	runscrID = open(Path['results'] + jobName,'w')
	runscrID.close()
    else:
      if (sim.steps['in_ct'] == 0):
	runscrID = open(Path['results'] + jobName,'w')
	runscrID.close()
    
  #if not 'Ref' in Path.keys():
    #PathRef = ''
  #else:
    #PathRef = Path['Ref']
  if not os.path.exists(Path['Out']):	# no parallel mode possible like that
    cmd = runstr + Path['Net'] + ' ' + Path['Topo'] + ' ' + Path['Sim'] + ' ' + Path['Out'] + ' ' + Path['Puppet']
    if sim.cluster['q']:
      if sim.state:#(sim.state and (sim.mode in ['eft','eft_tmp'])):
	cmd += '\n'
	
	# write it to file
	runscrID = open(Path['results'] + jobName,'a')
	runscrID.write(cmd)
	runscrID.close()
	
	if sim.steps['in_ct'] == sim.steps['in']-1:
	  batch_cmd = 'python ./Code/py_code/submit_batch.py %s' % (Path['results'] + jobName)
	  runIt(sim.cluster,batch_cmd,Path['txt'],Path['ID'],hide)
      else:
	runIt(sim.cluster,cmd,Path['txt'],Path['ID'],hide)
    else:
      runIt(sim.cluster,cmd,Path['txt'],Path['ID'],hide)
      #cmd = runstr + Path['Net'] + ' ' + Path['Topo'] + ' ' + Path['Sim'] + ' ' + Path['Out'] + ' ' + Path['Puppet'] + ' ' + PathRef
      #runIt(q,cmd,Path['txt'],hide)
      
#def run_script_cluster(Path,sim,CDB,loop,hide):
  
  #bk = 2
  #nP=bk*12	#if par=1, how many processors to use, nanna: 12 cores/node,12 nodes/block (fast connections only with block=> max 144 cores)
  
  #q = sim.cluster['q']
  
  #if loop == 'pd':
    #jobName = 'P_%02d_%02d_%02d' % (sim.steps['ps_ct'],sim.steps['co_ct'],sim.steps['in_ct'])
  #else:
    #jobName = 'W_%02d' % (sim.steps['co_ct'])
    
  #if sim.steps['%s_ct'%loop] == 0:	# remove old file
    #runscrID = open(Path['results'] + jobName,'w')
    #runscrID.close()
  
  ## construct run command
  #runstr= './Code/lequipe/LEquipe '
  #if not 'ID' in Path.keys():
    #Path['ID'] = ''  
  #if ((not 'Ref' in Path.keys()) or (not CDB)):
    #PathRef = ''
  #else:
    #PathRef = Path['Ref']
  #if not os.path.exists(Path['Out']):	# no parallel mode possible like that
    ##print sim.cluster
    #cmd = runstr + Path['Net'] + ' ' + Path['Topo'] + ' ' + Path['Sim'] + ' ' + Path['Out'] + ' ' + Path['Puppet'] + ' ' + PathRef
    #if (sim.cluster['q'] and sim.cluster['par']):
      #jobName = 'LEquipe'
      #runscrID = open(jobName,'w')
      #os.system('chmod u+x ' + jobName)
      #runscrID.write('#!/usr/bin/tcsh\n' + \
	      #'#$ -V\n' + \
	      #'#$ -N ' + jobName + '\n' + \
	      #'#$ -e ' + Path['txt'] + ' -o ' + Path['txt'] + '\n' + \
	      #'#$ -cwd\n' + \
	      #'#$ -b y\n' + \
	      #'#$ -pe mvapich2-nanna01 ' + str(nP) + '\n' + \
	      #'mpiexec -machinefile $TMPDIR/machines' + ' -n $NSLOTS ' + cmd + '\n'\
	      #'exit 0\n')
      #runscrID.close()
      #if hide:
	#devnull = open('/dev/null', 'w')
	#oldstdout_fno = os.dup(sys.stdout.fileno())
	#os.dup2(devnull.fileno(), 1)
      #os.system('qsub ./' + jobName)		# run SONet creation with suppressed output
      #if hide:
	#devnull.close()
	#os.dup2(oldstdout_fno, 1)
	#os.close(oldstdout_fno)
      ##runIt('qsub ./' + jobName,hide)
    #else:
      #cmd += '\n'
      
      ## write it to file
      #runscrID = open(Path['results'] + jobName,'a')
      #runscrID.write(cmd)
      #runscrID.close()
      
      #if sim.steps['%s_ct'%loop] == sim.steps[loop]-1:
	#batch_cmd = 'python ./Code/py_code/submit_batch.py %s' % (Path['results'] + jobName)
	#runIt(sim.cluster,Path['txt'],batch_cmd,Path['ID'],hide)
  
def run_script_stats(Path,sim,CDB,hide):
  
  bk = 2
  nP=bk*12	#if par=1, how many processors to use, nanna: 12 cores/node,12 nodes/block (fast connections only with block=> max 144 cores)
  
  q = sim.cluster['q']
  
  jobName = 'Stats_%d_%d' % (sim.steps['co_ct'],sim.steps['in_ct'])
  
  # construct run command
  runstr= './Code/lequipe/LEquipe '
  if not 'ID' in Path.keys():
    Path['ID'] = ''  
  if ((not 'Ref' in Path.keys()) or (not CDB)):
    PathRef = ''
  else:
    PathRef = Path['Ref']
  if not os.path.exists(Path['Out']):	# no parallel mode possible like that
    #print sim.cluster
    cmd = runstr + Path['Net'] + ' ' + Path['Topo'] + ' ' + Path['Sim'] + ' ' + Path['Out'] + ' ' + Path['Puppet'] + ' ' + PathRef
    if (sim.cluster['q'] and sim.cluster['par']):
      #jobName = 'LEquipe'
      runscrID = open(jobName,'w')
      os.system('chmod u+x ' + jobName)
      runscrID.write('#!/usr/bin/tcsh\n' + \
	      '#$ -V\n' + \
	      '#$ -N ' + jobName + '\n' + \
	      '#$ -e ' + Path['txt'] + ' -o ' + Path['txt'] + '\n' + \
	      '#$ -cwd\n' + \
	      '#$ -b y\n' + \
	      '#$ -pe mvapich2-nanna01 ' + str(nP) + '\n' + \
	      'mpiexec -machinefile $TMPDIR/machines' + ' -n $NSLOTS ' + cmd + '\n'\
	      'exit 0\n')
      runscrID.close()
      if hide:
	devnull = open('/dev/null', 'w')
	oldstdout_fno = os.dup(sys.stdout.fileno())
	os.dup2(devnull.fileno(), 1)
      os.system('qsub ./' + jobName)		# run SONet creation with suppressed output
      if hide:
	devnull.close()
	os.dup2(oldstdout_fno, 1)
	os.close(oldstdout_fno)
  ##### run flags
  ### scratch = cluster['scratch']
  ### par = cluster['par']
  #### scratch: 	save output data on scratch directory, otherwise on home
  #### par: 	run code in parallel, otherwise serially 
  #### q: 		submit the job to the queue, otherwise run on local machine (q=0,p=1, requires in homemade machinefiles indicating which processors should be used)
  ###bk = 12
  ###nP=bk*12	#if par=1, how many processors to use, nanna: 12 cores/node,12 nodes/block (fast connections only with block=> max 144 cores)
  
  ###nReal = 1	# what's that for?
  ####define run command calling correct version of compiled code
  
  ####runstr=' ../lequipe_compiled_par_nanna/LEquipe ' #should be compiled on the machine/cluster you want to run it on
  
  
  #####and associated cluster cmd
  
  ####if scratch==1:
    ####if (par==0 and q==0):
      ####clusterstr = ''
    ####elif (par==1 and q==0):
      ####clusterstr = 'mpirun -np ' + str(nP) + ' -machinefile ../lequipe_compiled_par_nanna/nannamachines'
    ####elif (par==0 and q==1):
      ####clusterstr = 'qsub -e ' + Path['txt'] + ' -o ' + Path['txt'] + ' -q %s.q -cwd -b y -V' % clustername
          #####'-v LD_LIBRARY_PATH=' + '"/usr/nld/netcdf-4.1.3/lib:/core/ge-6.2.1/lib/UNSUPPORTED-lx3.4.11-2.16-desktop-amd64:/usr/mpi/gcc/openmpi-1.4.3/lib64:/usr/nld/zeromq-4.0.4/lib64:/usr/nld/readline-6.2/lib:/usr/nld/intel-2013.0.028/mkl/lib/intel64:/usr/nld/python-2.7.6-freya/lib/python2.7/site-packages"'
    ####else: #par==1 && q==1
      ####print "that one's empty!" #current method generates shell script (See below). Would be nice to get this working out of matlab.
      #####return 0
  ####else:
    
    ####clusterstr = ''
  
  #run command
  
