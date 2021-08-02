execfile('Code/py_code/tubesize.py')
os.system('source /core/ge/NLD/common/settings.sh')

infoFile = sys.argv[1]
scratch = int(sys.argv[2])
net = network('r',infoFile)
#net.scriptName = 'SONET'

sim = simulation_data('r',infoFile,net,cluster={'scratch':scratch,'par':0,'q':1},mode='eft')
break_it = net.setup(sim,hide=1,call=1)

net.initDataSv(sim)

def read_it(time_now):
  for ps_ct in range(sim.steps['ps']):
    correct_ct = 0
    green_light = 0
    wait = 0
    while not green_light:
      if os.path.exists(net.Path['results_links'] + '_%d' % (ps_ct+1)) or os.path.exists(net.Path['results_links']):
	try:
	  ncid = netcdf.netcdf_file(net.Path['results_links'] + '_%d' % ps_ct,'r',mmap=False)
	  print "Now reading perturbation step %d" % ps_ct
	  green_light = 1
	except:
	  if not wait:
	    print 'Simulations are not yet ready (1)'
	  time.sleep(30)
	  wait = 1
      else:
	if not wait:
	  print 'Simulations are not yet ready (2)'
	time.sleep(30)
	wait = 1
    for co_ct in range(sim.steps['co']):
      
      for in_ct in range(sim.steps['in']):
	
	read = 0
	wait = 0
	
	link = ''.join(ncid.variables['links'][co_ct,in_ct,-1])
	ID_perturbed = ncid.variables['IDs'][co_ct,in_ct]
	
	while not read:
	  
	  try:
	    assert os.path.exists(link)
	    ncid_read = netcdf.netcdf_file(link,'r',mmap=False)
	    distances = ncid_read.variables['PTdistances'][:]
            net.DataSv['CDBTimes'][ps_ct][co_ct][in_ct] = ncid_read.variables['PTtimes'][:]
	    ncid_read.close()
	  
	    mask = np.invert(np.isnan(distances)) # mask entries, where no perturbation could be found
	    correct_ct += np.sum(mask)
	    net.DataSv['distances'][ps_ct][co_ct][in_ct][:] = distances
	    net.DataSv['div_track'][ps_ct][co_ct][in_ct][mask] = (distances >= 0.01)[mask].astype('int')
	    
	    read = 1
	    
	    # clean after reading
	    #os.remove(link)
	  except:
	    p = subprocess.Popen('qstat -u aschmidt | grep %s%s_P%d_%d | wc -l' % (net.scriptName[0],net.Hash[:2],ps_ct,ID_perturbed),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
	    out,err = p.communicate()
	    out = int(out)
	    if not out:
	      print 'Simulation data not there anymore - skip! (ID: %d, waited %d secs)' % (ID_perturbed,wait*60)
	      read = -1
	    else:
	      if not wait:
		print 'Simulation still in queue - waiting (ID: %d, out = %d)' % (ID_perturbed,out)
	      wait += 1
	      time.sleep(60)
    
    if correct_ct:
      net.DataSv['p'][ps_ct] = 1 - np.nansum(net.DataSv['div_track'][ps_ct])/correct_ct
    else:
      net.DataSv['p'][ps_ct] = np.nan
    
    net.DataSv['pert_size'][ps_ct] = ncid.variables['pert_size'].getValue()
    net.DataSv['correct_ct'][ps_ct] = correct_ct

    print 'probability of divergence at perturbation size ps = %5.3g is %5.3g (%d measurements)' % (net.DataSv['pert_size'][ps_ct],net.DataSv['p'][ps_ct],correct_ct)

    if os.path.exists(net.Path['results_links']) and not os.path.exists(net.Path['results_links'] + '_%d' % (ps_ct+1)):
      sim.steps['ps'] = ps_ct + 1
      break

    ncid.close()
  
  net.fluxSimSv(sim.steps)
  
  net.clean()	#clean up afterwards

### how to wait for the next perturbation steps to have finished?
print 'reading %d files from lists %s_xx...' % (sim.steps['ps']*sim.steps['total'],net.Path['results_links'])

time_now = time.strftime('%c')
read_it(time_now)
