import os, time, subprocess
import numpy as np
import matplotlib.pyplot as plt
#findPIDtxt = 'findPID.txt'
#writeRam = 'timeRam.txt'

#write = open(writeRam)
running = 1
running_ct = 0
storage = []
while running:
  handle = subprocess.Popen('ps aux | grep LEquipe',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
  out,err = handle.communicate()
  read = out.split('\n')
  PID = []
  
  for line in read:
    if ('./Code/lequipe/LEquipe' in line) and not ('sh -c' in line):
      
      for entry in line.split(' '):
	if not (entry == ''):
	  try:
	    PID.append(int(entry))
	    break
	  except:
	    1
      
      for t in range(10):
	command = "top -b -n 1 -d 1 -p %d | sed -n '8p'" % PID[-1]
	
	handle = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE)
	out,err = handle.communicate()
	if not (out == '\n'):
	  storage.append([])
	  idx = 0
	  for entry in out.split(' '):
	    if not (entry in ['','\n']):
	      storage[-1].append(entry)
  
  if not len(PID):
    running_ct += 1
    time.sleep(0.5)
  else:
    running_ct = 0
  #time.sleep(0.5)
        
  if running_ct >= 100:
    running = 0
  
  #print running_ct, running
  
#print storage

good = 1

Data = {}

for key in ['PID','time','CPU','mem','virt','res','shr']:
  Data[key] = []

PIDtmp = 0

for i in range(len(storage)):
  if not (PIDtmp == int(storage[i][0])):
    for key in ['PID','time','CPU','mem','virt','res','shr']:
      Data[key].append([])
  
  PIDtmp = int(storage[i][0])
  Data['PID'][-1].append(PIDtmp)
  
  timeTmp = storage[i][10].split(':')
  timeTmp2 = timeTmp[1].split('.')
  time = 60*float(timeTmp[0]) + float(timeTmp2[0]) + float(timeTmp2[1])/100.
  Data['time'][-1].append(time)
  Data['CPU'][-1].append(float(storage[i][8]))
  Data['mem'][-1].append(float(storage[i][9]))
  
  if storage[i][4][-1] == 'm':
    virtTmp = float(storage[i][4][:-1])
  elif storage[i][4][-1] == 'g':
    virtTmp = float(storage[i][4][:-1])*10**3
  else:
    virtTmp = float(storage[i][4])*10**(-3)
  Data['virt'][-1].append(virtTmp)

  if storage[i][5][-1] == 'm':
    RESTmp = float(storage[i][5][:-1])
  elif storage[i][5][-1] == 'g':
    RESTmp = float(storage[i][5][:-1])*10**3
  else:
    RESTmp = float(storage[i][5])*10**(-3)
  Data['res'][-1].append(RESTmp)

  if storage[i][6][-1] == 'm':
    SHRTmp = float(storage[i][6][:-1])
  elif storage[i][6][-1] == 'g':
    SHRTmp = float(storage[i][6][:-1])*10**3
  else:
    SHRTmp = float(storage[i][6])*10**(-3)
  Data['shr'][-1].append(SHRTmp)

plt.figure()
for i in range(len(Data['time'])):
  plt.plot(Data['time'][i],Data['mem'][i])
plt.ylim([0,50])
plt.show(block=False)

plt.figure()
for i in range(len(Data['time'])):
  plt.plot(Data['time'][i],Data['virt'][i],'k',label='VIRT')
  plt.plot(Data['time'][i],Data['res'][i],'r',label='RES')
  plt.plot(Data['time'][i],Data['shr'][i],'b',label='SHR')
plt.ylabel('MB')
plt.xlabel('time')
plt.legend()
plt.show(block=False)
