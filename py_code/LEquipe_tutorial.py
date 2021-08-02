import os, imp, inspect, math
import matplotlib.pylab as plt
import matplotlib.colors as mcolors

# assert os.getcwd().split('/')[-1] == 'Program', 'Please run the program in the directory of LEquipe to properly set up the paths.'

directory = os.getcwd()

imp.load_source('read_write_code','py_code/read_write.py')
imp.load_source('file_management_code','py_code/file_management.py')
imp.load_source('support_code', 'py_code/support.py')
from read_write_code import *
from file_management_code import *
from support_code import *

############### define the parameters of the network here ###########
def LEquipe_tutorial(scratch=0,directory = directory, suppressMessages = 1):

    print('Which mode would you like to run?')
    mode_dict = {1:'train', 2:'LE', 3:'CLV', 4:'ISI', 5:'ISI2', 6:'mixed'}
    print(mode_dict.keys())
    for key,val in mode_dict.items():
        print('%d - LEquipe_%s' % (key,val))
    mode_idx = int(input('choose the number: '))
    mode = mode_dict[mode_idx]

    print('Running tutorial LEquipe_%s.' % mode)

    scriptname = inspect.getfile(inspect.currentframe()).split('/')[-1].split('.')[0]

    #PathDict = set_path(scriptname + '_' + mode,scratch)

    ParaNet = {}
    ParaSim = {}
    ParaTopo = {}

    ParaNet,ParaTopo,ParaSim = eval(mode)(ParaNet,ParaTopo,ParaSim)


    HashNet, PathDict['Net'], ParaNet = writeNet(PathDict['para'], ParaNet)
    HashTopo, PathDict['Topo'], ParaTopo = writeTopo(PathDict['para'], ParaNet, ParaTopo)
    HashSim, PathDict['Sim'], ParaSim = writeSim(PathDict['para'], ParaNet, ParaSim)

    HashDataOut = hashlib.sha1(np.array([HashNet, HashTopo, HashSim])).hexdigest()
    PathDict['Out'] = PathDict['results'] + 'DataOut-' + HashDataOut + '.nc'

    PathDict = sv_script_para(PathDict,ParaNet,ParaTopo,ParaSim,[PathDict['Net'],PathDict['Sim'],PathDict['Topo'],PathDict['Out']],['N','K','rateWnt','init'])

    # run the script
    run_script(PathDict,scratch,par=0,q=0)

    # read the output file and plot the results
    Data = readDataOut(PathDict['Out'])

    eval(mode + '_plot')(Data,ParaNet,ParaTopo,ParaSim)


######################################## code for tutorials starts here ##########################################

# -------------------------------------- for spike train --------------------------------------------
def train(ParaNet,ParaTopo,ParaSim):
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    ParaNet['N'] = 200		#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = -1		#coupling strength J0
    ParaSim['rateWnt'] = 5.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaNet['rapidness'] = 1	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    ## set the random graph with K synapses per neuron on average
    ParaTopo = random_graph(ParaNet['K'],ParaNet['N'])

    #set synapstic coupling strength (sqrt(K) scaling for the balanced state)
    ParaTopo['J'] = ParaNet['J0']/math.sqrt(ParaNet['K'])

    ##the external currents that yield the wanted firing rate can be well approximated by the balance equation
    ##f = -I0/(J0*tauM) then with the balanced state scaling we end up with Iext = sqrt(K)*I0
    ##synaptic time constant in case of cLIF or twoDlinear tauS
    ParaNet['Iext'] = -ParaNet['J0']*ParaSim['rateWnt']/1000*ParaNet['tauM']*np.sqrt(ParaNet['K'])

    ParaSim['SW'] = 100		#number of spikes per neuron during warmup
    ParaSim['TC'] = 1		#time duration of the calculation in seconds
    ParaSim['train'] = range(ParaNet['N'])	#neurons, whose spike times are saved

    return ParaNet, ParaTopo, ParaSim

def train_plot(Data,ParaNet,ParaTopo,ParaSim):
    fig1 = plt.figure()
    plt.plot(Data['trainTime'], Data['trainNeuron'], '.', markersize = 7)
    plt.xlabel('time (s)')
    plt.ylabel('neurons')
    plt.xlim([0,Data['trainTime'].max()])
    plt.show()


# -------------------------------------- for CLV-spectra --------------------------------------------

def CLV(ParaNet,ParaTopo,ParaSim):
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    ParaNet['N'] = 100		#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = -1		#coupling strength J0
    ParaSim['rateWnt'] = 5.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaNet['rapidness'] = 10	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    ## set the random graph with K synapses per neuron on average
    ParaTopo['post'], ParaTopo['row_length'] = random_graph(ParaNet['K'],ParaNet['N'])

    #set synapstic coupling strength (sqrt(K) scaling for the balanced state)
    ParaTopo['J'] = ParaNet['J0']/math.sqrt(ParaNet['K'])

    ##the external currents that yield the wanted firing rate can be well approximated by the balance equation
    ##f = -I0/(J0*tauM) then with the balanced state scaling we end up with Iext = sqrt(K)*I0
    ##synaptic time constant in case of cLIF or twoDlinear tauS
    ParaNet['Iext'] = -ParaNet['J0']*ParaSim['rateWnt']/1000*ParaNet['tauM']*math.sqrt(ParaNet['K'])

    ParaSim['SW'] = 100		#number of spikes per neuron during warmup
    ParaSim['train'] = range(ParaNet['N'])	#neurons, whose spike times are saved

    #Lyapunov exponent parameters
    if ParaNet['NeuronType'] < 10:
        ParaSim['LyapunovExp'] = ParaNet['N']
    else:
        ParaSim['LyapunovExp'] = 2*ParaNet['N']

        ParaSim['SC'] = 15		#avg. number of spikes per neuron in the calculation
        ParaSim['SWONS'] = 10		#warmup of the ONSE
        ParaSim['ONstep'] = 1		#orthonormalization step size

    #covariant Lyapunov vectors
    ParaSim['CLV'] = 1		#calculate the CLVs
    ParaSim['SWCLV'] = 10		#warmup of the CLVs (must be < ParaSim['SC'])

    return ParaNet, ParaTopo, ParaSim

def CLV_plot(Data,ParaNet, ParaTopo, ParaSim):

    fig = plt.figure()

    ax1 = plt.subplot(121)
    axcb = plt.subplot(1,30,16)
    ax2 = plt.subplot(233)
    ax3 = plt.subplot(236)

    levs = range(80)
    assert len(levs) % 2 == 0, 'N levels must be even.'
    zero_val = abs(Data['localLE'].min())/(Data['localLE'].max()-Data['localLE'].min())
    rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors=[(0,(0,0,1)),(zero_val,(1,1.,1)),(1,(1,0,0))],N=len(levs)-1,)
    dim0 = Data['localLE'].shape[0]
    dim1 = Data['localLE'].shape[1]
    times = np.tile(Data['LEtimes'][0:dim0,np.newaxis], (1,dim1))	#times = np.tile(Data['LEtimes'][0:dim0], (dim1,1))
    index = np.tile(np.linspace(1/float(dim1),1,dim1), (dim0,1))	#index = np.tile(np.linspace(1/float(dim1),1,dim1)[:,np.newaxis], (1,dim0))

    pc = ax1.pcolor(times, index, Data['localLE'],cmap=rwb)
    ax1.set_xlabel('time (ms)')
    ax1.set_ylabel('i / N')
    ax1.set_xlim([0,np.max(Data['LEtimes'][0:dim0])])
    ax1.set_title('local Lyapunov exponents per spike')
    plt.colorbar(pc, cax=axcb)

    ax2.plot(Data['trainTime'], Data['trainNeuron'], '.', markersize = 5)
    ax2.set_xlabel('time')
    ax2.set_ylabel('neurons')
    ax2.set_xlim([0, max(Data['trainTime'])])
    ax2.set_title('spike train')

    ax3.plot(np.linspace(1/float(ParaSim['LyapunovExp'][0]),1,ParaSim['LyapunovExp'][0]), Data['LyapunovExponents'])
    ax3.plot(np.linspace(1/float(ParaSim['LyapunovExp'][0]),1,ParaSim['LyapunovExp'][0]), Data['LEclv'])
    ax3.set_ylabel('lambda_i ( s ^{ -1})')
    ax3.set_xlabel('i / N')
    ax3.set_title('Lyapunov spectra')
    ax3.set_ylim([np.min(Data['LyapunovExponents']),np.max(Data['LyapunovExponents'])])

    plt.show()


# -------------------------------------- for Lyapunov Exponents --------------------------------------------

def LE(ParaNet, ParaTopo, ParaSim):
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    ParaNet['N'] = 200		#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = -1		#coupling strength J0
    ParaSim['rateWnt'] = 5.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaNet['rapidness'] = 1	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    ## set the random graph with K synapses per neuron on average
    ParaTopo['post'], ParaTopo['row_length'] = random_graph(ParaNet['K'],ParaNet['N'])

    #set synapstic coupling strength (sqrt(K) scaling for the balanced state)
    ParaTopo['J'] = ParaNet['J0']/math.sqrt(ParaNet['K'])

    ##the external currents that yield the wanted firing rate can be well approximated by the balance equation
    ##f = -I0/(J0*tauM) then with the balanced state scaling we end up with Iext = sqrt(K)*I0
    ##synaptic time constant in case of cLIF or twoDlinear tauS
    ParaNet['Iext'] = -ParaNet['J0']*ParaSim['rateWnt']/1000*ParaNet['tauM']*math.sqrt(ParaNet['K'])

    ParaSim['SW'] = 100		#number of spikes per neuron during warmup
    ParaSim['train'] = range(ParaNet['N'])	#neurons, whose spike times are saved

    ParaSim['LyapunovExp'] = ParaNet['N']	# number of Lyapunov exponents
    ParaSim['SC'] = 10		# avg. number of spikes per neuron in the calculation

    return ParaNet, ParaTopo, ParaSim

def LE_plot(Data, ParaNet, ParaTopo, ParaSim):
    fig = plt.figure()

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)

    ax1.plot(Data['trainTime'],Data['trainNeuron'],'.',markersize = 5)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('neurons')
    ax1.set_xlim([0,np.max(Data['trainTime'])])

    ax2.plot(np.linspace(1/ParaSim['LyapunovExp'][0],1,ParaSim['LyapunovExp'][0]),Data['LyapunovExponents'])
    ax2.set_xlabel('i/N')
    ax2.set_ylabel('\lambda_i (s^{-1}')

    plt.show()


# -------------------------------------- for inter spike intervalls --------------------------------------------

def ISI(ParaNet, ParaTopo, ParaSim):
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    ParaNet['N'] = 200		#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = -1		#coupling strength J0
    ParaSim['rateWnt'] = 5.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaNet['rapidness'] = 1	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    ## set the random graph with K synapses per neuron on average
    ParaTopo['post'], ParaTopo['row_length'] = random_graph(ParaNet['K'],ParaNet['N'])

    #set synapstic coupling strength (sqrt(K) scaling for the balanced state)
    ParaTopo['J'] = ParaNet['J0']/math.sqrt(ParaNet['K'])

    ##the external currents that yield the wanted firing rate can be well approximated by the balance equation
    ##f = -I0/(J0*tauM) then with the balanced state scaling we end up with Iext = sqrt(K)*I0
    ##synaptic time constant in case of cLIF or twoDlinear tauS
    ParaNet['Iext'] = -ParaNet['J0']*ParaSim['rateWnt']/1000*ParaNet['tauM']*math.sqrt(ParaNet['K'])

    ParaSim['SW'] = 100		#number of spikes per neuron during warmup
    ParaSim['SC'] = 10		#average number of spikes per neuron

    ParaSim['ISIneurons'] = range(ParaNet['N'])	#neurons whose spike statistics are calculated
    ParaSim['ISIstats'] = 1			#kinda moments that are being calculated

    return ParaNet, ParaTopo, ParaSim

def ISI_plot(Data, ParaNet, ParaTopo, ParaSim):
    fig = plt.figure()
    plt.plot(Data['rateNeurons'], '.', markersize=7)
    plt.xlabel('neuron index')
    plt.ylabel('firing rate (Hz)')
    plt.show()


def ISI2(ParaNet, ParaTopo, ParaSim):
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    ParaNet['N'] = 1000		#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = -1		#coupling strength J0
    ParaSim['rateWnt'] = 5.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaNet['rapidness'] = 1	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    ## set the random graph with K synapses per neuron on average
    ParaTopo['post'], ParaTopo['row_length'] = random_graph(ParaNet['K'],ParaNet['N'])

    #set synapstic coupling strength (sqrt(K) scaling for the balanced state)
    ParaTopo['J'] = ParaNet['J0']/math.sqrt(ParaNet['K'])

    ##the external currents that yield the wanted firing rate can be well approximated by the balance equation
    ##f = -I0/(J0*tauM) then with the balanced state scaling we end up with Iext = sqrt(K)*I0
    ##synaptic time constant in case of cLIF or twoDlinear tauS
    ParaNet['Iext'] = -ParaNet['J0']*ParaSim['rateWnt']/1000*ParaNet['tauM']*math.sqrt(ParaNet['K'])

    ParaSim['SW'] = 100		#number of spikes per neuron during warmup
    ParaSim['SC'] = 1000		#average number of spikes per neuron

    ParaSim['ISIneurons'] = range(ParaNet['N'])	#neurons whose spike statistics are calculated
    ParaSim['ISIstats'] = 4			#kinda moments that are being calculated
    ParaSim['ISIbins'] = 20

    return ParaNet, ParaTopo, ParaSim

def ISI2_plot(Data, ParaNet, ParaTopo, ParaSim):
    fig = plt.figure()

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    ax1.plot(Data['rateDistX'],Data['rateDistY'],'.-')
    ax1.set_xlabel('rates')

    ax2.plot(Data['cvDistX'],Data['cvDistY'],'.-')
    ax2.set_xlabel('cv')

    ax3.plot(Data['skewnessDistX'],Data['skewnessDistY'],'.-')
    ax3.set_xlabel('skewness')

    ax4.plot(Data['kurtosisDistX'],Data['kurtosisDistY'],'.-')
    ax4.set_xlabel('kurtosis')

    plt.show()


def mixed(ParaNet, ParaTopo, ParaSim):

    #print "This seems to be bugged still - rather try something else"
    ParaNet['NeuronType'] = 1	#neuron model (viz. writeNet.py), AP onset rapidness in case of rapid theta neurons

    NI = 100
    NE = 400

    ParaNet['N'] = NI + NE	#number of neurons N
    ParaNet['K'] = 50 		#number of synapses per neuron K
    ParaNet['J0'] = 1		#coupling strength J0
    ParaSim['rateWnt'] = 10.	#network-averaged firing rate in Hz
    ParaNet['tauM'] = 10.		#membrane time constant

    ParaSim['train'] = range(ParaNet['N'])
    ParaSim['LyapunovExp'] = ParaNet['N']
    #ParaNet['rapidness'] = 10	#AP onset rapidness in case of rapid theta neurons

    #tauS: synaptic time constant in case of cLIF or twoDlinear
    ParaNet['twoDlinear'] = {'alpha':1., 'beta':0., 'gamma':0., 'delta':1., 'Cw':0., 'tauS': ParaNet['tauM']/2.}

    A_EE = 0.01*ParaNet['J0']/math.sqrt(ParaNet['K'])*(np.random.uniform(0,1,(NE,NE)) < ParaNet['K']/float(NE))
    A_II = -ParaNet['J0']/math.sqrt(ParaNet['K'])*(np.random.uniform(0,1,(NI,NI)) < ParaNet['K']/float(NI))

    A_IE = 0.01*ParaNet['J0']/math.sqrt(ParaNet['K'])*(np.random.uniform(0,1,(NE,NI)) < ParaNet['K']/float(NI))
    A_EI = -ParaNet['J0']/math.sqrt(ParaNet['K'])*(np.random.uniform(0,1,(NI,NE)) < ParaNet['K']/float(NE))

    A = np.concatenate((np.concatenate((A_EE,A_EI),axis=0),np.concatenate((A_IE,A_II),axis=0)),axis=1)

    ParaTopo['post'],ParaTopo['J'],ParaTopo['row_length'] = adjacencymatrix2postsynapticvector(A)

    return ParaNet, ParaTopo, ParaSim

def mixed_plot(Data, ParaNet, ParaTopo, ParaSim):
    fig = plt.figure()

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)

    ax1.plot(Data['trainTime'],Data['trainNeuron'],'.',markersize = 5)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('neurons')
    ax1.set_xlim([0,ParaSim['TC'][0]])

    ax2.plot(np.linspace(1/ParaSim['LyapunovExp'][1],1,ParaSim['LyapunovExp'][1]),Data['LyapunovExponents'])
    ax2.set_xlabel('i/N')
    ax2.set_ylabel('\lambda_i (s^{-1}')

    plt.show()


LEquipe_tutorial()
