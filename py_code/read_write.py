import numpy as np
import random as rnd

import hashlib, math, os, imp

from scipy.io import netcdf

imp.load_source('support_code', './py_code/support.py')
from support_code import *

class read_write:
  ################# write stuff for netcdf ##################

    def writeNet(self, suppressMessages=1):

        ## Neuron Type:
        # 1 = rapid theta
        # 2 = LIF
        # 3 = type1type2
        # 10 = twoDlinear (RIF->cLif)
        # 11 = cLIF (gamma = 1/3)
        # 12 = cLIF (gamma = 1/2)
        # 13 = cLIF (gamma = 2)
        # 14 = cLIF (gamma = 3)
        # (5-9 = other phase neurons)
        # 5 = poisson neuron
        # 51 = inhomogenous poisson neuron. Note how variables are used:
        # Lambda(t) = poissonrate + Iext * sin(2*pi*tauM*TC + rapidness)
        # so poissonrate := constant poissonrate
        # Iext := scaling of sine
        # tauM := frequency of oscillation in Hz,
        # so rapidness := phase shift of sine.
        # 6 = gamma neuron
        # 7 = puppet neuron -> fed with precalculated spike train
        # 15 = cQIF
        # 16 = rcQIF

        # check for network homogeneity and consistent datasizes (1 or N) and -types (np.array)
        set_default(self.ParaNet,'N',1000)
        self.ParaNet, HomogNetwork = network_check(self.ParaNet,self.Const['N'])
        ## set the default values if not defined in Para
        set_default(self.ParaNet,'K',100)						#number of synapses per neuron K
        set_default(self.ParaNet,'J0',-1)						#coupling strength J0
        set_default(self.ParaNet,'NeuronType',2)					#neuron model (see top)
        set_default(self.ParaNet,'tauM',10.)
        set_default(self.ParaNet,'Iext',1)
        set_default(self.ParaNet,'seedInit',1)
        set_default(self.ParaNet,'rapidness',1.)					#AP onset rapidness in case of rapid theta neurons
        set_default(self.ParaNet,'poissonrate',100)					#spikes/second from poisson neurons
        set_default(self.ParaNet,'type1type2Para',0)
        if not any([(item == 'twoDlinear') for item in self.ParaNet.keys()]):
            self.ParaNet['twoDlinear'] = {}
            set_default(self.ParaNet['twoDlinear'],'alpha',1)
            set_default(self.ParaNet['twoDlinear'],'beta',0)
            set_default(self.ParaNet['twoDlinear'],'gamma',0)
            set_default(self.ParaNet['twoDlinear'],'delta',1)
            set_default(self.ParaNet['twoDlinear'],'Cw',0)
            set_default(self.ParaNet['twoDlinear'],'tauS',self.ParaNet['tauM'][0]/2.)

        # find neuron Type and set according initial values
        if not ('init' in self.ParaNet.keys()):
            self.ParaNet['init'] = np.zeros(self.ParaNet['N'][0])
            #print self.ParaNet['N']
            for n in range(int(self.ParaNet['N'][0])):
                if len(self.ParaNet['NeuronType']) > 1:
                    initFlag = self.ParaNet['NeuronType'][n]
                else:
                    initFlag = self.ParaNet['NeuronType'][0]

        if initFlag in [1,3]:
            self.ParaNet['init'][n] = 2*math.pi*(rnd.random() - 0.5)
        elif initFlag == 2:
            self.ParaNet['init'][n] = rnd.random()
        elif initFlag in [10,11,12,13,14,15]:
            self.ParaNet['init'][n] = [rnd.random(), 0]	#in matlab: vector
        elif initFlag in [5,51]:
            self.ParaNet['init'][n] = rnd.randint(1,2**32-1)
        elif initFlag == 7:
            assert 'puppetTrain' in self.ParaNet.keys(), 'Puppet neurons are specified in the network though no spike train is.'
        if not suppressMessages:
            print('init is set to default value according to neuron model')

        #Since Lequipe expects vectors for parameters whenever network is inhomogeneous in any way,
        #we need to convert any scalar, i.e. homogeneous, fields to unfortunately redundant vectors.
        #We should obviate the need for this by writing the check into the cpp code.
        if not HomogNetwork:
            for item in self.ParaNet.items():
                if (not (item[0] in ['N','seedInit','init','J0']) and (len(item[1])==1 or type(item[1])==dict)):
                    if not (type(item[1])==dict):
                        self.ParaNet[item[0]] = item[1]*np.ones(self.ParaNet['N'][0])
                    else:
                        for sub_item in item[1].items():
                            if (len(sub_item[1]) == 1):
                                self.ParaNet[item[0]][sub_item[0]] = sub_item[1]*np.ones(self.ParaNet['N'][0])

        ## Generate Hash and filename
        Hash = hashing(self.ParaNet,['Iext','N','NeuronType','rapidness','tauM','twoDlinear','init'])
        self.Path['Net'] = self.Path['inpara'] + 'ParaNetw-' + Hash + '.nc'
        # ParaNeurons goes in front to avoid ncdump error (bug)
        # ncdump: name begins with space or control-character: 1

        ## Open netCDF files.
        if not os.path.exists(self.Path['Net']):

            if not suppressMessages:
                print('writing neuron netcdf file: %s' % self.Path['Net'])
            ncid = netcdf.netcdf_file(self.Path['Net'], 'w')

            ## Define the dimensions of the variables.
            ncid.createDimension('one', 1)

            # Define new variables in the neuron file.
            VarNeuron_N = ncid.createVariable('N', 'i',('one',))

            if not HomogNetwork:
            	ncid.createDimension('N', self.ParaNet['N'][0])
            	VarNeuron_NeuronType = ncid.createVariable('NeuronType','i', ('N',))
            	VarNeuron_tauM = ncid.createVariable('tauM', 'd', ('N',))
            	VarNeuron_Iext = ncid.createVariable('Iext', 'd', ('N',))
            	VarNeuron_rap = ncid.createVariable('rapidness', 'd', ('N',))
            	VarNeuron_poissonrate = ncid.createVariable('poissonrate', 'd', ('N',))
            	VarNeuron_type1type2 = ncid.createVariable('type1type2Para', 'd', ('N',))
            	VarNeuron_a = ncid.createVariable('twoDlinear_alpha', 'd', ('N',))
            	VarNeuron_b = ncid.createVariable('twoDlinear_beta', 'd', ('N',))
            	VarNeuron_c = ncid.createVariable('twoDlinear_gamma', 'd', ('N',))
            	VarNeuron_d = ncid.createVariable('twoDlinear_delta', 'd', ('N',))
            	VarNeuron_Cw = ncid.createVariable('twoDlinear_Cw', 'd', ('N',))
            	VarNeuron_tauS = ncid.createVariable('twoDlinear_tauS', 'd', ('N',))
            else:
            	VarNeuron_NeuronType = ncid.createVariable('NeuronType','i',('one',))
            	VarNeuron_tauM = ncid.createVariable('tauM', 'd', ('one',))
            	VarNeuron_Iext = ncid.createVariable('Iext', 'd', ('one',))
            	VarNeuron_rap = ncid.createVariable('rapidness', 'd', ('one',))
            	VarNeuron_poissonrate = ncid.createVariable('poissonrate', 'd', ('one',))
            	VarNeuron_type1type2 = ncid.createVariable('type1type2Para', 'd', ('one',))
            	VarNeuron_a = ncid.createVariable('twoDlinear_alpha', 'd', ('one',))
            	VarNeuron_b = ncid.createVariable('twoDlinear_beta', 'd', ('one',))
            	VarNeuron_c = ncid.createVariable('twoDlinear_gamma', 'd', ('one',))
            	VarNeuron_d = ncid.createVariable('twoDlinear_delta', 'd', ('one',))
            	VarNeuron_Cw = ncid.createVariable('twoDlinear_Cw', 'd', ('one',))
            	VarNeuron_tauS = ncid.createVariable('twoDlinear_tauS', 'd', ('one',))

            if self.topo == 'p':
                ncid.createDimension('Ndrive',self.ParaNet['N'][0] - self.Const['N'])

            state = []
            stateIdx = np.zeros(len(self.ParaNet['init']))
            idx = 0
            for i in range(len(self.ParaNet['init'])):
                stateIdx[i] = idx
                try:
                    if self.ParaNet['init'][i][0] == self.ParaNet['init'][i][1]:
                        state.append(self.ParaNet['init'][i][0])
                        idx += 1
                    else:
                        state.extend(self.ParaNet['init'][i]) ##ok<AGROW>
                        idx += len(self.ParaNet['init'][i])
                except:
                    state.append(self.ParaNet['init'][i]) ##ok<AGROW>
                    idx += 1

                    ncid.createDimension('initStatesSz', len(state))
                    ncid.createDimension('initStatesIdxSz', len(stateIdx))

                    VarNeuron_state = ncid.createVariable('initStates', 'd', ('initStatesSz',))
                    VarNeuron_stateIdx = ncid.createVariable('initStatesIdx', 'i', ('initStatesIdxSz',))

                    ## Write data to variable.
                    VarNeuron_N[:] = self.ParaNet['N']
                    VarNeuron_NeuronType[:] = self.ParaNet['NeuronType']
                    VarNeuron_tauM[:] = self.ParaNet['tauM']
                    VarNeuron_Iext[:] = self.ParaNet['Iext']
                    VarNeuron_state[:] = state
                    VarNeuron_stateIdx[:] = stateIdx
                    VarNeuron_rap[:] = self.ParaNet['rapidness']
                    VarNeuron_poissonrate[:] = self.ParaNet['poissonrate']
                    VarNeuron_type1type2[:] = self.ParaNet['type1type2Para']
                    VarNeuron_a[:] = self.ParaNet['twoDlinear']['alpha']
                    VarNeuron_b[:] = self.ParaNet['twoDlinear']['beta']
                    VarNeuron_c[:] = self.ParaNet['twoDlinear']['gamma']
                    VarNeuron_d[:] = self.ParaNet['twoDlinear']['delta']
                    VarNeuron_Cw[:] = self.ParaNet['twoDlinear']['Cw']
                    VarNeuron_tauS[:] = self.ParaNet['twoDlinear']['tauS']

                    ## close file
                    ncid.close()


    def writePuppet(self,suppressMessages=1):

        if not 'puppetTrain' in self.ParaNet.keys():
            self.Path['Puppet'] = 'none'
            return

        Hash = hashing(self.ParaNet,['puppetTrain'])
        self.Path['Puppet'] = self.Path['results'] + 'Refer/P-' + Hash + '.nc'

        ## Open netCDF files.
        if not os.path.exists(self.Path['Puppet']):
            if not suppressMessages:
                print('writing puppet netcdf file: %s' % self.Path['Puppet'])
            ncid = netcdf.netcdf_file(self.Path['Puppet'], 'w')

            ncid.createDimension('Ndrive',self.Const['N'])
            ncid.createDimension('puppetTrainSz',self.puppetTrainSz)
            VarNeuron_puppetTrain = ncid.createVariable('puppetTrain','d',('Ndrive','puppetTrainSz'))

            VarNeuron_puppetTrain[:] = self.ParaNet['puppetTrain']

            ncid.close()


    def writeTopo(self, suppressMessages=1):

        self.ParaTopo, trash = network_check(self.ParaTopo,self.Const['N'])
        self.ParaNet, trash = network_check(self.ParaNet,self.Const['N'])	#need to convert stuff to lists

        ## Find out whether synapses are homogenous or not:
        HomogSynapses = 1

        if len(self.ParaNet['J0']) > 1:
            HomogSynapses = 0
        if ('J' in self.ParaTopo.keys()):
            if (len(self.ParaTopo['J']) > 1):
                HomogSynapses = 0
        if ('pSyn' in self.ParaTopo.keys()):
            if (len(self.ParaTopo['pSyn']) > 1):
                HomogSynapses = 0

        if not ('pSyn' in self.ParaTopo.keys()):
            self.ParaTopo['pSyn'] = [1]		#0 = false, 1 = true
            if not HomogSynapses:
                self.ParaTopo['pSyn'] = self.ParaTopo['pSyn']*len(self.ParaTopo['postSyn'])
                if not suppressMessages:
                    print('pSyn is set to default value 1 with HomogSynapses = %d' % HomogSynapses)

        if not ('J' in self.ParaTopo.keys()):	# sanity check
            if HomogSynapses:
                self.ParaTopo['J'] = [self.ParaNet['J0'][0]/math.sqrt(self.Const['K'])]	# set synaptic coupling strength (sqrt(K) scaling for the balanced state)
            else:
                self.ParaTopo['J'] = [] #np.zeros(len(Para['postSyn']))
        for i in range(self.Const['N']):
            self.ParaTopo['J'].extend([self.ParaNet['J0'][i]/math.sqrt(self.Const['K'])]*self.ParaTopo['outDegree'][i])	# set post synapses for each neuron
            if not suppressMessages:
                print('J is set to default value J0/sqrt(K)')

        ## Generation of the Hash and the filename
        if len(self.ParaTopo['outDegree']) < 10000: # for small networks, take Hash of all Para.*
            Hash = hashing(self.ParaTopo,['outDegree','inDegree','postSyn','preSyn','J','pSyn'])

        else: # for networks larger 1000 neurons, only every nth entry of post, J and pSyn are included in the Hash.
            idx = np.linspace(0,len(self.ParaTopo['postSyn'])-1,1000).astype(int)
            self.ParaTopo['post_hash'] = list(np.array(self.ParaTopo['postSyn'])[idx])
            if not HomogSynapses:
            	self.ParaTopo['J_hash'] = list(np.array(self.ParaTopo['J'])[idx])
            	self.ParaTopo['pSyn_hash'] = list(np.array(self.ParaTopo['pSyn'])[idx])
            	Hash = hashing(self.ParaTopo,['outDegree','post_hash','J_hash','pSyn_hash'])
            else:
            	Hash = hashing(self.ParaTopo,['outDegree','post_hash','J','pSyn'])
        self.Path['Topo'] = self.Path['inpara'] + 'ParaTopo-' + Hash + '.nc'


        ## Open netCDF files.
        if not os.path.exists(self.Path['Topo']):
            if not suppressMessages:
                print('writing topology netcdf file: %s' % self.Path['Topo'])
            ncid = netcdf.netcdf_file(self.Path['Topo'], 'w') #'64BIT_OFFSET');

            ## Define the dimensions of the variables.
            ncid.createDimension('one', 1)
            ncid.createDimension('N', len(self.ParaTopo['outDegree']))
            ncid.createDimension('synapses', len(self.ParaTopo['postSyn']))

            ## Define new variables in the topology file.
            VarTopology_post = ncid.createVariable('postSyn', 'i',('synapses',))
            VarTopology_pre = ncid.createVariable('preSyn', 'i',('synapses',))
            VarTopology_outDeg = ncid.createVariable('outDegree', 'i', ('N',))
            VarTopology_inDeg = ncid.createVariable('inDegree', 'i', ('N',))

            if HomogSynapses:
            	VarTopology_J = ncid.createVariable('J', 'd', ('one',))
            	VarTopology_pSyn = ncid.createVariable('pSyn', 'd', ('one',))
            else:
            	VarTopology_J = ncid.createVariable('J', 'd', ('synapses',))
            	VarTopology_pSyn = ncid.createVariable('pSyn', 'd', ('synapses',))

            #NetCDFVariable = a netCDFVariable object
            #attName - the name of the attribute
            #attValue - the value of the attribute
            setattr(VarTopology_post, 'outDegree', 'outDegree')
            setattr(VarTopology_pre, 'inDegree', 'inDegree')
            if not HomogSynapses:
            	setattr(VarTopology_J, 'outDegree', 'outDegree')
            	setattr(VarTopology_pSyn, 'outDegree', 'outDegree')


            ## Write data to variable.
            VarTopology_post[:] = self.ParaTopo['postSyn']
            VarTopology_pre[:] = self.ParaTopo['preSyn']
            VarTopology_J[:] = self.ParaTopo['J']
            VarTopology_pSyn[:] = self.ParaTopo['pSyn']
            VarTopology_outDeg[:] = self.ParaTopo['outDegree']
            VarTopology_inDeg[:] = self.ParaTopo['inDegree']

            if 'phat' in self.ParaTopo.keys(): # write actual network statistics
            	VarTopology_p_hat = ncid.createVariable('p_hat', 'd', ('one',))
            	VarTopology_recip_hat = ncid.createVariable('alpha_recip_hat', 'd', ('one',))
            	VarTopology_conv_hat = ncid.createVariable('alpha_conv_hat', 'd', ('one',))
            	VarTopology_div_hat = ncid.createVariable('alpha_div_hat', 'd', ('one',))
            	VarTopology_chain_hat = ncid.createVariable('alpha_chain_hat', 'd', ('one',))

            	VarTopology_p_hat[:] = self.ParaTopo['phat']
            	VarTopology_recip_hat[:] = self.ParaTopo['alphahat_recip']
            	VarTopology_conv_hat[:] = self.ParaTopo['alphahat_conv']
            	VarTopology_div_hat[:] = self.ParaTopo['alphahat_div']
            	VarTopology_chain_hat[:] = self.ParaTopo['alphahat_chain']

            ## close file
            ncid.close()


    def writeSim(self, suppressMessages=1):

        # check for consistent lengths of input
        self.ParaSim, HomogNetwork = network_check(self.ParaSim,self.Const['N'])

        ## Set default simulation parameters

        ## time / spikes per neuron for rate finding
        set_default(self.ParaSim,'TR',0)
        set_default(self.ParaSim,'SR',0)
        ## time / spikes per neuron for simulation warmup
        set_default(self.ParaSim,'TW',0)		# minimum simulation time in s
        set_default(self.ParaSim,'SW',0)
        ## time / spikes per neuron for calculations
        set_default(self.ParaSim,'TC',0)		# provide either SC or TC and either SW or TW
        set_default(self.ParaSim,'SC',0)

        set_default(self.ParaSim,'pd',0)	# gives the number of perturbation directions, 0 = no perturbing in LEquipe

        ## rate data
        set_default(self.ParaSim,'rateWnt',0)	# set to 0, if the external currents just set the final rate
        set_default(self.ParaSim,'rateWntSubN',0)	# set to 0, if the external currents just set the final rate
        set_default(self.ParaSim,'pR',0.01)		# precision of the desired rate

        ## for calculation of Lyapunov Exponents
        set_default(self.ParaSim,'LyapunovExp',0)	# number of LEs to be calculated (starting from highest)
        set_default(self.ParaSim,'seedONS',1)	# seed for initial setup of ONS
        set_default(self.ParaSim,'SWONS',1)		# spikes per neuron during warmup of ONS for LEs
        set_default(self.ParaSim,'ONstep',1)		# stepsize for orthonormalization during LE calculation
        set_default(self.ParaSim,'pLE',0)			# precision of LE convergence
        set_default(self.ParaSim,'LyapunovExpConvergence',0)		# save convergence of LE (file might be large if on)
        set_default(self.ParaSim,'CLV',0)				# calculate covariant LE vectors
        set_default(self.ParaSim,'SWCLV',0)		# spikes per neuron during warmup of ONS for covariant LE vectors
        set_default(self.ParaSim,'randomizedLEspectra',0)		# don't really know... ???

        ## when calculating LEs only in a subnetwork
        set_default(self.ParaSim,'subN',0)				# size of subnetwork
        set_default(self.ParaSim,'subLyapunovExp',max(0,self.ParaSim['subN'][0]))		# number of LEs to be calculated in subN
        set_default(self.ParaSim,'subLyapunovExpConvergence',0)	# save convergence of subN - LE (file might be large if on)

        ## other data to be kept for analysis
        set_default(self.ParaSim,'saveFinalState',1)		# save final state of simulation
        set_default(self.ParaSim,'ISIneurons',-1)		# indices of neurons, whose ISI will be calculated (-1 = none)
        set_default(self.ParaSim,'ISIstats',0)		# moment of the ISI to be calculated (1=mean, 2=coeff. of variation, 3=skewness, 4=kurtosis)
        set_default(self.ParaSim,'ISIbins',0)		# number of bins for firing rate distribution
        set_default(self.ParaSim,'train',-1)			# indices of neurons whos spike train is recorded (-1 = none)
        set_default(self.ParaSim,'synchrony',0)		# measures synchrony chi
        set_default(self.ParaSim,'ISIdecorr',-1)		# as in train - neurons for which ISI decorr statistics are calculated

        ## adding external currents
        set_default(self.ParaSim,'addCur',0)
        # sanity checks for external current
        if self.ParaSim['addCur'][0]:
            set_default(self.ParaSim,'addCurHomo',1)
            CurVal = np.array(['addCurTime','addCurIext','addCurNeurons'])
            if self.ParaSim['addCurHomo'][0]:
                assert (sum([(item == CurVal[0:2]).any() for item in self.ParaSim.keys()]) == 2), 'Specifying an external current in a homogeneous network calls for variables "addCurTime" and "addCurIext"'	# check for completeness of input
            else:
                assert (sum([(item == CurVal[0:3]).any() for item in self.ParaSim.keys()]) == 3), 'Specifying an external current in a heterogeneous network calls for variables "addCurTime", "addCurIext" and "addCurNeurons"'	# check for completeness of input
                assert((len(self.ParaSim['addCurTime']) == len(self.ParaSim['addCurIext'][0,:])) and (len(self.ParaSim['addCurNeurons']) == len(self.ParaSim['addCurIext'][:,0]))), 'the dimensions of the current array are not correct'	# check for consistency of input

        set_default(self.ParaSim,'pertSize',0)		# perturbation size only one

        #if self.ParaSim['pertSize'][0]:
        if self.ParaSim['pertSize'][0] > 0:
            self.ParaSim['pertSeed'] = [rnd.randint(1,2**32-1)]

        #try:
    	#self.ParaSim['pertVector']			# vector should be normalized afaik
    	#assert len(self.ParaSim['pertVector']) == self.Const['N'], 'The perturbation vector should have the same size as the network'
          #except:
    	#print "No perturbation direction given. PertVector is created randomly"
    	#vector_tmp = np.random.uniform(-1,1,self.Const['N'])
    	#vector_tmp = vector_tmp - np.linalg.norm(vector_tmp)
    	#set_default(self.ParaSim,'pertVector',list(vector_tmp/np.linalg.norm(vector_tmp)))	#direction of perturbation if pertSize is given

        set_default(self.ParaSim,'pertSpike',0)		# perturbed spike (i think first spike will be skipped)	only one
        set_default(self.ParaSim,'pertSynapse',0)		# perturbed synapse (i think first synapse of first spiking neuron will be skipped)		only one

        set_default(self.ParaSim,'distances',0)		# distance measurement
        set_default(self.ParaSim,'measures',0)		# should measurements be made?

        set_default(self.ParaSim,'CDB',0)			# fast break at convergence or divergence disabled
        set_default(self.ParaSim,'measureSpikes',[])		# measurements are made at times of here specified spikes (only one works)
        set_default(self.ParaSim,'measureTimes',[])		# measurements are made at here specified times (only one works)

        if self.ParaSim['measures'][0] or self.ParaSim['distances'][0] or self.ParaSim['CDB'][0]:
            assert (len(self.ParaSim['measureSpikes']) or len(self.ParaSim['measureTimes'])), 'If "measures", "distances" or "CDB" is given, at least one of "measuresSpikes" or "measuresTimes" has to be given.'

        set_default(self.ParaSim,'instPopRateBinSize',0)

        if self.ParaSim['CDB'] == 2:
            set_default(self.ParaSim,'D_decorr',100)

        ## Some sanity checks on the ParaSim input
        if self.ParaSim['instPopRateBinSize'][0]:
            assert self.ParaSim['TC'] > 0 , 'TC required for pop measurement'

        ## generate hash and filename
        if ('addCurIext' in self.ParaSim.keys()):
            if (len(self.ParaSim['addCurIext']) > 1000000): # for many addCurIext, take Hash of every e*1000th addIext
            	addCurIextReset = np.copy(self.ParaSim['addCurIext'])
            	idx = np.floor(np.linspace(0,len(self.ParaSim['addCurIext']),10000))
            	self.ParaSim['addCurIext_hash'] = self.ParaSim['addCurIext'][idx]
            	Hash = hashing(self.ParaSim,['CDB','rateWnt','pR','seedONS','train','ISIneurons','ISIstats','LyapunovExp','randomizedLEspectra','pLE','subN','subLyapunovExp','rateWntSubN','LyapunovExpConvergence','subLyapunovExpConvergence','SWONS','addCur','addCurHomo','addCurTime','addCurIext_hash','addCurNeurons','pertSize','pertVector','pertSpike','pertSynapse','pertSeed','distances','measures','measureSpikes','measureTimes','saveFinalState','CLV','SWCLV','instPopRateBinSize'])
            else:
                Hash = hashing(self.ParaSim,self.ParaSim.keys())  # for few addCurIext, take Hash of all Para.*
        else:
          #Hash = hashing(self.ParaSim,self.ParaSim.keys())  # for no addCurIext, take Hash of all Para.*
          #Hash = hashing(self.ParaSim,['CDB','TW','SW','TC','SC','TR','SR','pertSize','pertVector','measureTimes','saveFinalState','pertSpike','pertSynapse','distances','measures','rateWnt','rateWntSubN'])
            Hash = hashing(self.ParaSim,self.ParaSim.keys())
        self.Path['Sim'] = self.Path['inpara'] + 'ParaSimu-' + Hash + '.nc'
        # ParaSimulation goes in front to avoid ncdump error (bug)
        # ncdump: name begins with space or control-character: 1

        ## Open netCDF files.
        if not os.path.exists(self.Path['Sim']):
            if not suppressMessages:
                print('writing simulation netcdf file: %s' % self.Path['Sim'])
            ncid = netcdf.netcdf_file(self.Path['Sim'], 'w')

            ## Define the dimensions of the variables.
            ncid.createDimension('one', 1)
            ncid.createDimension('szISI', len(self.ParaSim['ISIneurons']))
            ncid.createDimension('szTrain', len(self.ParaSim['train']))
            ncid.createDimension('szISIdecorr', len(self.ParaSim['ISIdecorr']))

            if self.ParaSim['addCur'][0]:
            	ncid.createDimension('szCurNeurons', len(self.ParaSim['addCurNeurons']))
            	ncid.createDimension('szCurTime', len(self.ParaSim['addCurTime']))

        if self.ParaSim['addCurHomo'][0]:
            ncid.createDimension('szCurIext', len(self.ParaSim['addCurIext']))
        else:
            # save the 2D array with the additional current for each
            # neuron, as a 1D array, where each neuron's current is
            # appended after each other
            ncid.createDimension('szCurIext', len(self.ParaSim['addCurIext'])*len(self.ParaSim['addCurNeurons']))
            self.ParaSim['addCurIext'] = np.reshape(self.ParaSim['addCurIext'], 1, len(self.ParaSim['addCurIext'])*len(self.ParaSim['addCurNeurons']))

            if (self.ParaSim['measures'][0] or self.ParaSim['distances'][0] or self.ParaSim['CDB']):
            	ncid.createDimension('szMeasureSpikes', len(self.ParaSim['measureSpikes']))
            	ncid.createDimension('szMeasureTimes', len(self.ParaSim['measureTimes']))

          #if self.ParaSim['pertSize'][0]:
    	#ncid.createDimension('szPertVector', len(self.ParaSim['pertVector']))

        VarSim_pd = ncid.createVariable('pd', 'i4', ('one',))
        VarSim02 = ncid.createVariable('train', 'i4', ('szTrain',))
        VarSim03 = ncid.createVariable('ISIstats', 'i4', ('one',))
        VarSim04 = ncid.createVariable('ISIneurons', 'i4', ('szISI',))
        VarSim05 = ncid.createVariable('ISIbins', 'i4', ('one',))
        VarSim06 = ncid.createVariable('LyapunovExp', 'i4', ('one',))
        VarSim07 = ncid.createVariable('seedONS', 'i4', ('one',))
        VarSim08 = ncid.createVariable('SR', 'f4', ('one',))
        VarSim09 = ncid.createVariable('SW', 'f4', ('one',))
        VarSim10 = ncid.createVariable('SC', 'f4', ('one',))
        VarSim11 = ncid.createVariable('TR', 'f4', ('one',))
        VarSim12 = ncid.createVariable('TW', 'f4', ('one',))
        VarSim13 = ncid.createVariable('TC', 'f4', ('one',))
        VarSim15 = ncid.createVariable('rateWnt', 'f4', ('one',))
        VarSim16 = ncid.createVariable('pR', 'f4', ('one',))
        VarSim17 = ncid.createVariable('ONstep', 'i4', ('one',))
        VarSim18 = ncid.createVariable('saveFinalState', 'i4', ('one',))
        VarSim19 = ncid.createVariable('LyapunovExpConvergence', 'i4', ('one',))
        VarSim20 = ncid.createVariable('SWONS', 'f4', ('one',))
        VarSim21 = ncid.createVariable('pLE', 'f4', ('one',))
        VarSim22 = ncid.createVariable('CLV', 'i4', ('one',))
        VarSim23 = ncid.createVariable('SWCLV', 'i4', ('one',))
        VarSim24 = ncid.createVariable('instPopRateBinSize', 'f4', ('one',))
        VarSim25 = ncid.createVariable('subN', 'i4', ('one',))
        VarSim25a = ncid.createVariable('subLyapunovExp', 'f4', ('one',))
        VarSim25b = ncid.createVariable('subLyapunovExpConvergence', 'i4', ('one',))
        VarSim26 = ncid.createVariable('rateWntSubN', 'f4', ('one',))
        VarSim39 = ncid.createVariable('randomizedLEspectra', 'f4', ('one',))
        VarSim40 = ncid.createVariable('CDB', 'i4', ('one',))
        VarSim41 = ncid.createVariable('synchrony', 'i4', ('one',))
        VarSim42 = ncid.createVariable('ISIdecorr', 'i4', ('szISIdecorr',))

        VarSim_Cur = ncid.createVariable('addCur', 'i4', ('one',))

        if self.ParaSim['addCur'][0]:
        	VarSim_CurHomo = ncid.createVariable('addCurHomo', 'i4', ('one',))
        	VarSim_CurNeurons = ncid.createVariable('addCurNeurons', 'i4', ('szCurNeurons',))
        	VarSim_CurTime = ncid.createVariable('addCurTime', 'f4', ('szCurTime',))
        	VarSim_CurIext = ncid.createVariable('addCurIext', 'f4', ('szCurIext',))

        if self.ParaSim['CDB'][0] == 2:
            VarSim_D_decorr = ncid.createVariable('D_decorr','f4',('one',))

        VarSim_pertSize = ncid.createVariable('pertSize', 'f4', ('one',))
        VarSim_pertSpike = ncid.createVariable('pertSpike', 'i4', ('one',))
        VarSim_pertSynapse = ncid.createVariable('pertSynapse', 'i4', ('one',))

        VarSim_dist = ncid.createVariable('distances', 'i4', ('one',))
        VarSim_measures = ncid.createVariable('measures', 'i4', ('one',))

        if (self.ParaSim['measures'][0] or self.ParaSim['distances'][0] or self.ParaSim['CDB']):
            if len(self.ParaSim['measureSpikes']):
                VarSim_measureSpikes = ncid.createVariable('measureSpikes', 'i4', ('szMeasureSpikes',))
            elif len(self.ParaSim['measureTimes']):
                VarSim_measureTimes = ncid.createVariable('measureTimes', 'f4', ('szMeasureTimes',))

        if self.ParaSim['pertSize'][0] > 0:
        	VarSim_pertSeed = ncid.createVariable('pertSeed', 'i', ('one',))	# is this now the right datatype?
        	VarSim_pertSeed[:] = self.ParaSim['pertSeed']
    	#VarSim_pertVector = ncid.createVariable('pertVector', 'f4', ('szPertVector',))

        ## Write data to variable.
        VarSim_pd[:] = self.ParaSim['pd']
        VarSim02[:] = self.ParaSim['train']
        VarSim03[:] = self.ParaSim['ISIstats']
        VarSim04[:] = self.ParaSim['ISIneurons']
        VarSim05[:] = self.ParaSim['ISIbins']
        VarSim06[:] = self.ParaSim['LyapunovExp']
        VarSim07[:] = self.ParaSim['seedONS']
        VarSim08[:] = self.ParaSim['SR']
        VarSim09[:] = self.ParaSim['SW']
        VarSim10[:] = self.ParaSim['SC']
        VarSim11[:] = [self.ParaSim['TR'][0]*1000]
        VarSim12[:] = [self.ParaSim['TW'][0]*1000]
        VarSim13[:] = [self.ParaSim['TC'][0]*1000]
        VarSim15[:] = [self.ParaSim['rateWnt'][0]/1000.]
        VarSim16[:] = self.ParaSim['pR']
        VarSim17[:] = self.ParaSim['ONstep']
        VarSim18[:] = self.ParaSim['saveFinalState']
        VarSim19[:] = self.ParaSim['LyapunovExpConvergence']
        VarSim20[:] = self.ParaSim['SWONS']
        VarSim21[:] = self.ParaSim['pLE']
        VarSim22[:] = self.ParaSim['CLV']
        VarSim23[:] = self.ParaSim['SWCLV']
        VarSim24[:] = self.ParaSim['instPopRateBinSize']
        VarSim25[:] = self.ParaSim['subN']
        VarSim25a[:] = self.ParaSim['subLyapunovExp']
        VarSim25b[:] = self.ParaSim['subLyapunovExpConvergence']
        VarSim26[:] = [self.ParaSim['rateWntSubN'][0]/1000.]
        VarSim39[:] = self.ParaSim['randomizedLEspectra']
        VarSim40[:] = self.ParaSim['CDB']
        VarSim41[:] = self.ParaSim['synchrony']
        VarSim42[:] = self.ParaSim['ISIdecorr']

        VarSim_Cur[:] = self.ParaSim['addCur']

        if self.ParaSim['addCur'][0]:
        	VarSim_CurHomo[:] = self.ParaSim['addCurHomo']
        	VarSim_CurNeurons[:] = self.ParaSim['addCurNeurons']
        	VarSim_CurTime[:] = [self.ParaSim['addCurTime'][0]*1000]
        	VarSim_CurIext[:] = self.ParaSim['addCurIext']

        if self.ParaSim['CDB'][0] == 2:
        	VarSim_D_decorr[:] = self.ParaSim['D_decorr']
        	#print self.ParaSim['D_decorr']

        VarSim_pertSize[:] = self.ParaSim['pertSize']
        VarSim_pertSpike[:] = self.ParaSim['pertSpike']
        VarSim_pertSynapse[:] = self.ParaSim['pertSynapse']
        VarSim_dist[:] = self.ParaSim['distances']
        VarSim_measures[:] = self.ParaSim['measures']

        if (self.ParaSim['measures'][0] or self.ParaSim['distances'][0] or self.ParaSim['CDB'][0] or self.ParaSim['synchrony'][0]):
            if len(self.ParaSim['measureSpikes']):
                VarSim_measureSpikes[:] = self.ParaSim['measureSpikes']
            elif len(self.ParaSim['measureTimes']):
                VarSim_measureTimes[:] = [time*1000 for time in self.ParaSim['measureTimes']]

        #if self.ParaSim['pertSize'][0]:
	       #VarSim_pertVector[:] = self.ParaSim['pertVector']

        ## close file
        ncid.close()


############# read stuff for netcdf ##############

def readDataOut(fileName,suppressMessages=1):

    assert os.path.exists(fileName), '%s does not exist!' % fileName

    if not suppressMessages:
        print('reading result netcdf file: %s' % fileName)

    #######################
    ## IMPORT OF NETcdf
    #######################

    Data = {}
    # Get file ID:
    ncid = netcdf.netcdf_file(fileName,'r',mmap=False)

    # Get the value of the variables
    Data['SW'] = ncid.variables['SW'].getValue()
    Data['SC'] = ncid.variables['SC'].getValue()
    Data['TW'] = ncid.variables['TW'].getValue()
    Data['TC'] = ncid.variables['TC'].getValue()
    Data['rateW'] = ncid.variables['rateW'].getValue()
    Data['rateC'] = ncid.variables['rateC'].getValue()

    #read spike train if defined
    dim = 0
    try:
        dim = ncid.dimensions['finalStatesSz']
    except:
        if not suppressMessages:
            print('no final states and currents stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading final states and currents')
        Data['finalCurrents'] = ncid.variables['finalCurrents'][:]
        #print "curr read"
        finalState = ncid.variables['finalStates'][:]
        finalStateIdx = ncid.variables['finalStatesIdx'][:]

        #reshape the final state
        Data['finalStates'] = np.zeros(len(finalStateIdx))	#what dimension?

        for i in range(len(finalStateIdx)):
            Data['finalStates'][i] = finalState[finalStateIdx[i]]

    #Data['finalStates'][-1] = finalState[finalStateIdx[-1]:]

    #read spike train if defined
    dim = 0
    try:
        dim = ncid.dimensions['spikeTrainSz']
    except:
        if not suppressMessages:
            print('no spike train stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading spike train')

        Data['train'] = np.zeros((dim,2))
        Data['train'][:,0] = ncid.variables['trainNeuron'][:]
        Data['train'][:,1] = ncid.variables['trainTime'][:]

        if 'trainTimePert' in ncid.variables.keys():
            Data['trainPert'] = np.zeros((ncid.dimensions['PTSz'],dim,2))
            Data['trainPert'][:,:,0] = ncid.variables['trainNeuronPert'][:]
            Data['trainPert'][:,:,1] = ncid.variables['trainTimePert'][:]


    #read Lyapunov exponents if defined
    dimLE = 0
    try:
        dimLE = ncid.dimensions['LEonsSz']
    except:
        if not suppressMessages:
            print('no Lyapunov exponents stored in the netcdf file')

    if (dimLE > 0):
        if not suppressMessages:
            print('reading Lyapunov exponents')
        Data['LyapunovExponents'] = ncid.variables['LEons'][:]

        #read times at which the orthonormalization were done
        dimTimes = 0
        try:
            dimTimes = ncid.dimensions['LEtimesSz']
        except:
            if not suppressMessages:
                print('no Lyapunov exponent times stored in the netcdf file')

        if (dimTimes > 0):
            if not suppressMessages:
                print('reading Lyapunov exponents times')
            Data['LEtimes'] = ncid.variables['LEtimes'][:]

        #read Lyapunov exponents convergence if defined
        dimConv = 0
        try:
            dimConv = ncid.dimensions['LEconvergenceSz']
        except:
            if not suppressMessages:
                print('no Lyapunov exponents convergence stored in the netcdf file')

        if (dimConv > 0):
            if not suppressMessages:
                print('reading Lyapunov exponents convergence')
            Data['LEconvergence'] = ncid.variables['LEconvergence'][:]
            Data['LEconvergence'] = Data['LEconvergence'].reshape((dimConv/dimLE, dimLE))


        #read backward iteration Lyapunov exponents if defined
        dim = 0
        try:
            dim = ncid.dimensions['LEclvSz']
        except:
            if not suppressMessages:
                print('no backward iteration Lyapunov exponents stored in the netcdf file')

        if (dim > 0):
            if not suppressMessages:
                print('reading backward iteration Lyapunov exponents')
            Data['LEclv'] = ncid.variables['LEclv'][:]

        #read local Lyapunov exponents if defined
        dim = 0
        try:
            dim = ncid.dimensions['localLESz']
        except:
            if not suppressMessages:
                print('no local Lyapunov exponents stored in the netcdf file')

        if (dim > 0):
            if not suppressMessages:
                print('reading local Lyapunov exponents')
            Data['localLE'] = ncid.variables['localLE'][:]
            Data['localLE'] = Data['localLE'].reshape((dim/dimLE,dimLE))

    #read subLyapunov exponents if defined
    dimLE = 0
    try:
        dimLE = ncid.dimensions['subLEonsSz']
    except:
        if not suppressMessages:
            print('no subLyapunov exponents stored in the netcdf file')

    if (dimLE > 0):
        if not suppressMessages:
            print('reading subLyapunov exponents')
        Data['subLyapunovExponents'] = ncid.variables['subLEons'][:]

        #read Lyapunov exponents convergence if defined
        dimConv = 0
        try:
            dimConv = ncid.dimensions['subLEconvergenceSz']
        except:
            if not suppressMessages:
                print('no Lyapunov exponents convergence stored in the netcdf file')

        if (dimConv > 0):
            if not suppressMessages:
                print('reading Lyapunov exponents convergence')
            Data['subLEconvergence'] = ncid.variables['subLEconvergence'][:]
            Data['subLEconvergence'] = Data['subLEconvergence'].reshape((dimConv/dimLE, dimLE))

    #read firing rates if defined
    dim = 0
    try:
        dim = ncid.dimensions['rateNeuronsSz']
    except:
        if not suppressMessages:
            print('no firing rates stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading firing rates')
        Data['rateNeurons'] = ncid.variables['rateNeurons'][:]

    #read rate distribution
    dim = 0
    try:
        dim = ncid.dimensions['rateDistSz']
    except:
        if not suppressMessages:
            print('no rate distribution stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading rate distribution')
        Data['rateDistX'] = ncid.variables['rateDistX'][:]
        Data['rateDistY'] = ncid.variables['rateDistY'][:]

    #read coefficients of variation if defined
    dim = 0
    try:
        dim = ncid.dimensions['cvNeuronsSz']
    except:
        if not suppressMessages:
            print('no coefficient of variation stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading coefficients of variation')
        Data['cvNeurons'] = ncid.variables['cvNeurons'][:]

    #read cv distribution
    dim = 0
    try:
        dim = ncid.dimensions['cvDistSz']
    except:
        if not suppressMessages:
            print('no coefficient of variation distribution stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading coefficient of variation distribution')
        Data['cvDistX'] = ncid.variables['cvDistX'][:]
        Data['cvDistY'] = ncid.variables['cvDistY'][:]

    #read skewness if defined
    dim = 0
    try:
        dim = ncid.dimensions['skewnessNeuronsSz']
    except:
        if not suppressMessages:
            print('no skewness stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading skewness')
        Data['skewnessNeurons'] = ncid.variables['skewnessNeurons'][:]

    #read skewness distribution
    dim = 0
    try:
        dim = ncid.dimensions['skewnessDistSz']
    except:
        if not suppressMessages:
            print('no skewness distribution stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading skewness distribution')
        Data['skewnessDistX'] = ncid.variables['skewnessDistX'][:]
        Data['skewnessDistY'] = ncid.variables['skewnessDistY'][:]

    #read kurtosis if defined
    dim = 0
    try:
        dim = ncid.dimensions['kurtosisNeuronsSz']
    except:
        if not suppressMessages:
            print('no kurtosis stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading kurtosis')
        Data['kurtosisNeurons'] = ncid.variables['kurtosisNeurons'][:]

    #read cv distribution
    dim = 0
    try:
        dim = ncid.dimensions['kurtosisDistSz']
    except:
        if not suppressMessages:
            print('no kurtosis distribution stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading kurtosis distribution')
        Data['kurtosisDistX'] = ncid.variables['kurtosisDistX'][:]
        Data['kurtosisDistY'] = ncid.variables['kurtosisDistY'][:]

    try:
        Data['chi'] = ncid.variables['chi'].getValue()
    except:
        if not suppressMessages:
            print('no synchrony measurement included')

    #read phase variables
    dim = 0
    try:
        dim = ncid.dimensions['phaseNeuronsSz']
    except:
        if not suppressMessages:
            print('no phases stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading phases')
        Data['phaseTimes'] = ncid.variables['phaseTimes'][:]
        Data['phaseNeurons'] = ncid.variables['phaseNeurons'][:]

        Data['phaseNeurons'] = Data['phaseNeurons'].reshape((len(Data['phaseTimes'],len(Data['phaseNeurons'])/len(Data['phaseTimes']))))

    #read distance
    dim = 0
    try:
        dim = ncid.dimensions['distancesSz']
    except:
        if not suppressMessages:
            print('no distances stored in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading distances')
        Data['measureTimes'] = ncid.variables['measureTimes'][:]		#attention! changed 'phaseTimes' -> 'measureTimes'. should phT be in the program somewhere?
        Data['distances'] = ncid.variables['distances'][:]

    # read states
    dim = 0
    try:
        dim = ncid.dimensions['measureTimesSz']
    except:
        if not suppressMessages:
            print('no measure Times specified in the netcdf file')

    if (dim > 0):
        if not suppressMessages:
            print('reading measurements')
        Data['measureTimes'] = ncid.variables['measureTimes'][:]		#attention! changed 'phaseTimes' -> 'measureTimes'. should phT be in the program somewhere?
        #try:
        Data['measureStates'] = ncid.variables['measureStates'][:]

        if 'measureStatesPert' in ncid.variables.keys():
            Data['measureStatesPert'] = ncid.variables['measureStatesPert'][:]

            #tmp_measures = ncid.variables['measure1stStateVarNeurons'][:]
            #Data['measurePhases'] = tmp_measures.reshape(dim,len(tmp_measures)/dim)
        #except:
            #if not suppressMessages:
        #print 'no measurements specified in the netcdf file'

        try:
            Data['PTtimes'] = ncid.variables['PTtimes'][:]
            Data['PTdistances'] = ncid.variables['PTdistances'][:]
        except:
            if not suppressMessages:
                print('no convergence/divergence break included')



    dim = 0
    try:
        Data['finalDistance'] = ncid.variables['finalDistance'].getValue()
    except:
        if not suppressMessages:
            print('no final distance stored')

    dim = 0
    try:
        dim = ncid.dimensions['ISIdecorrSz']
    except:
        if not suppressMessages:
            print('no ISI of potential decorrelation events stored')

    if (dim > 0):
        if not suppressMessages:
            print('reading ISIdecorr times')
        Data['ISIdecorr'] = ncid.variables['mean_ISIdecorr'][:].transpose()
        Data['cvISIdecorr'] = ncid.variables['cv_ISIdecorr'][:].transpose()


    ncid.close()

    return Data
