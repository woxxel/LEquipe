/*
 * LEinputoutput.cpp
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */
#include "LEinputoutput.h"

LE_input_output::LE_input_output()
{
	//! initialize parallel environment
	nP = 1;																	//!< number of processors involved
	myID = 0;																//!< local ID of the processors, for unique communications

#ifdef PAR
	//! Initialize the parallel environment.
	MPI_Comm_size(MPI_COMM_WORLD,&nP);
	MPI_Comm_rank(MPI_COMM_WORLD,&myID);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

LE_input_output::~LE_input_output()
{

}

void LE_input_output::read_netcdf(string fileNeurons, string fileTopology, string fileSimulation, string filePuppet, string fileRefTrajectory, st_in* in)
{
	//! This reads in the netcdf file\a netcdf_filename and assigns the values in the structure\a in to run the simulation with.

	/** ----------------------------
	 *   open file ParaInNeurons.nc
	 */

	cout << "reading " << fileNeurons << " ..." << endl;

	NcFile ParaInNeurons(fileNeurons.c_str(), NcFile::ReadOnly);

	//! Retrieve pointers to NEURON variables in netcdf file
	NcVar *NP = ParaInNeurons.get_var("N");
	NcVar *NeuronTypeP = ParaInNeurons.get_var("NeuronType");
	NcVar *tauMP = ParaInNeurons.get_var("tauM");
	NcVar *IextP = ParaInNeurons.get_var("Iext");

		
	//! Check whether the neuron parameters are homogenous throughout the network 
	//NOTE: it is enough to check whether the size of NeuronType is larger 1, since all neuron parameters are made same size in writeNet.m 
	int NeuronTypeSz =  NeuronTypeP->num_vals();
	int homogNet;
	NeuronTypeSz<2 ? homogNet = 1 : homogNet  = 0;
	in->homogNet = homogNet;
   	
	//! Read all the values from the variables into memory.
	int Nall, Nloc, Noffset;

	NP->get(&Nall, 1);

#ifdef PAR
	//! Number of neurons processed on each node. (integer division)
	Nloc = (int) floor(Nall/nP);

	//! The offset for each node for the global neuron index.
	Noffset = myID*Nloc;

	// If the number of all neurons cannot be divided equally on all processors, the last one will take the left over neurons.
	if (myID == nP - 1) Nloc = Nall - Noffset;
	
	

#else

	Nloc = Nall;
	Noffset = 0;

#endif

	in->Nall = Nall;
	in->Nloc = Nloc;
	in->Noffset = Noffset;
	in->Ndrive = 0;


	//! Depending on whether the network is homogeneous, the arrays describing the neurons in the network are either 1 (homogeneous) or N (heterogeneous)-dimensional
	int Nhomo = (in->homogNet) ? 1 : Nall;

	in->neuronType.resize(Nhomo);
	in->tauM.resize(Nhomo);
	in->Iext.resize(Nhomo);
	

	NeuronTypeP->get(&in->neuronType.front(), Nhomo);
	tauMP->get(&in->tauM.front(), Nhomo);
	vector<double> tmp=vector<double>(Nhomo);
	IextP->get(&tmp.front(), Nhomo);
	copy(tmp.begin(),tmp.end(),in->Iext.begin());
	
	//! check if at least one of the neurons is a poisson neuron, and if so read in the poisson rates for all neurons (even though, not all of them might be poisson neurons)
	bool poisson = false;
	for (int n=0; n<Nhomo; n++)
		if ((in->neuronType[n] == 5) || (in->neuronType[n] == 51)|| (in->neuronType[n] == 16))
		{
			poisson = true;
			break;
		}
		
	if (poisson)
	{
		NcDim* NdriveP = ParaInNeurons.get_dim("Ndrive");
		in->Ndrive = NdriveP->size();
		
		in->poissonrate.resize(Nhomo);
		NcVar *poissonrateP = ParaInNeurons.get_var("poissonrate");
		poissonrateP->get(&in->poissonrate.front(), Nhomo);
	}
	

	//! check if at least one of the neurons is a gamma neuron, and if so read in the gammashape rates for all neurons (even though, not all of them might be poisson neurons)
	bool gamma = false;
	for (int n=0; n<Nhomo; n++)
		if (in->neuronType[n] == 6)
		{
			gamma = true;
			break;
		}
  
	if (gamma)
	{
		in->gammashape.resize(Nhomo);
		NcVar *gammashapeP = ParaInNeurons.get_var("gammashape");
		gammashapeP->get(&in->gammashape.front(), Nhomo);
	}


	//! check if at least one of the neurons is a puppet neuron, and if so read in the given spike train
	bool puppet = false;			// added by Alex (puppet)
	for (int n=0; n<Nhomo; n++)
		if (in->neuronType[n] == 7)
		{
			puppet = true;
			break;
		}
		
	if (puppet)
	{
	  	/** -----------------------------
		  *   open file ParaInPuppet.nc
		*/

		cout << "reading " << filePuppet << " ..." << endl;

		NcFile ParaInPuppet(filePuppet.c_str(), NcFile::ReadOnly);

		NcDim* NdriveP = ParaInPuppet.get_dim("Ndrive");
		in->Ndrive = NdriveP->size();
		
		NcDim* pTrainSzP =  ParaInPuppet.get_dim("puppetTrainSz");
		int  pTrainSz = pTrainSzP->size();
		
		in->puppetTrain.resize(in->Ndrive,vector<reell> (pTrainSz));

		NcVar *pTrainP = ParaInPuppet.get_var("puppetTrain");

		for (int n=0; n<in->Ndrive; n++)
		{
			pTrainP->set_cur(n,0);
			pTrainP->get(&in->puppetTrain[n][0],1,pTrainSz);
		}
	}
	
	

	//! check if at least one of the neurons is a theta neuron, and if so read in the AP rapidness values for all neurons (even though, not all of them might be rapid theta neurons)
	bool rapidtheta = false;
	for (int n=0; n<Nhomo; n++)
		if ((in->neuronType[n] == 1) || (in->neuronType[n] == 51) || (in->neuronType[n] == 16))
		{
			rapidtheta = true;
			break;
		}

	if (rapidtheta)
	{
		in->rapidness.resize(Nhomo);
		NcVar *rapidnessP = ParaInNeurons.get_var("rapidness");
		rapidnessP->get(&in->rapidness.front(), Nhomo);
	}

	//! check if at least one of the neurons is a type1type2 neuron, and if so read in the parameter a for all neurons (even though, not all of them might be of this kind)
	bool type1type2 = false;
	for (int n=0; n<Nhomo; n++)
		if (in->neuronType[n] == 3)
		{
			type1type2 = true;
			break;
		}
	
	if (type1type2)
	{
		NcVar *paraP = ParaInNeurons.get_var("type1type2Para");
		in->type1type2Para.resize(Nhomo);
		paraP->get(&in->type1type2Para.front(), Nhomo);
	}


	//! check if at least one neuron is a twoDlinear model neuron, and if so read in the parameter a for all neurons (even though, not all of them might be of this kind)
	bool twoDlinear = false;
	for (int n=0; n<Nhomo; n++)
		if (in->neuronType[n] == 10)
		{
			twoDlinear = true;
			break;
		}
	
	if (twoDlinear)
	{
		vector<double> alpha(Nhomo);
		vector<double> beta(Nhomo);
		vector<double> gamma(Nhomo);
		vector<double> delta(Nhomo);
		vector<double> Cw(Nhomo);
		vector<double> tauS(Nhomo);

		NcVar *paraP = ParaInNeurons.get_var("twoDlinear_alpha");
		paraP->get(&alpha.front(), Nhomo);

		paraP = ParaInNeurons.get_var("twoDlinear_beta");
		paraP->get(&beta.front(), Nhomo);

		paraP = ParaInNeurons.get_var("twoDlinear_gamma");
		paraP->get(&gamma.front(), Nhomo);

		paraP = ParaInNeurons.get_var("twoDlinear_delta");
		paraP->get(&delta.front(), Nhomo);

		paraP = ParaInNeurons.get_var("twoDlinear_Cw");
		paraP->get(&Cw.front(), Nhomo);

		paraP = ParaInNeurons.get_var("twoDlinear_tauS");
		paraP->get(&tauS.front(), Nhomo);

		in->twoDlinearParas = vector<st_twoDlinear>(Nhomo);
		for (int n=0; n<Nhomo; n++)
		{
			in->twoDlinearParas[n].alpha = alpha[n];
			in->twoDlinearParas[n].beta = beta[n];
			in->twoDlinearParas[n].gamma = gamma[n];
			in->twoDlinearParas[n].delta = delta[n];
			in->twoDlinearParas[n].Cw = Cw[n];
			in->twoDlinearParas[n].tauS = tauS[n];
		}
	}
	
	in->Nrec = in->Nall - in->Ndrive;

	//! get the initial conditions of all neurons
	NcDim* stateSzP =  ParaInNeurons.get_dim("initStatesSz");
	int  stateSz = stateSzP->size();

	NcDim* stateIdxSzP =  ParaInNeurons.get_dim("initStatesIdxSz");
	int  stateIdxSz = stateIdxSzP->size();


	NcVar *stateP = ParaInNeurons.get_var("initStates");
	vector<double> state(stateSz);
	stateP->get(&state.front(), stateSz);

	NcVar *stateIdxP = ParaInNeurons.get_var("initStatesIdx");
	vector<double> stateIdx(stateIdxSz);
	stateIdxP->get(&stateIdx.front(), stateIdxSz);
	
	if (Nall != (int)stateIdx.size())
	{
		cout << "The initial conditions were not provided for all neurons!" << endl;
		throw(1);

	}

	in->init.resize(Nall);
	for (int n=0; n<Nall-1; n++)
		in->init[n].assign(state.begin() + stateIdx[n], state.begin() + stateIdx[n+1]);
	in->init[Nall-1].assign(state.begin() + stateIdx[Nall-1], state.end());
	

	//! if in parallel mode, only store the local neurons
#ifdef PAR
	if (!in->homogNet)
	{
		in->neuronType.assign(in->neuronType.begin() + Noffset, in->neuronType.begin() + Noffset + Nloc);
		in->tauM.assign(in->tauM.begin() + Noffset, in->tauM.begin() + Noffset + Nloc);
		in->Iext.assign(in->Iext.begin() + Noffset, in->Iext.begin() + Noffset + Nloc);

		if (rapidtheta)
			in->rapidness.assign(in->rapidness.begin() + Noffset, in->rapidness.begin() + Noffset + Nloc);

		if (type1type2)
			in->type1type2Para.assign(in->type1type2Para.begin() + Noffset, in->type1type2Para.begin() + Noffset + Nloc);

		if (twoDlinear)
			in->twoDlinearParas.assign(in->twoDlinearParas.begin() + Noffset, in->twoDlinearParas.begin() + Noffset + Nloc);
	}

	in->init.assign(in->init.begin() + Noffset, in->init.begin() + Noffset + Nloc);

#endif
	

	/** -----------------------------
	 *   open file ParaInTopology.nc
	 */

	cout << "reading " << fileTopology << " ..." << endl;

	NcFile ParaInTopology(fileTopology.c_str(), NcFile::ReadOnly);

	NcVar* postP = ParaInTopology.get_var("postSyn");
	NcVar* JP = ParaInTopology.get_var("J");
	NcVar* pSynP = ParaInTopology.get_var("pSyn");
	NcVar* row_lengthP = ParaInTopology.get_var("outDegree");

	NcDim* max_elementsP =  ParaInTopology.get_dim("synapses");
	int  max_elements = max_elementsP->size();

	vector<int> post(max_elements);
	vector<int> rowLength(Nall);
	vector<double> J(1);
	vector<double> pSyn(1);

	postP -> get(&post.front(), max_elements);
	row_lengthP -> get(&rowLength.front(), Nall);

	//! Check whether the synaptic couplings and release probabilities are homogeneous throughout the network.
	int homogSyn = 0;
	int pSynSz =  pSynP->num_vals();
	pSynSz<2? homogSyn = 1 : homogSyn  = 0;

	if (homogSyn)
	{
		JP -> get(&J.front(),1);
		pSynP -> get(&pSyn.front(),1);
	}
	else
	{
		J.resize(max_elements);
		JP -> get(&J.front(), max_elements);

		pSyn.resize(max_elements);
		pSynP -> get(&pSyn.front(), max_elements);
	}


	in->synapses = vector<vector<st_synapse> > (Nall);
	st_synapse synTmp;
	int ind = 0;
	for (int n=0; n<Nall; n++)
		for (int e=0; e<rowLength[n]; e++)
		{
			//! if in parallel mode, only store the local neurons
			if ((post[ind] >= Noffset) && (post[ind] < Noffset + Nloc))
			{
				synTmp.post = post[ind] - Noffset;
				synTmp.cplg = (homogSyn) ? J[0] : J[ind];
				synTmp.prob = (homogSyn) ? pSyn[0] : pSyn[ind];

				in->synapses[n].push_back(synTmp);
			}

			ind++;
		}



	/** -------------------------------
	 *   open file ParaInSimulation.nc
	 */

	cout << "reading " << fileSimulation << " ..." << endl;

	NcFile ParaInSimulation(fileSimulation.c_str(), NcFile::ReadOnly);
	
	//! read in the number of spikes (S) or the time in ms (T) for rate finding (R), the warmup (W) and the calculation (C)
	NcVar *SRP = ParaInSimulation.get_var("SR");
	NcVar *TRP = ParaInSimulation.get_var("TR");
	NcVar *SWP = ParaInSimulation.get_var("SW");
	NcVar *TWP = ParaInSimulation.get_var("TW");
	NcVar *SCP= ParaInSimulation.get_var("SC");
	NcVar *TCP= ParaInSimulation.get_var("TC");

	double spikesPerNeuron;

	TRP -> get(&in->TR, 1);
	SRP -> get(&spikesPerNeuron, 1);
	in->SR = ceil(spikesPerNeuron*Nall);

	TWP -> get(&in->TW, 1);
	SWP -> get(&spikesPerNeuron, 1);
	in->SW = ceil(spikesPerNeuron*Nall);

	TCP -> get(&in->TC, 1);
	SCP -> get(&spikesPerNeuron, 1);
	in->SC = ceil(spikesPerNeuron*Nall);
	
	NcVar *pdP = ParaInSimulation.get_var("pd");
	pdP -> get(&in->pertDirections,1);

	//! read the wanted average firing rate in the network, if greater than 0, this rate will be tried to be reached by varying the external currents
	NcVar *rateP = ParaInSimulation.get_var("rateWnt");
	rateP->get(&in->rateWnt, 1);
	NcVar *ratePsub = ParaInSimulation.get_var("rateWntSubN");
	ratePsub->get(&in->rateWntSubN, 1);

	//! and the precision of the wanted firing rate
	NcVar *pRP = ParaInSimulation.get_var("pR");
	pRP->get(&in->pR, 1);

	//! Number of Lyapunov exponents to be calculated
	NcVar *LyapunovExpP = ParaInSimulation.get_var("LyapunovExp");
	LyapunovExpP -> get(&in->LyapunovExponents, 1);
	
	if (in->LyapunovExponents)
	{
		//! initial seed for the ONS in the Lyapunov exponent calculation
		NcVar *seedONSP= ParaInSimulation.get_var("seedONS");
		seedONSP-> get(&in->seedONS, 1);

		//! number of spikes per neurons for ONS warmup
		NcVar *SWONSP= ParaInSimulation.get_var("SWONS");
		SWONSP -> get(&spikesPerNeuron, 1);
		in->SWONS = ceil(spikesPerNeuron*Nall);

		//! initial step size of the orthonormalizations in the Lyapunov exponent calculation
		NcVar *ONstepP= ParaInSimulation.get_var("ONstep");
		ONstepP-> get(&in->ONstep, 1);

		NcVar *pLEP= ParaInSimulation.get_var("pLE");
		pLEP-> get(&in->pLEONS, 1);

		//! Should the convergence of the Lyapunov exponents be saved?
		//! The result file will be quite large, so the default should be off. On the other hand it is important to check the convergence.
		NcVar *LEconvP= ParaInSimulation.get_var("LyapunovExpConvergence");
		LEconvP->get(&in->LyapunovExponentsConvergence, 1);


		//! Parameters to calculate the covariant Lyapunov vectors
		NcVar *CLVP = ParaInSimulation.get_var("CLV");
		CLVP -> get(&in->CLV, 1);

		if (in->CLV)
		{
			NcVar *SWCLVP = ParaInSimulation.get_var("SWCLV");
			SWCLVP -> get(&spikesPerNeuron, 1);
			in->SWCLV = ceil(spikesPerNeuron*Nall);
		}
	
	}
	
	//!< Number of randomized state sequences for spectra calculation
	NcVar *randomizedLEspectraP = ParaInSimulation.get_var("randomizedLEspectra");
	randomizedLEspectraP -> get(&in->randomizedLEspectra, 1);

	
	int subNall, subNloc, subNoffset;
	//! Subnetwork size on which to compute LE stuff 
	NcVar *subNet = ParaInSimulation.get_var("subN");
	subNet -> get(&subNall, 1);

#ifdef PAR

	//! Number of subNetwork neurons processed on each node.
	if (myID*Nloc < subNall)
	{
		subNloc = Nloc;
		subNoffset = myID*Nloc;
	}
	else if ( (myID*Nloc < subNall) & (myID*(Nloc+1) > subNall))
	{
		subNloc = subNall - myID*Nloc;
		subNoffset = myID*Nloc;
	}
	else	//higher processors with no subNetwork neurons
	{
		subNloc = 0;
		subNoffset = subNall;
	}
#else

	subNloc = subNall;
	subNoffset = 0;

#endif
	in->subNall = subNall;
	in->subNloc = subNloc;
	in->subNoffset = subNoffset;
	
	//! Number of subNetwork Lyapunov exponents to be calculated
	NcVar *subLyapunovExpP = ParaInSimulation.get_var("subLyapunovExp");
	subLyapunovExpP -> get(&in->subLyapunovExponents, 1);
	
	if (in->subLyapunovExponents)
	{
		//! initial seed for the ONS in the Lyapunov exponent calculation
		NcVar *seedONSP= ParaInSimulation.get_var("seedONS");
		seedONSP-> get(&in->seedONS, 1);

		//! number of spikes per neurons for ONS warmup
		NcVar *SWONSP= ParaInSimulation.get_var("SWONS");
		SWONSP -> get(&spikesPerNeuron, 1);
		in->SWONS = ceil(spikesPerNeuron*Nall);

		//! initial step size of the orthonormalizations in the Lyapunov exponent calculation
		NcVar *ONstepP= ParaInSimulation.get_var("ONstep");
		ONstepP-> get(&in->ONstep, 1);

		NcVar *pLEP= ParaInSimulation.get_var("pLE");
		pLEP-> get(&in->pLEONS, 1);

		//! Should the convergence of the subLyapunov exponents be saved?
		//! The result file will be quite large, so the default should be off. On the other hand it is important to check the convergence.
		NcVar *subLEconvP= ParaInSimulation.get_var("subLyapunovExpConvergence");
		subLEconvP->get(&in->subLyapunovExponentsConvergence, 1);

		//! Parameters to calculate the covariant Lyapunov vectors
		NcVar *CLVP = ParaInSimulation.get_var("CLV");
		CLVP -> get(&in->CLV, 1);

		if (in->CLV)
		{
			NcVar *SWCLVP = ParaInSimulation.get_var("SWCLV");
			SWCLVP -> get(&spikesPerNeuron, 1);
			in->SWCLV = ceil(spikesPerNeuron*Nall);
		}
		
		//!< Number of randomized state sequences for spectra calculation
		NcVar *randomizedLEspectraP = ParaInSimulation.get_var("randomizedLEspectra");
		randomizedLEspectraP -> get(&in->randomizedLEspectra, 1);
	
	}

	//! save the final state of the simulation for other calculations
	NcVar *finalP= ParaInSimulation.get_var("saveFinalState");
	finalP-> get(&in->saveFinalState, 1);
	
	//! The neurons to calculate the inter spike interval (ISI) of.
	NcDim* szISIP =  ParaInSimulation.get_dim("szISI");
	int  szISI = szISIP->size();
	NcVar *ISIneuronsP = ParaInSimulation.get_var("ISIneurons");
	in->ISIneurons.resize(szISI);
	ISIneuronsP -> get(&in->ISIneurons.front(), szISI);
	
	//! If the provided neurons have index -1, then clear the array
	if (in->ISIneurons[0] == -1) in->ISIneurons.clear();

	//! The 'kinda' moments of the inter spike intervals (1=mean, 2=coeff. of variation, 3=skewness, 4=kurtosis)
	NcVar *ISIstatsP = ParaInSimulation.get_var("ISIstats");
	ISIstatsP -> get(&in->ISIstats, 1);

	//! The number of bins for the firing rate distribution
	NcVar *ISIbinsP = ParaInSimulation.get_var("ISIbins");
	ISIbinsP -> get(&in->ISIbins, 1);
	
	//! The neurons to calculate ISIdecorr statistics of (time between successive post-/presynaptic spikes)
	NcDim* szISIdecorrP =  ParaInSimulation.get_dim("szISIdecorr");
	int  szISIdecorr = szISIdecorrP->size();

	NcVar *ISIdecorrP = ParaInSimulation.get_var("ISIdecorr");
	in->ISIdecorr.resize(szISIdecorr);
	ISIdecorrP -> get(&in->ISIdecorr.front(), szISIdecorr);

	//! If the provided neurons have index -1, then don't calculate the spike train
	in->ISIdecorr_stats = (in->ISIdecorr[0] == -1) ? 0 : 1;
	
	//! The neurons to calculate the spike train of (size of the train vector).
	NcDim* szTrainP =  ParaInSimulation.get_dim("szTrain");
	int  szTrain = szTrainP->size();

	NcVar *trainP = ParaInSimulation.get_var("train");
	in->train.resize(szTrain);
	trainP -> get(&in->train.front(), szTrain);
	in->calc_train = (in->train[0] == -1) ? false : true;
	//! If the provided neurons have index -1, then don't calculate the spike train
	if (!in->calc_train) in->train.clear();
	
	//! calculate the synchrony measure chi
	NcVar *syncP= ParaInSimulation.get_var("synchrony");
	syncP-> get(&in->synchrony, 1);


	//! Time varying additional current stuff:
	NcVar *addCurP = ParaInSimulation.get_var("addCur");
	addCurP -> get(&in->addCur, 1);

	if (in->addCur)
	{
		NcVar *addCurHomoP = ParaInSimulation.get_var("addCurHomo");
		addCurHomoP -> get(&in->addCurHomo, 1);

		NcDim* szCurNeuronsP =  ParaInSimulation.get_dim("szCurNeurons");
		int  szCurNeurons = szCurNeuronsP->size();

		if (szCurNeurons == 0)
		{
			cout << "No neurons were specified that receive an additional current!" << endl;
			throw(1);
		}

		NcVar *addCurNeuronsP = ParaInSimulation.get_var("addCurNeurons");
		in->addCurNeurons.resize(szCurNeurons);
		addCurNeuronsP -> get(&in->addCurNeurons.front(), szCurNeurons);


		NcDim* szCurTimeP =  ParaInSimulation.get_dim("szCurTime");
		int  szCurTime = szCurTimeP->size();

		if (szCurTime == 0)
		{
			cout << "No times for the additional currents were specified!" << endl;
			throw(1);
		}

		NcVar *addCurTimeP = ParaInSimulation.get_var("addCurTime");
		in->addCurTimes = vector<double> (szCurTime);
		addCurTimeP -> get(&in->addCurTimes.front(), szCurTime);

		NcDim* szCurIextP =  ParaInSimulation.get_dim("szCurIext");
		int  szCurIext = szCurIextP->size();

		if (szCurIext == 0)
		{
			cout << "No additional currents were specified!" << endl;
			throw(1);
		}

		NcVar *addCurIextP = ParaInSimulation.get_var("addCurIext");
  	  	vector<double> IextTmp(szCurIext);
		addCurIextP -> get(&IextTmp.front(), szCurIext);

		if (in->addCurHomo)
		{
			if (szCurIext != szCurTime)
			{
				cout << "The additional current is not saved with the correct size in the netcdf file!" << endl;
				throw(1);
			}

			in->addCurIext = vector<vector<double> > (szCurTime, vector<double> (1));
			for (int t=0; t<szCurTime; t++)
				in->addCurIext[t][0] = IextTmp[t];

#ifdef PAR
			//! extract the neurons that belong to this processor in parallel mode
			vector<int> tmp(0);

			for (unsigned n=0; n<in->addCurNeurons.size(); n++)
				if ((in->addCurNeurons[n] >= Noffset) && (in->addCurNeurons[n] < Noffset + Nloc))
					tmp.push_back(in->addCurNeurons[n]);

			in->addCurNeurons = tmp;
#endif

		}
		else
		{
			if (szCurIext != szCurTime*szCurNeurons)
			{
				cout << "The additional current was not saved with the correct size in the netcdf file!" << endl;
				throw(1);
			}

			in->addCurIext = vector<vector<double> > (szCurTime, vector<double> (szCurNeurons));
			for (int t=0; t<szCurTime; t++)
				for (int n=0; n<szCurNeurons; n++)
					in->addCurIext[t][n] = IextTmp[n + t*szCurNeurons];

#ifdef PAR
			//! extract the neurons that belong to this processor in parallel mode (erase the other ones)
			//! extract the neurons that belong to this processor in parallel mode
			vector<int> tmp(0);
			vector<vector<double> > tmp2 (szCurTime, vector<double> (0));

			for (unsigned n=0; n<in->addCurNeurons.size(); n++)
				if ((in->addCurNeurons[n] >= Noffset) && (in->addCurNeurons[n] < Noffset + Nloc))
				{
					tmp.push_back(in->addCurNeurons[n]);
					for (int t=0; t<szCurTime; t++)
						tmp2[t].push_back(in->addCurIext[t][n]);
				}

			in->addCurNeurons = tmp;
			in->addCurIext = tmp2;

			/* old version, shouldn't worl, because the size changes during the loop and then each n after a deletion is skipped
			for (unsigned n=0; n<in->addCurNeurons.size(); n++)
				if ((in->addCurNeurons[n] < Noffset) || (in->addCurNeurons[n] >= Noffset + Nloc))
				{
					in->addCurNeurons.erase(in->addCurNeurons.begin() + n);

					for (int t=0; t<szCurTime; t++)
						in->addCurIext[t].erase(in->addCurIext[t].begin() + n);
				}
			*/
#endif

		}
	}

// 	//! The perturbation measurement stuff.
// 	NcVar *pertSizeP = ParaInSimulation.get_var("pertSize");
// 	pertSizeP -> get(&in->pertSize, 1);
	//! get perturbation size
	NcVar *pertSizeP = ParaInSimulation.get_var("pertSize");
	pertSizeP -> get(&in->pertSize, 1);
	
	cout << "pertSize: " << in->pertSize << endl;
	if (in->pertSize > 0)
	{
		NcVar *pertSeedP = ParaInSimulation.get_var("pertSeed");
		pertSeedP -> get(&in->pertSeed, 1);
		
		//! initialize random number generator
		gsl_rng_env_setup();
		const gsl_rng_type *TYPE = gsl_rng_mt19937;
		in->rng = gsl_rng_alloc(TYPE);
		gsl_rng_set(in->rng, in->pertSeed);
// 		cout << "pertSeed: " << in->pertSeed << endl;
	}
	
	
	//! The spike perturbation.
	NcVar *pertSpikeP = ParaInSimulation.get_var("pertSpike");
	pertSpikeP -> get(&in->pertSpike, 1);
	
	//! The synapse perturbation.
	NcVar *pertSynP = ParaInSimulation.get_var("pertSynapse");
	pertSynP -> get(&in->pertSynapse, 1);
	
	if ((in->pertSize && in->pertSpike) || (in->pertSize && in->pertSynapse) || (in->pertSynapse && in->pertSpike))
	{
		cout << "Only one type of perturbation can be applied at a time!" << endl;
		throw(1);
	}
	
	
// 	if (in->pertSize)
// 	{
// 		//!read perturbation vector
// 		NcVar *pertVectorP = ParaInSimulation.get_var("pertVector");
// 		vector<double> pertstate(stateSz);
// 		pertVectorP->get(&pertstate.front(), stateSz);
// 		in->pertVector.resize(Nall);
// 		for (int n=0; n<Nall-1; n++)
// 			in->pertVector[n].assign(pertstate.begin() + stateIdx[n], pertstate.begin() + stateIdx[n+1]);
// 		in->pertVector[Nall-1].assign(pertstate.begin() + stateIdx[Nall-1], pertstate.end());
// 		
// #ifdef PAR
// 			//! extract the neurons that belong to this processor in parallel mode (erase the other ones)
// 			in->pertVector.assign(in->pertVector.begin() + Noffset, in->pertVector.begin() + Noffset + Nloc);
// #endif
// 	}
	

	//! Which if any measurable should be saved?
	NcVar *measuresP = ParaInSimulation.get_var("measures");
	measuresP -> get(&in->measures, 1);
// 	NcVar *distP = ParaInSimulation.get_var("distances");
// 	distP -> get(&in->distances, 1);
	NcVar *CDBP = ParaInSimulation.get_var("CDB");
	CDBP -> get(&in->CDB,1);
	
	if (in->measures || in->CDB || in->synchrony)
	{
		//!read the spikes at which the measureable will be measured
		NcDim* szMeasureSpikesP =  ParaInSimulation.get_dim("szMeasureSpikes");
		int  szMeasureSpikes = szMeasureSpikesP->size();

		if (szMeasureSpikes)
		{
			in->measureSpikes.resize(szMeasureSpikes);
			in->measureTimes.clear();

			NcVar *measureSpikesP = ParaInSimulation.get_var("measureSpikes");
			measureSpikesP -> get(&in->measureSpikes.front(), szMeasureSpikes);

			cout << "The measurable will be measured and saved at the spiketimes" << endl;
		}
		else
		{
			//!otherwise read the specific times at which measureable will be measured
			NcDim* szMeasureTimesP =  ParaInSimulation.get_dim("szMeasureTimes");
			int  szMeasureTimes = szMeasureTimesP->size();

			in->measureSpikes.clear();
			if (szMeasureTimes)
			{
				in->measureTimes.resize(szMeasureTimes);

				NcVar *measureTimesP = ParaInSimulation.get_var("measureTimes");
				measureTimesP -> get(&in->measureTimes.front(), szMeasureTimes);

// 				cout << "The times at which the measurable will be measured and saved are: " << endl;
// 
// 				for (unsigned t=0; t<in->measureTimes.size(); t++)
// 					cout << in->measureTimes[t] << " ms" << "\t";
// 				cout << endl;
			}
			else
			{
				cout << "You need to specify at which times the phases are to be saved!" << endl;
				throw(1);
			}
		}
		
		if ((in->CDB > 0) and (in->CDB%2==0))
		{
			NcVar *D_decorrP = ParaInSimulation.get_var("D_decorr");
			D_decorrP->get(&in->D_decorr, 1);
		}
	}
	
	//!> if pop firing rate time bin size is nonzero, read in its value
	if (in->instPopRateBinSize > 0)                                             
	{
		NcVar *instPopRateBinSize= ParaInSimulation.get_var("instPopRateBinSize");
		instPopRateBinSize -> get(&in->instPopRateBinSize, 1);
	}
	
	
	
	//	read reference trajectory if given

	if (!fileRefTrajectory.empty())
	{
  		cout << "reading " << fileRefTrajectory << " ..." << endl;
		
		NcFile NcRefTraj(fileRefTrajectory.c_str(), NcFile::ReadOnly);
		
		NcDim *szrefTrajP = NcRefTraj.get_dim("measureNeuronsSz");
		int szrefTraj = szrefTrajP->size();
		
		vector<reell> refTrajTmp(szrefTraj);
		
		NcVar *refTrajP = NcRefTraj.get_var("measure1stStateVarNeurons");
		refTrajP -> get(&refTrajTmp.front(),szrefTraj);
		in->refTraj = vector<vector<double> > (szrefTraj/in->Nall, vector<double> (in->Nall));

		for (int n=0; n<szrefTraj/in->Nall; n++)
			in->refTraj[n].assign(refTrajTmp.begin() + n*in->Nall,refTrajTmp.begin() + (n+1)*in->Nall);

	}


	// The netCDF file is automatically closed by the NcFile destructor
	cout << "*** Parameter import successful!" << endl;
}

void LE_input_output::write_netcdf(string netcdf_filename, st_out* out, string fileNeurons, string fileTopology, string fileSimulation, bool applyPerturbation)
{
	//! Gather the final states from all processors on root
	int finalStateSz = out->finalStates.size();
	vector<reell> states;
	states.reserve(finalStateSz);				// This is the minimal size if all neurons have state dimension one.
	vector<reell> statesIdx(finalStateSz);

	if (finalStateSz)
	{
		int index = 0;
		for (int s=0; s<finalStateSz; s++)
		{
			statesIdx[s] = index;
			
			for (unsigned t=0; t<out->finalStates[s].size(); t++)
				states.push_back(out->finalStates[s][t]);

			index += out->finalStates[s].size();
		}

		gather_Nall(&states);
		gather_Nall(&statesIdx);

	}



	//! Gather the measurable of all neurons at each time step on root
	//!extract 1st and 2nd component of state

	int measureTimesSz = out->measure1stStateVarNeurons.size();
	if (measureTimesSz && out->measures)
	{
		for (int t=0; t<measureTimesSz; t++)
		{
			gather_Nall(&out->measure1stStateVarNeurons[t]);
			gather_Nall(&out->measure2ndStateVarNeurons[t]);
		}
	}





	if (myID == 0)
	{
		//! We are writing the results into an netCDF file.
		//! This is done in double precision, even though the simulation might have been done with long double, because netcdf can only store double.

		//! Create the file. The Replace parameter tells netCDF to overwrite
		//! this file, if it already exists.

		cout << "writing results to: "<< netcdf_filename << endl;
		NcFile dataFile(netcdf_filename.c_str(), NcFile::Replace);

		// You should always check whether a netCDF file creation or open
		// constructor succeeded.
		if (!dataFile.is_valid())
		{
		  cout << "Couldn't open file!\n";
		}



		//! Save the input files for which these results were obtained

// 		NcDim* NetDimP = dataFile.add_dim("inputfile_neurons_Sz", fileNeurons.size());
		NcDim* TopoDimP = dataFile.add_dim("inputfile_topology_Sz", fileTopology.size());
// 		NcDim* SimDimP = dataFile.add_dim("inputfile_simulation_Sz", fileSimulation.size());

// 		NcVar *fileNetP = dataFile.add_var("inputfile_neurons", ncChar, NetDimP);
		NcVar *fileTopoP = dataFile.add_var("inputfile_topology", ncChar, TopoDimP);
// 		NcVar *fileSimP = dataFile.add_var("inputfile_simulation", ncChar, SimDimP);

// 		const char* c = fileNeurons.c_str();
// 		fileNetP->put(c, fileNeurons.size());

		const char* c = fileTopology.c_str();
		fileTopoP->put(c, fileTopology.size());

// 		c = fileSimulation.c_str();
// 		fileSimP->put(c, fileSimulation.size());



		//! Write the data from the simulation to the file

		NcVar *SWP = dataFile.add_var("SW", ncDouble);
		NcVar *SCP = dataFile.add_var("SC", ncDouble);
		NcVar *TWP = dataFile.add_var("TW", ncDouble);
		NcVar *TCP = dataFile.add_var("TC", ncDouble);
		NcVar *rateWP = dataFile.add_var("rateW", ncDouble);
		NcVar *rateCP = dataFile.add_var("rateC", ncDouble);
		NcVar *IextScalingP = dataFile.add_var("IextScaling", ncDouble);
		
		double spikesPerNeuron = out->SW/(out->N*1000);
		SWP->put(&spikesPerNeuron);

		spikesPerNeuron = out->SC/(out->N*1000);
		SCP->put(&spikesPerNeuron);

		double time = out->TW/1000;
		TWP->put(&time);

		time = out->TC/1000;
		TCP->put(&time);

		double rate = out->rateW*1000;
		rateWP->put(&rate);

		rate = out->rateC*1000;
		rateCP->put(&rate);

		//! write scaling of Iext during ratefinding to netcdf file
// 		if in->
		  double IextScaling = out->IextScaling;
		  IextScalingP-> put(&IextScaling);


		// Write the results to the file

		int LyapunovExponentSz = out->LyapunovExponentsONS.size();
		int subLyapunovExponentSz = out->subLyapunovExponentsONS.size();
		int randomLyapunovExponentSz = out->randomLyapunovExponentsONS.size();
		
		unsigned spikeTrainSz = (out->calc_train) ? out->spikeTrainRef.size() : 1;
		
		int rateNeuronsSz = out->rateNeurons.size();
		int rateDistSz = out->rateDist.size();

		int cvNeuronsSz = out->cvNeurons.size();
		int cvDistSz = out->cvDist.size();

		int skewNeuronsSz = out->skewnessNeurons.size();
		int skewDistSz = out->skewnessDist.size();

		int kurtNeuronsSz = out->kurtosisNeurons.size();
		int kurtDistSz = out->kurtosisDist.size();

		const int PTSz = out->PT_steps;
		
		int measureTimesSz = out->measureTimes.size();
		int PopRateSz = out->instPopRate.size();


		if (finalStateSz)
		{
			//! Save the final states in a 1D vector
			NcDim* finalDimP = dataFile.add_dim("finalStatesSz", states.size());
			NcVar *finalP = dataFile.add_var("finalStates", ncDouble, finalDimP);
			vector<double> tmp=vector<double>(states.begin(),states.end());
			finalP->put(&tmp.front(), tmp.size());
			
			//! Save the indices at which places in the vector a new neuron starts
			NcDim* indexDimP = dataFile.add_dim("finalStatesIdxSz", statesIdx.size());
			NcVar *indexP = dataFile.add_var("finalStatesIdx", ncInt, indexDimP);
			vector<double> tmp2=vector<double>(statesIdx.begin(),statesIdx.end());
			indexP->put(&tmp2.front(), tmp2.size());
			
			//! Save the final currents
// 			cout << "Currents: " << endl;
// 			for (unsigned i=0;i<out->finalCurrents.size();i++)
// 				cout << out->finalCurrents[i] << endl;
			NcDim* curDimP = dataFile.add_dim("finalCurrentsSz", out->finalCurrents.size());
			NcVar *curP = dataFile.add_var("finalCurrents", ncDouble, curDimP);
			vector<double> tmp3=vector<double>(out->finalCurrents.begin(),out->finalCurrents.end());
			curP->put(&tmp3.front(), tmp3.size());

		}
		
		if (LyapunovExponentSz > 0)
		{
			//! save the Lyapunov exponents
			vector<double> LE(LyapunovExponentSz);
			for (int n=0; n<LyapunovExponentSz; n++)
				LE[n] = out->LyapunovExponentsONS[n]*1000;

			NcDim* LyaDimP = dataFile.add_dim("LEonsSz", LyapunovExponentSz);
			NcVar *LyapunovExponentsP = dataFile.add_var("LEons", ncDouble, LyaDimP);

			LyapunovExponentsP->put(&LE.front(), LyapunovExponentSz);
			
			if (randomLyapunovExponentSz>0)
			{

				//!save the statesequence randomized exponents
				vector<double> randomLE(randomLyapunovExponentSz*LyapunovExponentSz);

				for (int n=0; n<randomLyapunovExponentSz; n++)
					for (int m=0; m<LyapunovExponentSz; m++)
						randomLE[n*LyapunovExponentSz+m] = out->randomLyapunovExponentsONS[n][m]*1000;

				NcDim* randomLyaDimP = dataFile.add_dim("randomLEonsSz", randomLyapunovExponentSz*LyapunovExponentSz);
				NcVar *randomLyapunovExponentP = dataFile.add_var("randomLEons", ncDouble, randomLyaDimP);

				randomLyapunovExponentP->put(&randomLE.front(), randomLyapunovExponentSz*LyapunovExponentSz);
			}  

			//! Save the number of spikes during the warmup
			NcVar *SWONSP = dataFile.add_var("SWONS", ncDouble);
			double spikesPerNeuron = out->SWONS/out->N;
			SWONSP->put(&spikesPerNeuron);



			NcVar *pcP = dataFile.add_var("pLEONS", ncDouble);
			pcP->put(&out->pLEONS);



			//! save the convergence data
			int LEconvSz = out->LEconvergence.size();
			if (LEconvSz > 0)
			{
				vector<double> LE(LEconvSz*LyapunovExponentSz);
				for (int t=0; t<LEconvSz; t++)
					for (int n=0; n<LyapunovExponentSz; n++)
						LE[t*LyapunovExponentSz + n] = out->LEconvergence[t][n]*1000;

				NcDim* LEConvDimP = dataFile.add_dim("LEconvergenceSz", LEconvSz*LyapunovExponentSz);
				NcVar *LEConvP = dataFile.add_var("LEconvergence", ncDouble, LEConvDimP);

				LEConvP->put(&LE.front(), LEconvSz*LyapunovExponentSz);
			}

			int LEtimesSz = out->LEtimes.size();
			if (LEtimesSz > 0)
			{
				//! Save the times at which the orthonormalizations were done.
				cout << "saving times at which the orthonormalizations were done ..." << endl;
				NcDim* LEConvTimesDimP = dataFile.add_dim("LEtimesSz", LEtimesSz);
				NcVar *LEConvTimesP = dataFile.add_var("LEtimes", ncDouble, LEConvTimesDimP);
				vector<double> tmp=vector<double>(out->LEtimes.begin(),out->LEtimes.end());
				LEConvTimesP->put(&tmp.front(), LEtimesSz);
			}

			int localLESz = out->localLyapunovExponents.size();
			if (localLESz)
			{

				//! Save the local Lyapunov exponents from the CLV calcualtion.
				vector<double> LE(localLESz*LyapunovExponentSz);
				for (int t=0; t<localLESz; t++)
					for (int n=0; n<LyapunovExponentSz; n++)
						LE[t*LyapunovExponentSz + n] = out->localLyapunovExponents[t][n];

				NcDim* LElocDimP = dataFile.add_dim("localLESz", localLESz*LyapunovExponentSz);
				NcVar *LElocP = dataFile.add_var("localLE", ncDouble, LElocDimP);

				LElocP->put(&LE.front(), localLESz*LyapunovExponentSz);
			}
			
			int LEclvSz = out->LyapunovExponentsCLV.size();
			if (LEclvSz)
			{
				//! Save the asymptotic Lyapunov exponents calculated during the clv calculation.
				vector<double> LE(LEclvSz);
				for (int n=0; n<LEclvSz; n++)
					LE[n] = out->LyapunovExponentsCLV[n]*1000;

				NcDim* LyaDimP = dataFile.add_dim("LEclvSz", LyapunovExponentSz);
				NcVar *LyapunovExponentsP = dataFile.add_var("LEclv", ncDouble, LyaDimP);

				LyapunovExponentsP->put(&LE.front(), LyapunovExponentSz);

			}

		}
		
		if (subLyapunovExponentSz>0)
		{
			//! save the subnetwork's Lyapunov exponents
			vector<double> subLE(subLyapunovExponentSz);
			for (int n=0; n<subLyapunovExponentSz; n++)
				subLE[n] = out->subLyapunovExponentsONS[n]*1000;

			NcDim* subLyaDimP = dataFile.add_dim("subLEonsSz", subLyapunovExponentSz);
			NcVar *subLyapunovExponentsP = dataFile.add_var("subLEons", ncDouble, subLyaDimP);

			subLyapunovExponentsP->put(&subLE.front(), subLyapunovExponentSz);	
			
			//! save the sub netork's convergence data
			int subLEconvSz = out->subLEconvergence.size();
			if (subLEconvSz > 0)
			{
				vector<double> subLE(subLEconvSz*subLyapunovExponentSz);
				for (int t=0; t<subLEconvSz; t++)
					for (int n=0; n<subLyapunovExponentSz; n++)
						subLE[t*subLyapunovExponentSz + n] = out->subLEconvergence[t][n]*1000;

				NcDim* subLEConvDimP = dataFile.add_dim("subLEconvergenceSz", subLEconvSz*subLyapunovExponentSz);
				NcVar *subLEConvP = dataFile.add_var("subLEconvergence", ncDouble, subLEConvDimP);

				subLEConvP->put(&subLE.front(), subLEconvSz*subLyapunovExponentSz);

			}
			
			int sublocalLESz = out->sublocalLyapunovExponents.size();
			if (sublocalLESz)
			{

				//! Save the local Lyapunov exponents from the CLV calcualtion.
				vector<double> subLE(sublocalLESz*subLyapunovExponentSz);
				for (int t=0; t<sublocalLESz; t++)
					for (int n=0; n<subLyapunovExponentSz; n++)
						subLE[t*subLyapunovExponentSz + n] = out->sublocalLyapunovExponents[t][n];

				NcDim* subLElocDimP = dataFile.add_dim("sublocalLESz", sublocalLESz*subLyapunovExponentSz);
				NcVar *subLElocP = dataFile.add_var("sublocalLE", ncDouble, subLElocDimP);

				subLElocP->put(&subLE.front(), sublocalLESz*subLyapunovExponentSz);
			}
			
						
			int subLEclvSz = out->subLyapunovExponentsCLV.size();
			if (subLEclvSz)
			{
				//! Save the asymptotic Lyapunov exponents calculated during the clv calculation.
				vector<double> subLE(subLEclvSz);
				for (int n=0; n<subLEclvSz; n++)
					subLE[n] = out->subLyapunovExponentsCLV[n]*1000;

				NcDim* subLyaDimP = dataFile.add_dim("subLEclvSz", subLyapunovExponentSz);
				NcVar *subLyapunovExponentsP = dataFile.add_var("subLEclv", ncDouble, subLyaDimP);

				subLyapunovExponentsP->put(&subLE.front(), subLyapunovExponentSz);

			}
			
		}
		
		
		NcDim* PTDimP = dataFile.add_dim("PTSz", PTSz);
		
		cout << "train?" << endl;
		
		if (out->calc_train)
		{
			cout << "\t\tlet's see..." << endl;
			if (!applyPerturbation)
			{
				vector<double> trainTime(spikeTrainSz);
				vector<int> trainNeuron(spikeTrainSz);

				for (unsigned t=0; t<spikeTrainSz; t++)
				{
					trainTime[t] = out->spikeTrainRef[t].time/1000;
					trainNeuron[t] = out->spikeTrainRef[t].neuron;
				}
				
				NcDim* spTrDimP = dataFile.add_dim("spikeTrainSz", spikeTrainSz);
				NcVar *trainTimeP = dataFile.add_var("trainTime", ncDouble, spTrDimP);
				NcVar *trainNeuronP = dataFile.add_var("trainNeuron", ncInt, spTrDimP);

				trainTimeP->put(&trainTime.front(), spikeTrainSz);
				trainNeuronP->put(&trainNeuron.front(), spikeTrainSz);
			}
			
			if (applyPerturbation)
			{
				for (int p=0;p<PTSz;p++)
				{
					spikeTrainSz = (out->spikeTrainPert[p].size() > (unsigned) spikeTrainSz) ? out->spikeTrainPert[p].size() : spikeTrainSz;
				}
				
				vector<double> trainTime(spikeTrainSz);
				vector<int> trainNeuron(spikeTrainSz);
				
				for (unsigned t=0; t<spikeTrainSz; t++)
				{
					if (t<out->spikeTrainRef.size())
					{
						trainTime[t] = out->spikeTrainRef[t].time/1000;
						trainNeuron[t] = out->spikeTrainRef[t].neuron;
					}
					else
					{
						trainTime[t] = GSL_NAN;
						trainNeuron[t] = (int) GSL_NAN;
					}
				}
				NcDim* trainDimP = dataFile.add_dim("spikeTrainSz", spikeTrainSz);
				
				NcVar *trainTimeRefP = dataFile.add_var("trainTime", ncDouble, trainDimP);
				NcVar *trainNeuronRefP = dataFile.add_var("trainNeuron", ncDouble, trainDimP);
				
				trainTimeRefP->put(&trainTime.front(), spikeTrainSz);
				trainNeuronRefP->put(&trainNeuron.front(), spikeTrainSz);
				
				NcVar *trainTimePertP = dataFile.add_var("trainTimePert", ncDouble, PTDimP, trainDimP);
				NcVar *trainNeuronPertP = dataFile.add_var("trainNeuronPert", ncDouble, PTDimP, trainDimP);
				
				for (int p=0;p<PTSz;p++)
				{
					reell writeTrainTime[spikeTrainSz];
					reell writeTrainNeuron[spikeTrainSz];
// 					
					for (unsigned t=0;t<spikeTrainSz;t++)
					{
						if (t<out->spikeTrainPert[p].size())
						{
							writeTrainTime[t] = out->spikeTrainPert[p][t].time/1000;
							writeTrainNeuron[t] = out->spikeTrainPert[p][t].neuron;
						}
						else
						{
							writeTrainTime[t] = GSL_NAN;
							writeTrainNeuron[t] = GSL_NAN;
						}
						
					}
					trainTimePertP->put_rec(&writeTrainTime[0], p);
					trainNeuronPertP->put_rec(&writeTrainNeuron[0],p);
				}
			}
		}
		cout << "nope!" << endl;		
		//! write firing rates to netcdf file

		if (rateNeuronsSz > 0)
		{
			vector<double> rates(rateNeuronsSz);

			for (int i=0; i<rateNeuronsSz; i++)
				rates[i] = out->rateNeurons[i]*1000;

			NcDim* rateNDimP = dataFile.add_dim("rateNeuronsSz", rateNeuronsSz);
			NcVar *rateNeuronsP = dataFile.add_var("rateNeurons", ncDouble, rateNDimP);

			rateNeuronsP->put(&rates.front(), rateNeuronsSz);
		}

		if (rateDistSz > 0)
		{
			vector<double> x(rateDistSz);
			vector<double> y(rateDistSz);

			for (int n=0; n<rateDistSz; n++)
			{
			  x[n] = out->rateDist[n][0]*1000;
			  y[n] = out->rateDist[n][1]/1000;
			}

			NcDim* rateDDimP = dataFile.add_dim("rateDistSz", rateDistSz);
			NcVar *rateDistX = dataFile.add_var("rateDistX", ncDouble, rateDDimP);
			NcVar *rateDistY = dataFile.add_var("rateDistY", ncDouble, rateDDimP);

			rateDistX->put(&x.front(), rateDistSz);
			rateDistY->put(&y.front(), rateDistSz);
		}

		//! write coefficient of variation to netcdf file

		if (cvNeuronsSz > 0)
		{
			vector<double> cv(cvNeuronsSz);
			for (int i=0; i<cvNeuronsSz; i++)
				cv[i] = out->cvNeurons[i];

			NcDim* isiDimP = dataFile.add_dim("cvNeuronsSz", cvNeuronsSz);
			NcVar *isiP = dataFile.add_var("cvNeurons", ncDouble, isiDimP);

			isiP->put(&cv.front(), cvNeuronsSz);
		}

		if (cvDistSz > 0)
		{
			vector<double> x(cvDistSz);
			vector<double> y(cvDistSz);

			for (int n=0; n<cvDistSz; n++)
			{
			  x[n] = out->cvDist[n][0];
			  y[n] = out->cvDist[n][1];
			}

			NcDim* distDimP = dataFile.add_dim("cvDistSz", cvDistSz);
			NcVar *distX = dataFile.add_var("cvDistX", ncDouble, distDimP);
			NcVar *distY = dataFile.add_var("cvDistY", ncDouble, distDimP);

			distX->put(&x.front(), cvDistSz);
			distY->put(&y.front(), cvDistSz);
		}


		//! write skewness to netcdf file

		if (skewNeuronsSz > 0)
		{
			vector<double> skew(cvNeuronsSz);
			for (int i=0; i<skewNeuronsSz; i++)
				skew[i] = out->skewnessNeurons[i];


			NcDim* isiDimP = dataFile.add_dim("skewnessNeuronsSz", skewNeuronsSz);
			NcVar *isiP = dataFile.add_var("skewnessNeurons", ncDouble, isiDimP);

			isiP->put(&skew.front(), skewNeuronsSz);
		}

		if (skewDistSz > 0)
		{
			vector<double> x(skewDistSz);
			vector<double> y(skewDistSz);

			for (int n=0; n<skewDistSz; n++)
			{
				x[n] = out->skewnessDist[n][0];
				y[n] = out->skewnessDist[n][1];
			}

			NcDim* distDimP = dataFile.add_dim("skewnessDistSz", skewDistSz);
			NcVar *distX = dataFile.add_var("skewnessDistX", ncDouble, distDimP);
			NcVar *distY = dataFile.add_var("skewnessDistY", ncDouble, distDimP);

			distX->put(&x.front(), skewDistSz);
			distY->put(&y.front(), skewDistSz);
		}


		//! write kurtosis of variation to netcdf file

		if (kurtNeuronsSz > 0)
		{
			vector<double> kurt(cvNeuronsSz);
			for (int i=0; i<kurtNeuronsSz; i++)
				kurt[i] = out->kurtosisNeurons[i];


			NcDim* isiDimP = dataFile.add_dim("kurtosisNeuronsSz", kurtNeuronsSz);
			NcVar *isiP = dataFile.add_var("kurtosisNeurons", ncDouble, isiDimP);

			isiP->put(&kurt.front(), kurtNeuronsSz);
		}

		if (kurtDistSz > 0)
		{
			vector<double> x(kurtDistSz);
			vector<double> y(kurtDistSz);

			for (int n=0; n<kurtDistSz; n++)
			{
				x[n] = out->kurtosisDist[n][0];
				y[n] = out->kurtosisDist[n][1];
			}

			NcDim* distDimP = dataFile.add_dim("kurtosisDistSz", kurtDistSz);
			NcVar *distX = dataFile.add_var("kurtosisDistX", ncDouble, distDimP);
			NcVar *distY = dataFile.add_var("kurtosisDistY", ncDouble, distDimP);

			distX->put(&x.front(), kurtDistSz);
			distY->put(&y.front(), kurtDistSz);
		}
		
		if (out->synchrony)		// added by Alex (sync)
		{
			NcVar *chiP = dataFile.add_var("chi", ncDouble);
			chiP->put(&out->chi);
		}
		
		
		if (out->ISIdecorr_stats > 0)	// added by Alex (ISIdecorr)
		{
			int decorrSz = out->ISIdecorr.size();
			
			NcDim* decorrDimP = dataFile.add_dim("ISIdecorrSz", decorrSz);
			NcDim* permutP = dataFile.add_dim("permut", 2);
			
			NcVar *ISIdecorrP = dataFile.add_var("mean_ISIdecorr", ncDouble, decorrDimP, permutP);
			NcVar *cvISIdecorrP = dataFile.add_var("cv_ISIdecorr", ncDouble, decorrDimP, permutP);
			
			for (int rec=0;rec<decorrSz;rec++)
			{
				ISIdecorrP->put_rec(&out->ISIdecorr[rec][0], rec);
				cvISIdecorrP->put_rec(&out->cvISIdecorr[rec][0], rec);
			}
		}
		
		cout << "measures?" << endl;
		//!save the state variable and the times in double precision to the netcdf file
		if ((measureTimesSz && (out->measures || out->distances.size())))
		{
			cout << "\t\tlet's see..." << endl;
// 			int measureTimes_dim = out->CDB_state ? (out->CDB_idx+1) : measureTimesSz;
			int measureTimes_dim = measureTimesSz;
			NcDim* timesDimP = dataFile.add_dim("measureTimesSz", measureTimes_dim);
			vector<double> times(measureTimes_dim);
			for (int t=0; t<measureTimes_dim; t++)
				times[t] = out->measureTimes[t]/1000;
			NcVar *timesP = dataFile.add_var("measureTimes", ncDouble, timesDimP);
			timesP->put(&times.front(), measureTimes_dim);
			
			const int timesSz = out->measureTimes.size();
			const int neuronsSz = out->measureStates[0].size(); // make independent of array structure
			
// 			NcDim* timesDimP = dataFile.add_dim("timesSz", timesSz);
			NcDim* neuronsDimP = dataFile.add_dim("neuronsSz", neuronsSz);
			
			//! Save the measured measurables at the measuretimes.
			if (out->measures)
			{
				
				NcVar *statesP = dataFile.add_var("measureStates", ncDouble, timesDimP,neuronsDimP);
				
				cout << "timesDim: " << timesSz << ", neuronsDim: " << neuronsSz << endl;
				for (int rec=0;rec<timesSz;rec++)
					statesP->put_rec(&out->measureStates[rec][0], rec);
				
				cout << "\tpert" << endl;
				if (applyPerturbation)
				{
					cout << "\t\tlet's see..." << endl;
				
					NcVar *statesPertP = dataFile.add_var("measureStatesPert", ncDouble, PTDimP, timesDimP, neuronsDimP);
					
					for (int p=0;p<PTSz;p++)
					{
						reell writeDat[timesSz][neuronsSz];
						
						for (int t=0;t<timesSz;t++)
							for (int n=0;n<neuronsSz;n++)
								writeDat[t][n] = out->measureStatesPert[p][t][n];
						statesPertP->put_rec(&writeDat[0][0], p);
					}
				}
				cout << "\tnope!" << endl;
			}
			
// 			cout << "\t distances?" << endl;
			//! Save the distance at the measuretimes
			if (out->distances.size())
			{
				NcVar *distP = dataFile.add_var("distances",ncDouble,PTDimP,timesDimP);
				
				for (int p=0;p<PTSz;p++)
				{
					reell writeDat[timesSz];
					
					for (int t=0;t<timesSz;t++)
					{
						if (t >= out->distances.size())
							writeDat[t] = GSL_NAN;
						else
							writeDat[t] = out->distances[p][t];
					}
					distP->put_rec(&writeDat[0],p);
				}
			}
			cout << "\tnope!" << endl;
		}
		
// 		cout << "nope!" << endl;
		
// 		cout << "PT stuff" << endl;
// 		cout << out->CDB << endl;
		//! write results of comparing perturbed and reference trajectory
		if ((out->CDB > 0) && applyPerturbation)
		{
			if (out->CDB<=2)
			{
				NcVar *PTdistP = dataFile.add_var("PTdistances", ncDouble, PTDimP);
				PTdistP->put(&out->PTdistance.front(), out->PT_steps);
				
				NcVar *PTtimesP = dataFile.add_var("PTtimes", ncDouble, PTDimP);
				PTtimesP->put(&out->PTtime.front(), out->PT_steps);
			}
		}
// 		cout << "nope" << endl;
		
		//! write inst pop rate vector to file
		if (PopRateSz > 0) 
		{      
		  
			vector<double> PopRate(PopRateSz);
			for (int t=0; t<PopRateSz; t++)
				PopRate[t] = out->instPopRate[t]*1000;
			
			NcDim* instPopRateDim = dataFile.add_dim("instPopRateSz", PopRateSz);
			NcVar *instPopRate = dataFile.add_var("instPopRate", ncDouble, instPopRateDim);

			instPopRate->put(&PopRate.front(), PopRateSz);
		}

		// The netCDF file is automatically closed by the NcFile destructor
		cout << "*** Parameter export successful" << endl;
	}

}

void LE_input_output::display_in(st_in* in)
{

}

void LE_input_output::display_out(st_out* out)
{

}

void LE_input_output::gather_Nall(vector<reell>* vec)
{
#ifdef PAR
	//! Collects the data on root, which is stored locally in the vector \avec
	int Nloc = vec->size();

	//! Get local number of elements from all nodes
	vector<int> NlocV(nP);
	MPI_Gather(&Nloc, 1, MPI_INT, &NlocV.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	int Nall = 0;
	vector<int> indices(nP, 0);

	if (myID == 0)
	{
		for (int i=0; i<nP; i++)
		{
			indices[i] = Nall;
			Nall += NlocV[i];
		}
	}

	vector<double> vecAll(Nall);		//only needed on root, Nall = 0 on slaves

	MPI_Gatherv(&vec->front(), Nloc, MPI_DOUBLE, &vecAll.front(), &NlocV.front(), &indices.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//overwrite input vector on root
	if (myID == 0)
		*vec = vecAll;
#endif
}

void LE_input_output::gather_Nall(vector<int>* vec)
{
#ifdef PAR
	//! Collects the data on root, which is stored locally in the vector \avec
		int Nloc = vec->size();

		//! Get local number of elements from all nodes
		vector<int> NlocV(nP);
		MPI_Gather(&Nloc, 1, MPI_INT, &NlocV.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);

		int Nall = 0;
		vector<int> indices(nP, 0);

		if (myID == 0)
		{
			for (int i=0; i<nP; i++)
			{
				indices[i] = Nall;
				Nall += NlocV[i];
			}
		}

	vector<int> vecAll(Nall);		//only needed on root, Nall = 0 on slaves

	MPI_Gatherv(&vec->front(), Nloc, MPI_INT, &vecAll.front(), &NlocV.front(), &indices.front(), MPI_INT, 0, MPI_COMM_WORLD);

	//overwrite input vector on root
	if (myID == 0)
		*vec = vecAll;
#endif
}
