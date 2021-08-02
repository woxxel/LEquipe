/*
 * structures.h
 *
 *  Created on: 31.01.2012
 *      Author: mik
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include "reell.h"
#include "gsl/gsl_rng.h"

using namespace std;

struct st_spike
{
	reell time; 					// time until next spike (time interval)
	int neuron;						// spiking neuron
};

struct st_synapse{
	int post;						// postsynaptic neuron
	double cplg;					// coupling strength
	double prob;					// release probability
};

struct st_double_int
{
	double d;
	int i;
};

struct st_twoDlinear
{
	double alpha;
	double beta;
	double gamma;
	double delta;
	double Cw;
	double tauS;
};

struct st_in
{
	//The netcdf file stores doubles, therefore doubles need to be read in, even when reell is set to long double

	//! Either the number of spikes \a SC or the simulation time \a TC must be given. If both are given, then time \a TC is chosen.
	long long SR, SW, SC;											//!< spikes in rate finding(R), warmup(W) and calculation(C)
	double TR, TW, TC;												//!< time in rate finding(R), warmup(W) and calculation(C)
	
	int pertDirections;
	
	double rateWnt, rateWntSubN, pR;												//!< wanted average firing rate in the network with precision pR

	int Nall; 														//!< number of neurons
	int Nloc;														//!< number of local neurons on one node
	int Ndrive;
	int Nrec;
	int Noffset;													//!< offset between local and global neuron index
	int subNall;											//!< array of (Nall,Nloc,and Noffset) for the subNetwork
	int subNloc;
	int subNoffset;

	int homogNet;													//!< homogeneous network or not

	vector<int> neuronType;											//!< neuron types in network
	vector<double> rapidness;										//!< AP onset rapidness of the neurons
	vector<double> poissonrate;										//!< poisson rate
	vector<double> gammashape;										//!< poisson rate
	vector<double> type1type2Para;									//!< parameter to interpolate btw. type 1 and type 2 neurons
	vector<double> reset;
	vector<double> threshold;
	vector<st_twoDlinear> twoDlinearParas;							//!< parameters for the twoDlinear model

	vector<double> tauM;											//!< membrane time constants of neurons
	vector<reell> Iext;											//!< initial external current to the neurons
	vector<vector<reell> > init;									//!< initial states of the neurons
	vector<vector<st_synapse> > synapses;							//!< postsynaptic connections of the neurons

	vector<int> train;												//!< calculate the spike train of these neurons
	bool calc_train;

	vector<int> ISIneurons;											//!< calculate the inter spike interval (ISI) statistics of these neurons
	int ISIstats;													//!< calculate this number of 'kinda' moments (1=mean, 2=cv, 3=skewness, 4=kurtosis)
	int ISIbins;													//!< number of bins for the ISI statistics distribution
	
	int ISIdecorr_stats;			// added by Alex (ISIdecorr)
	vector<int> ISIdecorr;
	
	int LyapunovExponents;											//!< calculate this number of Lyapunov exponents
	int subLyapunovExponents;										//!< calculate this number of subNetwork Lyapunov exponents
	int subLyapunovExponentsConvergence;								//!< Save the subLyapunov exponents at each orthonormalization step to a netcdf file.
	int randomizedLEspectra;												//!< Number of randomized state sequences for spectra calculation
	int long long SWONS;											//!< spikes of the warmup of the orthonormal system before the Lyapunov exponents calculation
	int seedONS;													//!< seed for random number generator when creating the initial orthonormal system
	int ONstep;														//!< step size of orthonormalizations
	int LyapunovExponentsConvergence;								//!< Save the Lyapunov exponents at each orthonormalization step to a netcdf file.
	double pLEONS;														//!< precision of Lyapunov exponents
	int CLV;														//!< true false to calculate the covariant Lyapunov vectors
	int long long SWCLV;											//!< spikes of the warmup of the covariant Lyapunov vectors

	int saveFinalState;							//!< save the final state of the network
	
	int addCur;														//!< Add an additional time varying current? 0=false, 1=true
	int addCurHomo;													//!< Is this current the same for all neurons? 0=flase, 1=true
	vector<int> addCurNeurons;										//!< These neurons will receive the additional current.
	vector<double> addCurTimes;										//!< vector with the times at which the additional current is applied
	vector<vector<double> > addCurIext;								//!< vector with the additional currents

	int measures;														//! export the states of all neurons?
// 	int distances;													//! export the distance between reference and perturbed trajectory (int value represents the norm)?
	vector<vector<reell> > refTraj;	// added by Alex (CDB)
	double D_decorr;			// added by Alex (CDB)
	vector<vector<reell> > puppetTrain;	// added by Alex (puppet)
	
	int pertSpike;													//!< skip one spike
	int pertSynapse;												//!< skip one synaptic transmission
	double pertSize;												//!< perturbation size of given perturbation
	long pertSeed;
	gsl_rng* rng;
	vector<vector<double> > pertVector;										//!< perturbation vector

	vector<double> measureTimes;										//!< times at which the phases are saved
	vector<int> measureSpikes;										//!< spikes at which the phases are saved

	double instPopRateBinSize;										//!< bin size with which to compute the instantaneous population firing rateC
	
	int synchrony;				// added by Alex (sync)
	int CDB;				// added by Alex (CDB)
	
};

struct st_out
{
	//Here reell should stay reell, because the variables might be used in the analysis.
	int N;															//!< number of neurons
	int subN;														//!<number of neurons in subnetwork
	int long long SW, SC;											//!< spikes in warmup and calculation
	reell TW, TC;													//!< time in warmup and calculation													//! biological time of the simulation

	reell rateW, rateC;												//!< network-averaged firing rate in warmup and simulation

	vector<reell> finalCurrents;									//!< final external currents of the neurons
	vector<vector<reell> > finalStates;								//!< final states of the neurons

	vector<st_spike> spikeTrain;									//!< spike time and neuron, defined in LEnetwork.h
	vector<st_spike> spikeTrainRef;
	vector<vector<st_spike> > spikeTrainPert;
	bool calc_train;
	
	int long long SWONS;											//!< spikes of the warmup of the orthonormal system before the Lyapunov exponents calculation
	int ONstep;														//!< step size of orthonormalizations
	vector<reell> LyapunovExponentsONS;								//!< Array of the asymptotic Lyapunov exponents.
	vector<reell> subLyapunovExponentsONS;								//!< Array of the asymptotic Lyapunov exponents.
	double pLEONS;													//!< precision of Lyapunov exponents
	vector<vector<reell> > LEconvergence;							//!< Save the Lyapunov exponents at each orthonormalization step,
	vector<vector<reell> > subLEconvergence;						//!< Save the subnetork's Lyapunov exponents at each orthonormalization step,
	vector<reell> LEtimes;											//!< and save the times at which the exponents are saved.

	int long long SWCLV;											//!< spikes of the warmup of the covariant Lyapunov vectors
	vector<reell> LyapunovExponentsCLV;								//!< Array of the asymptotic Lyapunov exponents.
	vector<reell> subLyapunovExponentsCLV;								//!< Array of the subnetwork's asymptotic Lyapunov exponents.

	vector<vector<reell> > localLyapunovExponents;							//!< Array of the local Lyapunov exponents at each orthonormalization (times stored in\a LEtimes).
	vector<vector<reell> > sublocalLyapunovExponents;							//!< Array of the subNetwork's local Lyapunov exponents at each orthonormalization (times stored in\a LEtimes).
	vector< vector<reell> > randomLyapunovExponentsONS;								//!< Array of the state-sequence randomized asymptotic Lyapunov exponents.
	int LyapunovExponentsConvergence;

	vector<reell> rateNeurons;										//!< the kinda moments of the inter spike interval distribution
	vector<reell> cvNeurons;
	vector<reell> skewnessNeurons;
	vector<reell> kurtosisNeurons;

	vector<vector<reell> > rateDist;								//!< and their distributions
	vector<vector<reell> > cvDist;
	vector<vector<reell> > skewnessDist;
	vector<vector<reell> > kurtosisDist;

	vector<vector<reell> > addCurRate;								//!< rates of individually chosen neurons

	int measures;														//! export the phases/states of all neurons?
// 	vector<reell> distances;										//! the distances between reference and perturbed trajectory

	vector<reell> measureTimes;										//!< times at which the phases are saved
	vector<vector<reell> > measure1stStateVarNeurons;							//!< measureable of the neurons at measureTimes
	vector<vector<reell> > measure2ndStateVarNeurons;							//!< measureable of the neurons at measureTimes
	vector<vector<reell> > measureStates;
	vector<vector<vector<reell> > > measureStatesPert;							//!< measureable of the neurons at measureTimes
	
	int CDB;				// added by Alex (CDB)
	bool CDB_state;				// added by Alex (CDB)
	int CDB_idx;				// added by Alex (CDB)
	reell CDB_dist;				// added by Alex (CDB)
	reell distInit;				// added by Alex (CDB)
	
	int PT_ct;
	int PT_steps;
	vector<double> PTdistance;		// added by Alex (perturbed trajectories)
	vector<vector<reell> > distances;
	vector<reell> PTtime;			// added by Alex (perturbed trajectories)
	
	int synchrony;				// added by Alex (sync)
	reell chi;				// added by Alex (sync)
	
	int ISIdecorr_stats;			// added by Alex (ISIdecorr)
	
	vector<vector<int> > count_ISIdecorr;	// added by Alex (ISIdecorr)
	vector<vector<reell> > ISIdecorr;	// added by Alex (ISIdecorr)
	vector<vector<reell> > cvISIdecorr;	// added by Alex (ISIdecorr)
	
	vector<reell> instPopRate;							//!<time series of instaneous population firing rate
	reell IextScaling;
};

#endif /* STRUCTURES_H_ */
