/*
 * LEanalysis.h
 *
 *  Created on: 04.01.2012
 *      Author: mik
 */

#ifndef LEANALYSIS_H_
#define LEANALYSIS_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include "reell.h"
#include "structures.h"
#include "LEnetwork.h"
#include "LEons.h"
// #include <gsl_math>
#include <gsl/gsl_randist.h>
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;

class LE_analysis
{
	public:
		LE_analysis(LE_network*);
		virtual ~LE_analysis();

		void setRate_warmup(st_in*, st_out*);

		void multitask(st_in*, st_out*, bool);							//!< all calculations at each iteration,
																		//!< this allows for an easy computation of, e.g., the Lyapunov exponents and the corresponding spike train
		void singletask(st_in*, st_out*);								//!< one calculation at a time in subsequent simulations

		inline void calc_randomLyapunovExponents(st_in*, st_out*);
		
	private:

		void preprocessing(st_in*, st_out*);
		void postprocessing(st_in*, st_out*);

		inline void prespike(st_in*, st_out*);
		inline void updatespike(st_in*, st_out*);

		inline void simple_iterations(st_in*, long long*, reell*);
		inline void simple_iterationsSubNet(long long*, long long*, reell*,st_in*);
		void bisection_externalCurrent(st_in*, st_out*);
		void bisection_externalCurrentSubNet(st_in*, st_out*);

		inline void calc_spikeTrain(st_in*, st_out*);

		void pre_ISIstatistics(st_in*, st_out*);
		inline void calc_ISIstatistics(st_in*, st_out*);
		void post_ISIstatistics(st_in*, st_out*);
		
		void pre_CDB(st_in*, st_out*);					// added by Alex (CDB)
		
		void pre_ISIdecorr_statistics(st_in*, st_out*);		// added by Alex (ISIdecorr)
		void calc_ISIdecorr_statistics(st_in*, st_out*);		// added by Alex (ISIdecorr)
		void post_ISIdecorr_statistics(st_in*, st_out*);		// added by Alex (ISIdecorr)

		void post_synchrony(st_in*, st_out*);				// added by Alex (sync)
		
		void pre_measurement(st_in*, st_out*);

		void warmup_ONS(st_out*);
		void pre_LyapunovExponents(st_in*, st_out*);
		inline void calc_LyapunovExponents(st_out*);
		void post_LyapunovExponents(st_in*, st_out*);

		void calc_LyapunovVectors(st_out*);

		int myID, nP;
		LE_network* net;												//!< the neural network

		bool calcJacobian;
		LE_ons *ons;													//!< pointer to the orthonormal system used for the calculation of the Lyapunov exponents
		LE_ons *sub_ons;													//!< pointer to the orthonormal system used for the calculation of the Lyapunov exponents

		vector<reell> lastSpike;										//!< save time of the last spike for the inter spike interval (isi) calculation
		vector<int> spikesOut;

		vector<vector<reell> > histogram_uniform(vector<reell>, int);

		vector<vector<reell> > stateOriginal;							//!< the original states for the perturbation calculations
		vector<reell> IextOriginal;										//!< the original currents for the perturbation calculations

		//! for perturbation calculations
		bool applyPerturbation;
		vector<vector<reell> > stateTmp;
		reell get_distOrth(st_in*, vector<reell> state1, vector<reell> state2);
		unsigned measureTimeCount, measureSpikeCount, curCount;

		//! for varying external currents calculations
		vector<reell> currentsOrg, currentsTmp;

		bool notConverged;
		reell LEmaxOld, LEminOld;
		
		vector<vector<vector<reell> > > stateTrajectory;									//! trajectory for random LE computation
	
};

#endif /* LEANALYSIS_H_ */
