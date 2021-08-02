/*
 * LEanalysis.cpp
 *
 *  Created on: 04.01.2012
 *      Author: mik
 */

#include "LEanalysis.h"

LE_analysis::LE_analysis(LE_network* network) : net(network)
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

	ons = new LE_ons(1, 1, 1, 1);

	sub_ons = new LE_ons(1, 1, 1, 1);
	
// 	stateTmp.resize(net->Nrec);
}

LE_analysis::~LE_analysis()
{
	delete ons;
	delete sub_ons;
}



void LE_analysis::multitask(st_in* in, st_out* out, bool applyPerturbation)
{
	//! This method runs the simulation and analysis. The different results requested are all calculated from the same trajectory.
	//! An time-varying external current can be applied, although in this case, the Lyapunov spectrum cannot be calculated.
	if (myID == 0)
	{
		cout << endl << "on the " << ((!applyPerturbation) ? "reference" : "perturbed") << " trajectory, gonna calculate: " << endl;
		if (in->measures > 0)
		{
			string representation_str = (in->measures == 1) ? "voltage" : "phase";
			cout << "\t * the states in the " << representation_str << " representation at the provided:" << endl;
			
			if (in->measureTimes.size()) cout << "\t * times" << endl;
			if (in->measureSpikes.size()) cout << "\t * spikes" << endl;
		}
		if (in->pertSize > 0)
		{
			if (!applyPerturbation) cout << "\t * reference states" << endl;
			if (applyPerturbation) cout << "\t * distance to reference trajectory" << endl;
		}
		if (in->LyapunovExponents > 0)
		{
			cout << "\t * " << in->LyapunovExponents << " Lyapunov exponents" << endl;
			if (in->LyapunovExponentsConvergence > 0) cout << "\t * and its convergence" << endl;
			if (in->randomizedLEspectra > 0) cout << "\t * " << in->randomizedLEspectra << " copies of 'state-sequence randomized' spectra " << endl;
			if (in->CLV > 0) cout << "\t * " << in->LyapunovExponents << " covariant Lyapunov vectors" << endl;
		}
		if (in->train.size() > 0) cout << "\t * the spike train of " << in->train.size() << " neurons" << endl;
		if ((in->ISIneurons.size() > 0) && (in->ISIstats > 0)) cout << "\t * the inter spike interval statistics of " << in->ISIneurons.size() << " neurons, the first " << in->ISIstats << " kinda moments" << endl;
		if ((in->ISIdecorr.size() > 0) && (in->ISIdecorr_stats > 0)) cout << "\t * the inter spike interval statistics of connected neurons " << in->ISIdecorr.size() << " neurons" << endl;
		if (in->instPopRateBinSize > 0)  cout << "\t * the instantaneous population firing rate using " << in->instPopRateBinSize/1000 << "sec. bins" <<endl;
		if ((in->subNall > 0)&&(in->LyapunovExponents > 0))  cout << "\t * any Lyapunov stuff also for a " << in->subNall << "-neuron subnetwork" <<endl;
	}

	this->applyPerturbation = applyPerturbation;
	
	out->LyapunovExponentsConvergence = in->LyapunovExponentsConvergence;
	out->PT_steps = (applyPerturbation) ? in->pertDirections : 1;
	out->CDB = in->CDB;
	out->calc_train = in->calc_train;
	cout << " calc train: " << out->calc_train << endl;
	cout << " calc measures: " << out->measures << endl;
	cout << " CDB: " << out->CDB << endl;
	cout << "TC/SC: " << in->TC << ", " << in->SC << endl;
	
	if ((in->TC > 0) || (in->SC > 0))
	{
		for (out->PT_ct=0; out->PT_ct<out->PT_steps; out->PT_ct++)
		{
		  	out->TC = 0;
			out->SC = 0;
			
			out->CDB_state = false;
			
			//! run preprocessing first
			preprocessing(in, out);
			if (in->Ndrive > 0) net->reset_poisson(in);
			
			if (myID == 0) cout << endl << "simulation with SC = " << in->SC << " or TC = " << in->TC << " ms ..." << endl;

			time_t start, end, startAll;
			time(&start);
			time(&startAll);


			int CPUtime;
			int Tdone = 1, Sdone = 1;
			
			notConverged = false;
			
			while (((out->TC < in->TC) || (out->SC < in->SC) || notConverged) && !out->CDB_state)
			{
				prespike(in, out);
				
				if (!out->CDB_state)
					updatespike(in, out);

				//! Display how much is done already. If longer simulations are needed too reach the desired precision, it will exceed 100%.
				if ((in->TC > 0) && (out->TC >= in->TC/10.*Tdone))
					if (myID == 0)
					{
						time(&end);
						CPUtime = difftime(end, start);
						cout << "\t" << Tdone*10 << "% of TC done ... " << out->SC << " spikes @ t = " << out->TC/1000 << "s @ CPU time: " << CPUtime << "s. Remaining: "<< (10-Tdone)*CPUtime/Tdone << "s" << endl;
						Tdone++;

					}
				if ((in->SC > 0) && (out->SC >= in->SC/10.*Sdone))
					if (myID == 0)
					{
						time(&end);
						CPUtime = difftime(end, start);
						cout << "\t" << Sdone*10 << "% of TC done ... " << out->SC << " spikes @ t = " << out->TC/1000 << "s @ CPU time: " << CPUtime << "s. Remaining: "<< (10-Sdone)*CPUtime/Sdone << "s" << endl;
						Sdone++;
					}
			}
			
			out->rateC = out->SC/out->TC/in->Nrec;//net->get_N();
			if (myID == 0) cout << "\tspikes: " << out->SC << "\ttime: " << out->TC/1000 << "s\t -> avg. rate: " << out->rateC*1000 << " Hz" << endl;
			
			if (out->CDB == 0 && applyPerturbation)
			{
				out->PTdistance.push_back(GSL_NAN);
				out->PTtime.push_back(out->TC/1000);
			}
			
			cout << "measures: " << out->measures << endl;
			
			
			//! hand over measurements to out-array
			if ((out->measures > 0) || (out->CDB > 0))
			{
				if (applyPerturbation)
				{
// 					cout << "number of measurements: " << out->measure1stStateVarNeurons.size() << endl;
		
					out->measureStatesPert[out->PT_ct] = out->measure1stStateVarNeurons;
				}
				else
					out->measureStates = out->measure1stStateVarNeurons;
				
				out->measure1stStateVarNeurons.clear();
			}
			
			if (out->calc_train)
			{
				if (applyPerturbation)
				{
					out->spikeTrainPert.push_back(vector<st_spike>());
					for (unsigned t=0;t<out->spikeTrain.size();t++)
						out->spikeTrainPert[out->PT_ct].push_back(out->spikeTrain[t]);
				}
				else
					for (unsigned t=0;t<out->spikeTrain.size();t++)
						out->spikeTrainRef.push_back(out->spikeTrain[t]);
				out->spikeTrain.clear();
			}
		}
	}
	
	//!run postprocessing last!
	postprocessing(in, out);
}



void LE_analysis::preprocessing(st_in* in, st_out* out)
{

	if (myID == 0) cout << endl << "preprocessing ... " << endl << endl;

	//***************
	//! check  stuff
	//***************

	//! The analysis runs until the maximum of (spikes\a SC or time\a TC is reached)
	if ((in->TC <= 0) && (in->SC <= 0))
	{
		if (myID == 0) cout << "SC = " << in->SC << "\tTC = " << in->TC << endl;
		if (myID == 0) cout << "Neither the number of spikes SC or the simulation time TC were provided! -->> exit" << endl;
		throw(1);
	}

	//! Check, whether the Jacobian needs to be calculated, e.g. for the calculation of the Lyapunov exponents
	calcJacobian = false;
	if (in->LyapunovExponents > 0) //&& others
		calcJacobian = true;

	if (calcJacobian)
	{
		if (!net->get_allNeuronSame())
		{
			if (myID == 0) cout << "Cannot calculate the single spike Jacobian because the neurons in the network are of different types." << endl;
			throw(1);
		}

		if (in->addCur)
		{
			if (myID == 0) cout << "If the external current is changed within a spike interval, what does that mean for the Jacobian?" << endl;
			throw(1);
		}
	}

	//********************
	//! initialize  stuff
	//********************
	
	//! Pointer to the orhonormal system, e.g. for the calculation of the Lyapunov exponents.
	if (in->LyapunovExponents > 0)
		pre_LyapunovExponents(in, out);
	
	if (in->ISIstats > 0)
		pre_ISIstatistics(in, out);
	
	out->ISIdecorr_stats = in->ISIdecorr_stats;	// hand over information
	
	if (in->ISIdecorr_stats > 0)
		pre_ISIdecorr_statistics(in, out);
	
// 	if (in->CDB >0)
// 		pre_CDB(in, out);
	
	out->synchrony = in->synchrony;		// hand over information to save chi
	//! This only works with the phase models so far.
	if ( (in->measures > 0) || in->pertSize || in->pertSpike || in->pertSynapse || in->CDB || in->synchrony)
		pre_measurement(in, out);

	//! get the current external currents of the neurons that receive an additional current
	if (in->addCur)
	{
		curCount = 0;

		currentsOrg.resize(in->addCurNeurons.size());
		currentsTmp.resize(in->addCurNeurons.size());

		currentsOrg = net->get_externalCurrents_neurons(in->addCurNeurons);
	}

	//! initialize population firing rate times series vector
	if (in->instPopRateBinSize > 0)
		out->instPopRate = vector<reell> (int(ceil(in->TC/in->instPopRateBinSize)));

}



void LE_analysis::prespike(st_in* in, st_out* out)
{
	
	//! if a perturbation was or will be applied (e.g. for distance measurements), save the phase of the neurons at the given times
	//! if an additional current is provided, change the external currents
	if (in->measureTimes.size() || in->addCur || in->CDB)
	{
		
		//! The precision seems to be a problem, since there are a couple of difference calculations that influence the network evolution.
		reell timeToNextSpike = net->find_nextSpikes();
		reell timeToNextMeasurement = (in->measureTimes.size() && (measureTimeCount < in->measureTimes.size())) ? in->measureTimes[measureTimeCount] - out->TC : timeToNextSpike + 1;
		reell timeToNextCurrentStep = (in->addCur && (curCount < in->addCurTimes.size())) ? in->addCurTimes[curCount] - out->TC : timeToNextSpike + 1;
		
		while(((timeToNextSpike > timeToNextMeasurement) || (timeToNextSpike > timeToNextCurrentStep)) && !out->CDB_state)
		{
		
			//! if the next spike time is after the next time to measure the phase, then evolve the network to the next \aphaseTime.
			//! if the next spike time is after the next current change, then evolve to the next \aaddCurTime.
			
			//! initialize the vectors
			vector<vector<reell> > stateTmp;
					
			if (timeToNextMeasurement <= timeToNextCurrentStep)
			{
				if (timeToNextSpike > timeToNextMeasurement)
				{
					//! evolve the network until measureTime
					reell dt = timeToNextMeasurement;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;
					
					//PUT CASES HERE: get_state (phase or (V,I)) or get_phase

// 					if ( in->measures == 1 )
// 					{
// 						stateTmp = net->get_VRepState_Nloc();
// 						
// 						for (int n=0; n<net->get_Nloc(); n++)
// 						{
// 							out->measure1stStateVarNeurons[measureTimeCount][n]=stateTmp[n][0];
// 							out->measure2ndStateVarNeurons[measureTimeCount][n]=stateTmp[n][1];
// 						}
// 					}
					if ( (in->measures == 2) || (in->synchrony == 1) || ((in->pertSize > 0) && !applyPerturbation))		// need states for synchrony measure
					{
						stateTmp = net->get_PhRepState_Nrec();
						
						for (unsigned n=0; n<stateTmp.size(); n++)
						{
							out->measure1stStateVarNeurons[measureTimeCount][n]=stateTmp[n][0];
							out->measure2ndStateVarNeurons[measureTimeCount][n]=stateTmp[n][1];
						}
					}
					
// 					cout << in->CDB <<", " << applyPerturbation << endl;
					if ((in->CDB > 0) && applyPerturbation)
					{
						vector<reell> stateTest (in->Nrec);
						stateTmp = net->get_PhRepState_Nrec();
						
						for (unsigned n=0; n<stateTmp.size(); n++)
							stateTest[n] = stateTmp[n][0];
						
						reell distOrth = get_distOrth(in,stateTest,out->measureStates[measureTimeCount]);
// 						cout << "distance: " << distOrth << endl;
						
						if (in->CDB>2)
						{
							out->distances[out->PT_ct].push_back(distOrth);
// 							cout << "save distance: " << out->distances[out->PT_ct][measureTimeCount] << endl;
						}
						
						// break if converged, or at last timestep
						if (((distOrth <= pow(10,-8)) || (measureTimeCount >= in->measureTimes.size()-1)) || ((in->CDB%2 == 0) and (distOrth >= in->D_decorr*0.7)))
						{
							cout << "Break at: " << distOrth << " at time t = " << out->TC << endl;
							out->CDB_state = true;
							out->CDB_idx = measureTimeCount;
							
							if (in->CDB<=2)
							{
								// save the results of this simulation
								out->PTdistance.push_back(distOrth);
								out->PTtime.push_back(out->TC/1000);
							}
						}
					}

					if ((myID == 0) and !applyPerturbation) out->measureTimes[measureTimeCount] = out->TC;
					
					//! go to the next time step
					measureTimeCount++;
				}

				if ((timeToNextSpike > timeToNextCurrentStep) && !out->CDB_state)
				{
					//! evolve the neurons to the time at which the current is changed
					reell dt = timeToNextCurrentStep - timeToNextMeasurement;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;

					//! change the current of the specified neurons
					for (unsigned n=0; n<in->addCurNeurons.size(); n++)
					{
						int nHomo = (in->addCurHomo) ? 0 : n;
						currentsTmp[n] = currentsOrg[n] + in->addCurIext[curCount][nHomo];
				
					}
					net->set_externalCurrents_neurons(in->addCurNeurons, currentsTmp);

					//! go to the next time step
					curCount++;
				}
			}
			else	//the order is the other way around
			{
				if (timeToNextSpike > timeToNextCurrentStep)
				{
					//! evolve the neurons to the time at which the current is changed
					reell dt = timeToNextCurrentStep;
					net->evolve_dt(dt);
					out->TC += dt;

					//! change the current of the specified neurons
					for (unsigned n=0; n<in->addCurNeurons.size(); n++)
					{
						int nHomo = (in->addCurHomo) ? 0 : n;
						currentsTmp[n] = currentsOrg[n] + in->addCurIext[curCount][nHomo];
					}
					net->set_externalCurrents_neurons(in->addCurNeurons, currentsTmp);
					timeToNextSpike = net->find_nextSpikes();
					//! go to the next time step
					curCount++;
				}

				if (timeToNextSpike > timeToNextMeasurement)
				{
					//! evolve the network until measureTime
					reell dt = timeToNextMeasurement - timeToNextCurrentStep;
					net->evolve_dt(dt);
					out->TC += dt;
					timeToNextSpike -= dt;
					
// 					if ( in->measures == 1 )
// 					{
// 						stateTmp = net->get_VRepState_Nloc();
// 						for (int n=0; n<net->get_Nloc(); n++)
// 						{
// 							out->measure1stStateVarNeurons[measureTimeCount][n]=stateTmp[n][0];
// 							out->measure2ndStateVarNeurons[measureTimeCount][n]=stateTmp[n][1];
// 						}
// 					}
					if ( in->measures == 2 )
					{
						stateTmp = net->get_PhRepState_Nrec();
						for (unsigned n=0; n<stateTmp.size(); n++)
						{
							out->measure1stStateVarNeurons[measureTimeCount][n]=stateTmp[n][0];
							out->measure2ndStateVarNeurons[measureTimeCount][n]=stateTmp[n][1];
						}
					}
					if ((myID == 0) and !applyPerturbation) out->measureTimes[measureTimeCount] = out->TC;

					//! go to the next time step
					measureTimeCount++;
				}
			}
			timeToNextMeasurement = (in->measureTimes.size() && (measureTimeCount < in->measureTimes.size())) ? in->measureTimes[measureTimeCount] - out->TC : timeToNextSpike + 1;
			timeToNextCurrentStep = (in->addCur && (curCount < in->addCurTimes.size())) ? in->addCurTimes[curCount] - out->TC : timeToNextSpike + 1;
		}
	}

	if (!out->CDB_state)
	{
		
		//! The normal evolution of all neurons without any external changes
		reell timeToNextSpike = net->find_nextSpikes();

		net->evolve_dt(timeToNextSpike);

		out->TC += timeToNextSpike;
	}		
	if(in->randomizedLEspectra)
		stateTrajectory.push_back(net->get_simState_Nloc());

}


reell LE_analysis::get_distOrth(st_in* in, vector<reell> state1, vector<reell> state2)
{
	reell projTmp = 0;
	reell distOrth = 0;
	
	reell norm = in->Nrec;
	
	vector<reell> distTmp(in->Nrec);
	
	for (int n=0; n<in->Nrec; n++)
	{
		distTmp[n] = state1[n] - state2[n];
		projTmp += distTmp[n]/norm;
	}
	
	for (int n=0; n<in->Nrec; n++)
		distOrth += gsl_pow_2(distTmp[n] - projTmp);

	return sqrt(distOrth);
}


void LE_analysis::updatespike(st_in* in, st_out* out)
{
	
	//! The normal update of all neurons at spike reception.
	net->evolve_spike(calcJacobian);
// 	cout << in->Nrec << endl;
	vector<st_spike> spikeTmp = net->get_spikes();
	if (spikeTmp[0].neuron < in->Nrec)
		out->SC++;

	if (in->calc_train)
		calc_spikeTrain(in, out);

	if (in->LyapunovExponents > 0)
		calc_LyapunovExponents(out);

	if (in->ISIstats > 0)
		calc_ISIstatistics(in, out);
	
	if (in->ISIdecorr_stats > 0)
		calc_ISIdecorr_statistics(in, out);
	
// 	vector<vector<reell> > stateTmp;
	
	//! Save the states of the neurons at the provided spikes in the network.
	if (in->measureSpikes.size())
		if ((measureSpikeCount < in->measureSpikes.size()) && (out->SC == in->measureSpikes[measureSpikeCount]))
		{
// 			if ( in->measures == 1 )
// 			{
// 				stateTmp = net->get_VRepState_Nloc();
// 				for (int n=0; n<in->Nrec; n++)
// 				{
// 					out->measure1stStateVarNeurons[measureSpikeCount][n]=stateTmp[n][0];
// 					out->measure2ndStateVarNeurons[measureSpikeCount][n]=stateTmp[n][1];
// 				}
// 			}
			if ( in->measures == 2 )
			{
				stateTmp = net->get_PhRepState_Nrec();
				for (unsigned n=0; n<stateTmp.size(); n++)
				{
					out->measure1stStateVarNeurons[measureSpikeCount][n]=stateTmp[n][0];
					out->measure2ndStateVarNeurons[measureSpikeCount][n]=stateTmp[n][1];
				}
			}
			if ((myID == 0) and !applyPerturbation) out->measureTimes[measureSpikeCount] = out->TC;

			//! go to the next provided spike number
			measureSpikeCount++;
		}
// 	stateTmp.clear();
	
	//! Add the spike to the appropriate bin
	if ( (out->TC <= in->TC ) && ( in->instPopRateBinSize > 0 ) )
		out->instPopRate[int(out->TC / in->instPopRateBinSize)] += 1;
	
}



void LE_analysis::postprocessing(st_in* in, st_out* out)
{
	//! Save the final state of the network.
	if (in->saveFinalState)
	{
		cout << "saving the final state and external currents of the network ..." << endl;

		//! Get the currents of all N neurons.
		out->finalCurrents = net->get_externalCurrents_Nloc();

		//! Get the states of the neurons.
		out->finalStates = net->get_simState_Nloc();
	}
  
	if (myID == 0)
	{
		cout << endl << "postprocessing ... " << endl;

		if (in->LyapunovExponents)
		{
			post_LyapunovExponents(in, out);

			if (in->CLV)
				calc_LyapunovVectors(out);
		}

		if (in->ISIstats > 0)
			post_ISIstatistics(in, out);

		if (in->synchrony)
			post_synchrony(in, out);
		
		if (in->ISIdecorr_stats > 0)
			post_ISIdecorr_statistics(in, out);
		
		//! normalize the instantaneous population spike rate by bin size and neuron number
		if (in->instPopRateBinSize > 0)
			for(unsigned n=0; n<out->instPopRate.size(); n++)
				out->instPopRate[n] /= in->instPopRateBinSize*net->get_N();
	}
	if ((in->randomizedLEspectra) && (in->LyapunovExponents))
	      calc_randomLyapunovExponents(in, out);
		  
}



void LE_analysis::setRate_warmup(st_in* in, st_out* out)
{
	//! setRate_warmup finds the wanted rate and does a warmup of the network.

	//! if wanted and possible, find the right external currents to yield a specific average firing rate\a rateWnt
	out->IextScaling = 1; //! Outputs Scaling of initial external current.
	if ( ((in->rateWntSubN > 0) || (in->rateWnt > 0)) && ((in->TR>0) || (in->SR>0)) )
	{
		if (myID == 0) cout << "adapt the external currents to yield the wanted average firing rate:"<< endl;
		cout << "Last Neurons Type ";
		int lastIdx = (in->homogNet) ? 0 : (in->Nall - 1);
		cout << in->neuronType[lastIdx] << endl;
		
		if ((in->subNall>0) && ((in->neuronType[lastIdx]==5) || (in->neuronType[lastIdx]==51) || (in->neuronType[lastIdx]==7)))
			bisection_externalCurrentSubNet(in, out);

		//  bisection_externalCurrentOnlySubNet(in, out);
		else
			if (in->subNall>0)
				bisection_externalCurrentSubNet(in, out);
			else
				bisection_externalCurrent(in, out);
	}
	else
		if (myID == 0) cout << "The provided external currents Iext will be used for the simulation." << endl;

	if ((in->TW > 0) || (in->SW > 0))
	{
		//! run the warmup
		if (myID == 0) cout << endl << "warmup with SW = " << in->SW << " or TW = " << in->TW << " ms ..." << endl;

		out->TW = in->TW;
		out->SW = in->SW;
		simple_iterations(in, &out->SW, &out->TW);
	
		out->rateW = out->SW/out->TW/in->Nrec;//net->get_N();
		if (myID == 0) cout << "\tspikes: " << out->SW << "\ttime: " << out->TW/1000 << "s\t -> avg. rate: " << out->rateW*1000 << " Hz" << endl;
	}
	else
	  	if (myID == 0) cout << "Network will be run without a warmup." << endl;


}

void LE_analysis::simple_iterations(st_in* in, long long* spikes, reell* time)
{

	//! simple_iterations simply iterates the network from the current state for at least\a spikes spikes and \atime time.
	//! The given arguments\aspikes and\a time will be overwritten by the actual number of spikes and time of these iterations.

	long long spikesMin = *spikes;
	reell timeMin = *time;

	*spikes = 0;
	*time = 0;

	while (!((*spikes >= spikesMin) && (*time >= timeMin)))
	{
		//! run one iteration of the network
		reell dt = net->find_nextSpikes();
		if (isnan(dt))
		{
			cout << "no spike time found - set rate found to 0" << endl;
			break;
		}
		if (*time + dt > timeMin)	// only warm up until specified points
		{
			dt = timeMin - *time;
			net->evolve_dt(dt);
		}
		else
		{
			net->evolve_dt(dt);
			net->evolve_spike(false);
			vector<st_spike> spikeTmp = net->get_spikes();
			if (spikeTmp[0].neuron < in->Nrec)
				*spikes += 1;
		}
			
		*time += dt;
	}
}

void LE_analysis::bisection_externalCurrent(st_in* in, st_out* out)
{
	if (myID == 0) cout << "rateWnt = " << in->rateWnt*1000 << " Hz with precision pR = " << in->pR << " and SR = " << in->SR << " or TR = " << in->TR << " ms ..." << endl;
	
	//! find the external currents to match the wanted firing rate by bisection, the maximal number of iterations is 30
	int iterMax = 30;			//iterMax is OK to be hard coded I think
	vector<reell> Itmp = in->Iext;

	reell rateTmp = 0, rateUp = 0, rateDown = 0;
	reell factorTmp = 1, factorUp = 1, factorDown = 0;
	reell factorTmpTmp=0,rateTmpTmp=0;
	
	//! run the bisection with linear guess, since the firing rate is linear to the external current in the balanced state
	bool calc = true;
	int iter = 0;

	while ((calc) && (iter < iterMax))
	{
		if (in->homogNet)
			in->Iext[0] = factorTmp*Itmp[0];
		else
			for (int n=0; n<net->get_Nloc(); n++)
				in->Iext[n] = factorTmp*Itmp[n];

		net->set_externalCurrents_Nloc(in);
		net->set_simState_Nloc(in->init);				//!The reset to the initial conditions might not be too important here and could be dropped for a speed up in large networks.

		long long spikesTmp = in->SR;
		reell timeTmp = in->TR;

		simple_iterations(in, &spikesTmp, &timeTmp);

		if (spikesTmp == 0)
			rateTmp = 0;
		else
			rateTmp = spikesTmp/timeTmp/net->get_N();
		
		if (myID == 0) cout << "\t" << factorTmp << "*Iext yielded f = " << rateTmp*1000 << " Hz" << endl;

		calc = (abs(1 - rateTmp/in->rateWnt) > in->pR);

		if (calc)		//works going down to negative currents and up from positive currents
		{
			if (rateTmp > in->rateWnt)
			{
				factorUp = factorTmp;
				rateUp = rateTmp;
				if (iter == 0)
				{
					factorTmpTmp=factorUp;
					rateTmpTmp=rateUp;
				}
			}
			else
			{
				factorDown = factorTmp;
				rateDown = rateTmp;
			}
			
			if ( (iter==1) && !(factorTmpTmp==0) && (rateTmp > in->rateWnt))
				factorDown = -rateUp*(factorUp - factorTmpTmp)/(rateUp - rateTmpTmp) + factorUp;   //!adjusts lower bound of current to 0-rate value (otherwise it is assumed to be 0)
				
			//! guess the next factor for the next current (assuming a linear slope)
			factorTmp = factorDown + (in->rateWnt - rateDown)/(rateUp - rateDown)*(factorUp - factorDown);

			//! At the beginning the current rate can be below the wanted rate, then the following is necessary to find the initial upper bound for the bisection.
			if (rateUp < in->rateWnt)
			{
				factorUp *= 2; //to put it above desired current (factorUp here can be negative)
				factorTmp = factorUp;
			}
			
		}
		iter++;
	}

	if (calc)
		if (myID == 0)
			{
			cout << "the external current couldn't be found to yield the desired firing rate rateWnt = ";
			cout << in->rateWnt*1000 << " Hz with precision pR = " << in->pR << endl;
			}

	if (myID == 0) cout << "the external currents are set to " << factorTmp << "*Iext" << endl;
}

void LE_analysis::bisection_externalCurrentSubNet(st_in* in, st_out* out)
{
	//! find the external currents for each subnetwork to match the wanted firing rates for each by bisection, the maximal number of iterations is 30
	int iterMax = 100;			//iterMax is OK to be hard coded I think

	reell rateTmp, rateUp, rateDown;
	reell factorTmp, factorUp, factorDown;
	bool calc;
	int iter;
	vector<reell> Itmp;
	
	//!loop parameters over the super and sub network
	reell rateWntvec[2] =   {in->rateWnt	,in->rateWntSubN};
	string networkname[2] = {"SuperNet"	,"SubNet"};
	int lowind[2]	=	{in->subNloc	,0};
	int highind[2] = 	{net->get_Nloc(),in->subNloc};
	int lowindGLOBAL[2] =	{in->subNall	,0	};
	int highindGLOBAL[2] =	{net->get_N()	,in->subNall};

	long long spikesHere;
	long long spikesGlobal;
	long long spikeLimit;
	
	reell time;
	
	//!!Fixing the complement of subNetwork, then subnetwork, i.e. only works networks for no connections from sub to the rest. 
	
	for (int l=0; l<2; l++) 
	{
		if (in->rateWnt ==0)
			l++;
		cout << "No wanted rate for SuperNet defined."<<endl;

		if (myID == 0) cout << "rateWnt for "<< networkname[l] <<"= " << rateWntvec[l]*1000 << " Hz with precision pR = " << in->pR << " and SR = " << in->SR << " or TR = " << in->TR << " ms ..." << endl;
		
		Itmp = in->Iext;  //store all external currents 
				
		rateTmp = 0, rateUp = 0, rateDown = 0;
		factorTmp = 1, factorUp = 1, factorDown = 0;


		//! run the bisection with linear guess, since the firing rate is linear to the external current in the balanced state
		calc = true;
		iter = 0;

		while ((calc) && (iter < iterMax))
		{

			for (int n=lowind[l]; n<highind[l]; n++) 
				in->Iext[n] = factorTmp*Itmp[n];  //set all relevant external currents
			
			net->set_externalCurrents_Nloc(in);	//set local external currents
			net->set_simState_Nloc(in->init);				//!The reset to the initial conditions might not be too important here and could be dropped for a speed up in large networks.
			//simple_iterationsSubNet(&spikesTmpSubN,&spikesTmpSupN, &timeTmp, in);
			spikesHere = 0;
			spikesGlobal = 0;
			spikeLimit = 2147483647-1;//! sets 2^31-1 as limit to the number of spikes during rate finding. Otherwise the simulation might be stuck, because of a silent subnetwork. 
			time = 0;
			while (((spikesHere < in->SR) || (time < in->TR)) && (spikesGlobal < spikeLimit))
				//while (!( (spikes >= in->SR) && (time >= in->TR)) )
			{ 
				//! run one iteration of the network
				reell dt = net->find_nextSpikes();
// 				cout << "iterate: " << dt << ",  ";
				net->evolve_dt(dt);
				net->evolve_spike(false);
				spikesGlobal +=1;
				time += dt;
				vector<st_spike> spikeTmp = net->get_spikes();
				for(unsigned s=0; s<spikeTmp.size(); s++)
				{
					if ( (spikeTmp[s].neuron >= lowindGLOBAL[l]) & (spikeTmp[s].neuron< highindGLOBAL[l] ) )
						spikesHere += 1;
				}
			}
		 //   cout << "spikesGlobal "<< spikesGlobal<<endl;

		   // if ((spikesGlobal == spikeLimit) && (in->TR==0)) cout << "Stopped after " << spikesGlobal<< " spikes. Please use TR for rate finding in subnetwork, as the subnetwork might be silent!"<< endl;
			if ((spikesGlobal == spikeLimit)) cout << "Stopped after " << spikesGlobal<< " spikes. Please use TR for rate finding in subnetwork, as the subnetwork might be silent!"<< endl;

			rateTmp = spikesHere/time/(highindGLOBAL[l]-lowindGLOBAL[l]);

			if (myID == 0) cout << "\t" << factorTmp << "*Iext in "<< networkname[l]<< " yielded f = " << rateTmp*1000 << " Hz" << endl;

			calc = (abs(1 - rateTmp/rateWntvec[l]) > in->pR);

			if (calc)
			{
				if (rateTmp > rateWntvec[l])
				{
					factorUp = factorTmp;
					rateUp = rateTmp;
				}
				else
				{
					factorDown = factorTmp;
					rateDown = rateTmp;
				}

				//! guess the next factor for the next current (assuming a linear slope)
				factorTmp = factorDown + (rateWntvec[l] - rateDown)/(rateUp - rateDown)*(factorUp - factorDown);

				//! At the beginning the current rate can be below the wanted rate, then the following is necessary to find the initial upper bound for the bisection.
				if (rateUp < rateWntvec[l])
				{
					factorUp *= 2;
					factorTmp = factorUp;
				}

			}
			iter++;
		}

		if (calc)
			if (myID == 0)
				{
				cout << "the external current for the" << networkname[l] << " couldn't be found to yield the desired firing rate rateWnt = ";
				cout << rateWntvec[l]*1000 << " Hz with precision pR = " << in->pR << endl;
			}
out->IextScaling = factorTmp; 
		if (myID == 0) cout << "the external currents for " << networkname[l] << " are set to " << factorTmp << "*Iext" << endl << endl;
	}
}


void LE_analysis::calc_spikeTrain(st_in* in, st_out* out)
{
	if (myID == 0)
	{
		vector<st_spike> spikeTmp = net->get_spikes();

		//! Check if the spiking neurons are part of the wanted neurons for the spike train.
		//! Don't know if this could be done more efficiently than going through the train vector all the time.

		for(unsigned s=0; s<spikeTmp.size(); s++)
			for(vector<int>::iterator t=in->train.begin(); t<in->train.end(); t++)
			{
				if (spikeTmp[s].neuron == *t)
				{
					//!replace the time in spikeTmp, which was length of this interval, with the actual time\a time
					spikeTmp[s].time = out->TC;

					//!< add the currently spiking neurons\a spikeTmp to the output spike vector\a out.spikes
					out->spikeTrain.push_back(spikeTmp[s]);

					//! don't need to check for another occurence of the same neuron in\a in->train
					break;
				}
			}
	}
}


// void LE_analysis::pre_CDB(st_in* in, st_out* out)
// {	
// 	vector<vector<reell> > stateTmp = net->get_PhRepState_Nloc();
// 	
// 	reell distTmp = 0;
// 	
// 	for (int n=0;n<in->Nrec;n++)
// 		distTmp += gsl_pow_2(in->refTraj[0][n] - stateTmp[n][0]);
// 	
// 	out->distInit = sqrt(distTmp);
// }

void LE_analysis::pre_ISIstatistics(st_in* in, st_out* out)
{
	if (!(in->ISIneurons.size() > 0))
		in->ISIstats = 0;

	if (myID == 0)
	{
		//! initialize the vectors
		switch (in->ISIstats)
		{
			case 4:
				out->kurtosisNeurons = vector<reell> (in->ISIneurons.size());
				// no break, since all lower moments need to be calculated as well

			case 3:
				out->skewnessNeurons = vector<reell> (in->ISIneurons.size());
				// no break

			case 2:
				out->cvNeurons = vector<reell> (in->ISIneurons.size());
				// no break

			default:
				out->rateNeurons = vector<reell> (in->ISIneurons.size());
				spikesOut = vector<int> (in->ISIneurons.size());
				lastSpike = vector<reell> (in->ISIneurons.size());
				// no break
		}
	}
}

void LE_analysis::calc_ISIstatistics(st_in* in, st_out* out)
{
	if (myID == 0)
	{
		//! get the spike vector
		vector<st_spike> spikes = net->get_spikes();

		//! calculate the moments of the inter spike interval distribution
		for(vector<st_spike>::iterator s = spikes.begin(); s < spikes.end(); s++)
		{
			//! Check whether the spiking neuron is one of the considered ones.
			unsigned n = 0;
			while(n < in->ISIneurons.size())
			{
				if ((*s).neuron == in->ISIneurons[n])
					break;
				else
					n++;
			}

			// If the neurons wasn't found, it's not considered. Continue with the next spike.
			if (n == in->ISIneurons.size())
				continue;

			//! count the spikes to calculate the neuron's rates
			spikesOut[n]++;

			//! save the time of the first spike in the rate vector
			if (spikesOut[n] == 1)
				out->rateNeurons[n] = out->TC;

			//! calculate the isi moments
			if (in->ISIstats > 0)
			{
				if (spikesOut[n] > 0)			// if there was one spike before, we can calculate the isi
				{
					reell isi = out->TC - lastSpike[n];
					reell isi2 = isi*isi;

					switch (in->ISIstats)
					{
						case 4:
							out->kurtosisNeurons[n] += isi2*isi2;
							// no break, since all lower moments need to be calculated as well

						case 3:
							out->skewnessNeurons[n] += isi2*isi;
							// no break
						case 2:
							out->cvNeurons[n] += isi2;
							// no break
					}
				}

				lastSpike[n] = out->TC;
			}
		}
	}
}

void LE_analysis::post_ISIstatistics(st_in* in, st_out* out)
{
	for (unsigned n=0; n<in->ISIneurons.size(); n++)
	{
		//! Calculate the statistics only if there were more than 6 spikes.
		if (spikesOut[n] > 6)
		{
			//! Calculate the rates.

			// out->rateNeurons[n] currently holds the time of the first spike
			// lastSpike[n] currently holds the time of the last spike
			out->rateNeurons[n] = (spikesOut[n] - 1)/(lastSpike[n] - out->rateNeurons[n]);

			// Calculate the moments by dividing with the number of spikes of this neuron which where considered (the first spike wasn't considered => -1).
			switch (in->ISIstats)
			{
				case 4 :
					out->kurtosisNeurons[n] /= spikesOut[n] - 1;
					// no break, since all lower moments need to be calculated as well

				case 3 :
					out->skewnessNeurons[n] /= spikesOut[n] - 1;
					//no break

				case 2 :
					out->cvNeurons[n] /= spikesOut[n] - 1;
					//no break
			}

			switch (in->ISIstats)
			{
				case 4 :
					//! Calculate the 4th standardized moment (kurtosis)
					out->kurtosisNeurons[n] = out->kurtosisNeurons[n] -
											4*out->skewnessNeurons[n]/out->rateNeurons[n] +
											6*out->cvNeurons[n]/gsl_pow_2(out->rateNeurons[n]) -
											3/gsl_pow_4(out->rateNeurons[n]);

					// standardize by divding with the 4th power of the standard deviation
					out->kurtosisNeurons[n] /= gsl_pow_2(out->cvNeurons[n] - gsl_pow_2(1./out->rateNeurons[n]));
					// no break, since all lower moments need to be calculated as well

				case 3 :
					//! Calculate the 3rd standardized moment (skewness)
					out->skewnessNeurons[n] += -3*out->cvNeurons[n]/out->rateNeurons[n] + 2/gsl_pow_3(out->rateNeurons[n]);

					// standardize by divding with the 3rd power of the standard deviation
					out->skewnessNeurons[n] /= pow(out->cvNeurons[n] - gsl_pow_2(1./out->rateNeurons[n]), (reell)3./2);
					//no break

				case 2 :
					//! Calculate the coefficient of variation
					out->cvNeurons[n] = sqrt(out->cvNeurons[n]*gsl_pow_2(out->rateNeurons[n]) - 1);
					//no break
			}

		}
		else
		{
			switch (in->ISIstats)
			{
				case 4 :
					out->kurtosisNeurons[n] = GSL_NAN;
					// no break, since all lower moments need to be calculated as well

				case 3:
					out->skewnessNeurons[n] = GSL_NAN;
					// no break

				case 2 :
					out->cvNeurons[n] = GSL_NAN;
					// no break

				default :
					out->rateNeurons[n] = GSL_NAN;
					// no break
			}
		}
	}

	//calculate the distributions
	if (in->ISIbins > 0)
		switch (in->ISIstats)
		{
			case 4 :
				out->kurtosisDist = histogram_uniform(out->kurtosisNeurons, in->ISIbins);
				// no break, since all lower moments need to be calculated as well

			case 3 :
				out->skewnessDist = histogram_uniform(out->skewnessNeurons, in->ISIbins);
				// no break

			case 2 :
				out->cvDist = histogram_uniform(out->cvNeurons, in->ISIbins);
				// no break

			default :
				out->rateDist = histogram_uniform(out->rateNeurons, in->ISIbins);
				// no break
		}

}


void LE_analysis::pre_ISIdecorr_statistics(st_in* in, st_out* out)
{
	if (in->ISIdecorr_stats)
	{
		out->ISIdecorr = vector<vector<reell> > (in->ISIdecorr.size(),vector<reell> (2,0));
		out->cvISIdecorr = vector<vector<reell> > (in->ISIdecorr.size(),vector<reell> (2,0));
		
		out->count_ISIdecorr = vector<vector<int> > (in->ISIdecorr.size(), vector<int> (2,0));
	}
	else
	{
		out->ISIdecorr.clear();
		out->cvISIdecorr.clear();
	}
}


void LE_analysis::calc_ISIdecorr_statistics(st_in* in, st_out* out)
{
	net->put_spiketimes(in, out);
}

void LE_analysis::post_ISIdecorr_statistics(st_in* in, st_out* out)
{
	for (unsigned n=0;n<in->ISIdecorr.size();n++)
	{
		for (int i=0;i<2;i++)
		{
			if (out->count_ISIdecorr[n][i]>0)
			{
				out->ISIdecorr[n][i] /= out->count_ISIdecorr[n][i];
				out->cvISIdecorr[n][i] /= out->count_ISIdecorr[n][i];
// 				cout << out->count_ISIdecorr[n][i] << endl;
				
				out->cvISIdecorr[n][i] = sqrt(out->cvISIdecorr[n][i]/gsl_pow_2(out->ISIdecorr[n][i]));
			}
			else
			{
				out->ISIdecorr[n][i] = GSL_NAN;
				out->cvISIdecorr[n][i] = GSL_NAN;
			}
		}
	}
}


void LE_analysis::post_synchrony(st_in* in, st_out* out)
{
	reell STD_network;
	reell STD_single = 0;
	reell sum_network_time = 0;
	reell sum2_network_time = 0;
	
	vector<reell> sum_single(in->Nrec,0);
	vector<reell> sum2_single(in->Nrec,0);
	vector<reell> sum_network(in->measureTimes.size(),0);
	
	for (unsigned t=0;t<in->measureTimes.size();t++)
	{
		for (int n=0;n<in->Nrec;n++)
		{
			sum_network[t] += out->measureStates[t][n];
			
			sum_single[n] += out->measureStates[t][n];
			sum2_single[n] += gsl_pow_2(out->measureStates[t][n]);
		}
		
		sum_network_time += sum_network[t]/in->Nrec;
		sum2_network_time += gsl_pow_2(sum_network[t]/in->Nrec);		
	}
	
	STD_network = sqrt(sum2_network_time/in->measureTimes.size()-gsl_pow_2(sum_network_time/in->measureTimes.size()));

	for (int n=0;n<in->Nrec;n++)
		STD_single += sqrt(sum2_single[n]/in->measureTimes.size() - gsl_pow_2(sum_single[n]/in->measureTimes.size()));
	
	STD_single /= in->Nrec;
	
	out->chi = STD_network/STD_single;
}


void LE_analysis::pre_measurement(st_in* in, st_out* out)
{
	if (in->calc_train)
		out->spikeTrain.clear();		// reset train vector at beginning of measurements
	
	if (in->pertSize > 0)
	{
// // 		//! evolve all neurons until after a spike
// // 		net->evolve_dt(net->find_nextSpikes());
// // 		net->evolve_spike(false);
		
		if (!applyPerturbation)
		{
			if (myID == 0) cout << endl << "saving all neurons states ... " << endl;

			//! save the original state
			stateOriginal = net->get_PhRepState_Nrec();
			IextOriginal = net->get_externalCurrents_Nloc();
		}
		else
		{
			//! reset to the original state
			net->set_State_Nrec(stateOriginal);
			net->set_externalCurrents_Nrec(IextOriginal);
			unsigned pert_try = 0;
			int pert_fail = 1;
			
			vector<vector<reell> > pertState(in->Nrec,vector<reell>(2,0));

			while ((pert_fail==1) && (pert_try<100))
			{
				reell projTmp = 0;
				reell normTmp = 0;
				vector<reell> pertTmp (in->Nrec);
				
				for (int n=0; n<in->Nrec; n++)
				{
					pertTmp[n] = gsl_ran_flat(in->rng, -1,1);
					projTmp += pertTmp[n];
				}
				
				for (int n=0; n<in->Nrec; n++)
				{
					pertTmp[n] -= projTmp * 1/in->Nrec;			// project on normalized trajectory
					normTmp += gsl_pow_2(pertTmp[n]);
				}
				normTmp = sqrt(normTmp);
				
				for (int n=0; n<in->Nrec; n++)
				{
					pertTmp[n] *= in->pertSize/normTmp;
					pertState[n][0] = stateOriginal[n][0] + pertTmp[n];
					
					if (pertState[n][0] > 1)
					{
						pert_try++;
						break;			// when perturbed above threshold
					}
					if (n == in->Nrec-1) pert_fail = 0;			// when all perturbations are applied and no errors
				}
			}
			
			if (pert_fail > 0)
			{
				out->CDB_state = true;
				out->PTdistance.push_back(GSL_NAN);
				out->PTtime.push_back(GSL_NAN);
			}
			else
				net->set_State_Nrec(pertState);
			
			if (in->CDB>2)
			{
				vector<reell> vect_tmp (0);
				out->distances.push_back(vect_tmp);
				cout << "vector size: " << out->distances.size() << endl;
			}
			
			pertState.clear();
			
			//! test for proper perturbation
			vector<vector<reell> > stateTmp = net->get_PhRepState_Nrec();
			
			reell deltaState;
			reell projTmp = 0;
			reell distance = 0;
			for (unsigned n=0; n<stateOriginal.size(); n++)
			{
				deltaState = stateTmp[n][0] - stateOriginal[n][0];
// 				cout << "delta phase neuron " << n << ": " << deltaState << endl;
				projTmp += deltaState;
				distance += gsl_pow_2(deltaState);
			}
			cout << "total projection on 1,1,1,...: " << projTmp << endl;
			cout << "distance at beginning: " << sqrt(distance) << endl;
// 			
		}
	}
	
	if (in->pertSpike || in->pertSynapse)
	{
		//! evolve all neurons until the next // recurrent // spike

		if (!applyPerturbation)
		{
			bool rec_spike = false;
			while (!rec_spike)
			{
				net->evolve_dt(net->find_nextSpikes());
				vector<st_spike> spikes = net->get_spikes();
				if (spikes[0].neuron < in->Nrec)
				{
					rec_spike = true;
					cout << "evolved up to spike from neuron " << spikes[0].neuron << endl;
				}
				else
				{
					net->evolve_spike(false);
					cout << "spike from neuron " << spikes[0].neuron << " passed" << endl;
				}
			}
		
			if (myID == 0) cout << endl << "saving all neurons states ... " << endl;

			//! save the original state after the spike transmission
			stateOriginal = net->get_PhRepState_Nrec();
			IextOriginal = net->get_externalCurrents_Nloc();
// 			stateOriginal = net->get_simState_Nloc();
// 			IextOriginal = net->get_externalCurrents_Nloc();
			vector<st_spike> spikes = net->get_spikes();
			cout << "spiking neuron: " << spikes[0].neuron << endl;
			//! Let the next neuron spike
			net->evolve_spike(false);

		}
		else
		{
			//!reset to the original state
			net->set_State_Nrec(stateOriginal);
			net->set_externalCurrents_Nrec(IextOriginal);
// 			net->set_simState_Nloc(stateOriginal);
// 			net->set_externalCurrents_Nloc(IextOriginal);
			net->evolve_dt(net->find_nextSpikes());


			//! Don't let the neuron spike but reset the spiking neuron.

			//! Get the next spiking neuron and save its synapses.
			vector<st_spike> spikes = net->get_spikes();
			vector<st_synapse> postOrg = in->synapses[spikes[0].neuron];

			if (in->pertSpike)
			{
				if (myID == 0) cout << endl << "skipping one spike ... " << endl;
				cout << "of neuron " << spikes[0].neuron << endl;
				//! Delete all synapses of this neuron temporarily.
				in->synapses[spikes[0].neuron] = vector<st_synapse> (0);
			}
			else //if (in->pertSynapse)
			{
				if (myID == 0) cout << endl << "skipping one synaptic transmission ... " << endl;

				//! Delete the first synapse of this neuron temporarily.
				in->synapses[spikes[0].neuron].erase(in->synapses[spikes[0].neuron].begin());
			}

			//! Run the spike update. Since there are no postsynaptic neurons for the spiking neuron stored, this spike will be skipped, but the spiking neuron reset.
			net->evolve_spike(false);

			//! Reset the postsynaptic connections to the original state.
			in->synapses[spikes[0].neuron] = postOrg;
		}
	}

	//! Prepare stuff to measure
	measureTimeCount = 0;
	measureSpikeCount = 0;

	if (in->measureSpikes.size())
	{
		in->measureTimes.clear();
		out->measure1stStateVarNeurons.resize(in->measureSpikes.size());
		out->measure2ndStateVarNeurons.resize(in->measureSpikes.size());
		
		if (applyPerturbation)
		{
			out->measureStatesPert.resize(out->PT_steps);
			for (int p=0; p<out->PT_steps; p++)
				out->measureStatesPert[p].resize(in->measureSpikes.size());
		}
		else 
		{
			out->measureStates.resize(in->measureSpikes.size());
			if (myID == 0) out->measureTimes.resize(in->measureSpikes.size());	//times need only be saved on root}
		}
		
		for (unsigned s=0; s<in->measureSpikes.size(); s++)
		{
			out->measure1stStateVarNeurons[s].resize(in->Nrec);
			out->measure2ndStateVarNeurons[s].resize(in->Nrec);
			if (applyPerturbation)
				for (int p=0; p<out->PT_steps; p++)
					out->measureStatesPert[p][s].resize(in->Nrec);
			else out->measureStates[s].resize(in->Nrec);
		}
	}
		
	else if (in->measureTimes.size())
	{
		out->measure1stStateVarNeurons.resize(in->measureTimes.size());
		out->measure2ndStateVarNeurons.resize(in->measureTimes.size());
		
		if (applyPerturbation)
		{
			out->measureStatesPert.resize(out->PT_steps);
			for (int p=0; p<out->PT_steps; p++)
				out->measureStatesPert[p].resize(in->measureTimes.size());
		}
		else 
		{
			out->measureStates.resize(in->measureTimes.size());
			if (myID == 0) out->measureTimes.resize(in->measureTimes.size());	//times need only be saved on root
		}
		
		for (unsigned t=0; t<in->measureTimes.size(); t++)
		{
			out->measure1stStateVarNeurons[t].resize(in->Nrec);
			out->measure2ndStateVarNeurons[t].resize(in->Nrec);
			if (applyPerturbation)
				for (int p=0; p<out->PT_steps; p++)
					out->measureStatesPert[p][t].resize(in->Nrec);
			else out->measureStates[t].resize(in->Nrec);
		}
	}
	else
	{
		out->measure1stStateVarNeurons.clear();
		out->measure2ndStateVarNeurons.clear();
		out->measureStates.clear();
		out->measureStatesPert.clear();
		out->measureTimes.clear();
	}
	
	//! Save the first phase/state if measureSpikes == 0.
	vector<vector<reell> > stateTmp;
	if ((in->measureSpikes.size()) && (in->measureSpikes[measureSpikeCount] == 0))
	{
// 		if ( in->measures == 1 )
// 		{
// 			stateTmp = net->get_VRepState_Nloc();
// 			for (int n=0; n<in->Nrec; n++)
// 			{
// 				out->measure1stStateVarNeurons[measureSpikeCount][n]=stateTmp[n][0];
// 				out->measure2ndStateVarNeurons[measureSpikeCount][n]=stateTmp[n][1];
// 			}
// 		}
		if ( in->measures == 2 )
		{
			stateTmp = net->get_PhRepState_Nrec();
			for (int n=0; n<in->Nrec; n++)
			{
				out->measure1stStateVarNeurons[measureSpikeCount][n]=stateTmp[n][0];
				out->measure2ndStateVarNeurons[measureSpikeCount][n]=stateTmp[n][1];
			}
		}
		
		if (myID == 0 and !applyPerturbation) out->measureTimes[measureSpikeCount] = out->TC;

		//! go to the next provided spike number
		if (measureSpikeCount < in->measureSpikes.size() - 1)
			measureSpikeCount++;
	}
}



void LE_analysis::warmup_ONS(st_out* out)
{
	//! Find a good step size for the orthonormalizations and warmup the ONS with SW spikes


	reell condNo = 0;
	reell condNoSub = 0;
	reell condMax = 42;			// 8-)

	int NLE = out->LyapunovExponentsONS.size();
	int subNLE = 0;
	if (out->subN>0)
		subNLE = out->subLyapunovExponentsONS.size();

	if (((NLE > 1) || (subNLE > 1) ) && (out->ONstep > 1))
	{
		if (myID ==0 )
			cout << endl << "optimize ON step size (maximal conditional number = " << condMax << ")..." << endl;

		do
		{
			//! Find the local ONS by reiterating the orthonormalization three times
			reell dt = net->find_nextSpikes();
			net->evolve_dt(dt);
			net->evolve_spike(true);

			for(int iter=0; iter<3; iter++)
			{
				net->multiply_JacobianONS(ons, NLE);
				ons->orthonormalizeONS(NLE);
				
				if (out->subN>0)
				{
					net->multiply_JacobianONS(sub_ons, subNLE);
					sub_ons->orthonormalizeONS(subNLE);
				}
			}

			//! Evaluate the condition number of the Jacobian matrix from 6 calculations with a certain step size
			//! If the condition number is too large (here larger than 42) than the steps size needs to be smaller
			condNo = 0;
			condNoSub = 0;	
			for (int bla=0; bla<6; bla++)
			{
				for(int st=0; st<out->ONstep; st++)
				{
					reell dt = net->find_nextSpikes();
					net->evolve_dt(dt);
					net->evolve_spike(true);

					net->multiply_JacobianONS(ons, NLE);
					if (out->subN>0)
						net->multiply_JacobianONS(sub_ons, subNLE);
				}

				ons->orthonormalizeONS(NLE);
				if (out->subN>0)
					sub_ons->orthonormalizeONS(out->subN);

				condNo += ons->get_normONS(0)/ons->get_normONS(NLE-1);
			}
			
			condNo /= 6;
			if (out->subN>0)
			{
				condNoSub /= 6;
				condNo=max(condNo,condNoSub);
			}
			
			if (myID == 0) cout << "\tthe average condition number with ON steps = " << out->ONstep << " is " << condNo << endl;

			if (condNo > condMax)
				out->ONstep = (int)floor(out->ONstep/2);

		}
		while ((out->ONstep > 1) && (condNo > condMax));

		if (myID ==0 ) cout << "\tsetting ON step size to " << out->ONstep << endl;
	}



	//! warmup of the ONS with the found step size with 1 spike per neuron
	if (myID == 0) cout << endl << "warmup of the orthonormal system with " << out->SWONS << " spikes ..." << endl;

	int long long s = 0;
	while (s < out->SWONS)
	{
		for (int st=0; st<out->ONstep; st++)
		{
			reell dt = net->find_nextSpikes();
			net->evolve_dt(dt);
			net->evolve_spike(true);

			net->multiply_JacobianONS(ons, NLE);
			if (out->subN>0)
				net->multiply_JacobianONS(sub_ons, subNLE);

			s++;
		}

		ons->orthonormalizeONS(NLE);
		if (out->subN>0)
			sub_ons->orthonormalizeONS(subNLE);
	}

}

void LE_analysis::pre_LyapunovExponents(st_in* in, st_out* out)
{
	//! initialize the orthonormal system
	//! The state dimension is the same for all neurons, since net->get_allNeuronSame must be true for this calculation. This is checked above.
	delete ons;
	ons = new LE_ons(in->LyapunovExponents, net->get_Nloc(), net->get_stateDim(0), in->seedONS + myID);
	if (out->subN>0)
	{
		delete sub_ons;
		sub_ons = new LE_ons(in->subLyapunovExponents, in->subNloc, net->get_stateDim(0), in->seedONS + myID);	
	}
	
	//! zero vector for Lyapunov exponents in the out structure (on all nodes including the slaves, because NLE is extract from LyapunovExponents.size() later on)
	out->LyapunovExponentsONS = vector<reell> (in->LyapunovExponents);
	if (out->subN>0)
		out->subLyapunovExponentsONS = vector<reell> (in->subLyapunovExponents);

	//! initial values for the convergence calculation
	LEmaxOld = 0;
	LEminOld = 0;
	out->pLEONS = in->pLEONS;

	//! the orthonormalization step size will be adapted and the ONS warmuped
	out->ONstep = (in->ONstep > 0) ? in->ONstep : 1;
	out->SWONS = in->SWONS;

	if (in->LyapunovExponents > 0)
		warmup_ONS(out);

	//! Important: The result during the warmup are excluded in the following calculations including the convergence vector.

	//! Should the convergence of the Lyapunov exponents be saved (necessary only on root)?
	if ((myID == 0) && in->LyapunovExponentsConvergence)
	{
		//! initialize and assign the first output vector element with 0
		out->LEconvergence = vector<vector<reell> > (1, out->LyapunovExponentsONS);
		if (out->subN>0)
			out->subLEconvergence = vector<vector<reell> > (1, out->subLyapunovExponentsONS);
		//! Reserve the minimal amount of memory for the number of orthonormalizations in case the number of spikes per neuron were provided.
		//! If the time was provided this can't be done due to the lack of knowledge of the minimal number of orthonormalizations.
		if (in->SC)
		{
			out->LEconvergence.reserve(in->SC/out->ONstep);
			if (out->subN>0)
				out->subLEconvergence.reserve(in->SC/out->ONstep);
		}

	}
	else
	{
		out->LEconvergence = vector<vector<reell> > (0);
		if (out->subN)
			out->subLEconvergence = vector<vector<reell> > (0);
		out->LEtimes = vector<reell> (0);
	}

	//! Store the times of orthonormalization if the convergence is saved or the covariant Lyapunov exponents are computed
	if ((myID == 0) && (in->LyapunovExponentsConvergence || in->CLV))
	{
		out->LEtimes = vector<reell> (1);
		if (in->SC)
			out->LEtimes.reserve(in->SC/out->ONstep);
	}


	//! If the covariant Lyapunov vectors are to be calculated the projection matrix R needs to be stored for each orthonormalization.
	//! This should happen after the warmup of the ons, since this should not be used for the covariant Lyapunov vector calculation.
	if ((myID == 0) && in->CLV)
	{
		ons->saveProjections = true;
		if (out->subN>0)
			sub_ons->saveProjections = true;
		
		//! Initialize the covariant Lyapunov vectors.
		ons->initCLV();
		if (out->subN>0)
			sub_ons->initCLV();
		
		out->SWCLV = in->SWCLV;
		out->LyapunovExponentsCLV = vector<reell> (in->LyapunovExponents);
		if (out->subN>0)
			out->subLyapunovExponentsCLV = vector<reell> (in->subLyapunovExponents);
		
		out->localLyapunovExponents = vector<vector<reell> > (1, out->LyapunovExponentsCLV);
		if (out->subN>0)
			out->sublocalLyapunovExponents = vector<vector<reell> > (1, out->subLyapunovExponentsCLV);
		
		if (in->SC)
		{
			out->localLyapunovExponents.reserve((in->SC - in->SWCLV)/out->ONstep);
			if (out->subN>0)
				out->sublocalLyapunovExponents.reserve((in->SC - in->SWCLV)/out->ONstep);
		}
		
		
			
	}
	else
	{
		ons->saveProjections = false;
		if (out->subN>0)
			sub_ons->saveProjections = false;
	}
}

void LE_analysis::calc_LyapunovExponents(st_out* out)
{
	int NLE = out->LyapunovExponentsONS.size();
	int subNLE = 0;
	if (out->subN>0)
		subNLE = out->subLyapunovExponentsONS.size();
	
	net->multiply_JacobianONS(ons, NLE);
	if (out->subN>0)
		net->multiply_JacobianONS(sub_ons, subNLE);
	
	//! Reorthonormalize only every other step to save computation time
	if(out->SC % out->ONstep == 0)
	{
		ons->orthonormalizeONS(NLE);
		
		for(int n=0; n<NLE; n++)
			if (myID == 0) out->LyapunovExponentsONS[n] += log(ons->get_normONS(n));
		
		if (out->subN>0)
		{
			sub_ons->orthonormalizeONS(subNLE);
		
			for(int n=0; n<subNLE; n++)
				if (myID == 0) out->subLyapunovExponentsONS[n] += log(sub_ons->get_normONS(n));
		}

		//! Should the convergence be monitored?
		if ((myID == 0) && (out->LEtimes.size()))
		{
			//save the current times
			out->LEtimes.push_back(out->TC);

			//save the accumulated log(norms) (divide by the time in postprocessing, then the element index is save to use)
			if (out->LyapunovExponentsConvergence)
			{
				cout << "why save? LEconv: " << out->LyapunovExponentsONS[0] << endl;
				out->LEconvergence.push_back(out->LyapunovExponentsONS);
				if (out->subN>0)
					out->subLEconvergence.push_back(out->subLyapunovExponentsONS);
			}
		}


		//! Check the convergence of the largest and the smallest Lyapunov exponent (3 ONsteps back in time)
		if ((out->pLEONS > 0) && (int(out->SC/out->ONstep) % 3 == 0) )
		{
			reell LEmax = out->LyapunovExponentsONS[0]/out->TC;
			reell LEmin = out->LyapunovExponentsONS[NLE-1]/out->TC;

			bool notConvMax = (fabs(LEmax) > 1e-10) ? (fabs(1 - LEmaxOld/LEmax) > out->pLEONS) : false;
			bool notConvMin = (fabs(LEmin) > 1e-10) ? (fabs(1 - LEminOld/LEmin) > out->pLEONS) : false;

			notConverged = notConvMax || notConvMin;

#ifdef PAR
			int sum=(notConverged)? 1:0;
			MPI_Reduce(&notConverged, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  //mpi has no bool type :(
			notConverged=(sum>0);
#endif
			
			LEmaxOld = LEmax;
			LEminOld = LEmin;
		}
	}
}

void LE_analysis::post_LyapunovExponents(st_in* in, st_out* out)
{
	for(unsigned n=0; n<out->LyapunovExponentsONS.size(); n++)
		out->LyapunovExponentsONS[n] /= out->TC;
	if (out->subN>0)
	{
		for(unsigned n=0; n<out->subLyapunovExponentsONS.size(); n++)
			out->subLyapunovExponentsONS[n] /= out->TC; 
	}

	for (unsigned t=1; t<out->LEconvergence.size(); t++)
	{
		// divide the stored log(norms) by the time to yield the Lyapunov exponents (not the first element, it's zero)
		for(unsigned n=0; n<out->LEconvergence[t].size(); n++)
			out->LEconvergence[t][n] /= out->LEtimes[t];
		
		if (out->subN>0)
		{
			for(unsigned n=0; n<out->subLEconvergence[t].size(); n++)
				out->subLEconvergence[t][n] /= out->LEtimes[t];
		}
	}
	
// 	ons->printFinalONS();
// 	if (out->subN>0)
// 		sub_ons->printFinalONS();

}

void LE_analysis::calc_LyapunovVectors(st_out* out)
{
	if (myID == 0)
	{
		//! Get the number of orthonormalizations for which the projections were stored.
		unsigned SC = ons->get_numberSavedProjections();
		int NLE = out->LyapunovExponentsONS.size();
		int subNLE = 0;
		if (out->subN>0)
			subNLE = out->subLyapunovExponentsONS.size();
		
		unsigned SW = out->SWCLV/out->ONstep;

		cout << "calculating the covariant Lyapunov vectors from " << SC << " orthonormalization steps (this always runs on only 1 processor) ..." << endl;

		int Sdone = 1;
		//! Run the covariant Lyapunov vector calculation.
		for (unsigned s=1; s<=SC; s++)
		{

			ons->backwardIterationCLV(SC-s, NLE);
			if (out->subN>0)
				sub_ons->backwardIterationCLV(SC-s, subNLE);
			//! After the warmup of the CLVs start the analysis
			if (s > SW)
			{
				//! Get the log norms of the current orthonormalization of the CLVs
				vector<reell> logNorms(NLE);
				for(int n=0; n<NLE; n++)
					logNorms[n] = -log(ons->get_normCLV(n));		// -(minus) because time is reversed
				vector<reell> sublogNorms(subNLE);
				if (out->subN>0)
				{
					for(int n=0; n<subNLE; n++)
						sublogNorms[n] = -log(sub_ons->get_normCLV(n));		// -(minus) because time is reversed
				}
				//! store the local Lyapunov exponents.
				//! They are in reverse order, the 0 element is the first after the warmup of the CLVs. Change order later.
				out->localLyapunovExponents.push_back(logNorms);
				if (out->subN>0)
					out->sublocalLyapunovExponents.push_back(sublogNorms);
				
				//! store the backward Lyapunov exponents.
				for(int n=0; n<NLE; n++)
					out->LyapunovExponentsCLV[n] += logNorms[n];
				if (out->subN>0)
				{
					for(int n=0; n<subNLE; n++)
						out->subLyapunovExponentsCLV[n] += sublogNorms[n];
				}
				  
			}

			if (s == SW)
				cout << "\twarmup done, step: " << s << endl;

			if  (s >= SC/10.*Sdone)
				if (myID == 0) cout << "\t" << Sdone++*10 << "% done ... step: " << s << endl;


		}

		//! Postprocessing of the results

		//! Calculate the backward LEs.
		reell timeCLV = out->LEtimes[SC - (SW+1)] - out->LEtimes[0];
		for(int n=0; n<NLE; n++)
			out->LyapunovExponentsCLV[n] /= timeCLV;
		if (out->subN>0)
		{
			for(int n=0; n<subNLE; n++)
				out->subLyapunovExponentsCLV[n] /= timeCLV;
		}
		
		//! Swap order of the local LEs.
		for (unsigned s=0; s<out->localLyapunovExponents.size()/2; s++)
			out->localLyapunovExponents[s].swap(out->localLyapunovExponents[out->localLyapunovExponents.size() - (s+1)]);
		if (out->subN>0)
		{
			for (unsigned s=0; s<out->sublocalLyapunovExponents.size()/2; s++)
				out->sublocalLyapunovExponents[s].swap(out->sublocalLyapunovExponents[out->sublocalLyapunovExponents.size() - (s+1)]);
		}
		
	}
	else
	{
		cout << "The covariant Lyapunov vector calculation is completely done on the root process, but was called on a slave here" << endl;
		throw(2);
	}
}

void LE_analysis::calc_randomLyapunovExponents(st_in* in, st_out* out)
{

	//! Initialize the variables
	int NLE = out->LyapunovExponentsONS.size();
	int SC = stateTrajectory.size();

	out->randomLyapunovExponentsONS.resize(in->randomizedLEspectra);
	for (int i=0; i<in->randomizedLEspectra;i++)
		out->randomLyapunovExponentsONS[i].resize(NLE);
	
	//!< permutation environment
	gsl_rng_env_setup();
	const gsl_rng_type *TYPE = gsl_rng_mt19937;
	gsl_rng * rng = gsl_rng_alloc(TYPE);
	gsl_rng_set(rng, 0);


	int permutation[SC];          
        for (int i = 0; i < SC; i++)
		permutation[i] = i;
	
	//! Compute the random LE
	for( int i=0;i< in->randomizedLEspectra;i++)
	{	
		if (myID==0)
		{
			cout<< endl << "computing " << i << " of " << in->randomizedLEspectra << " state sequence randomized spectra" <<endl;
		
			gsl_ran_shuffle (rng, permutation, SC, sizeof (int));
		}
		
#ifdef PAR
		MPI_Bcast(&permutation, SC, MPI_INT, 0, MPI_COMM_WORLD);
#endif
		
		
		int Sdone = 1;
		
		ons->initONS(in->seedONS + myID);
		
		for (int s=0;s< SC;s++)
		{
			if (s >= SC/10.*Sdone)
				if (myID == 0) cout << "\t" << Sdone++*10 << "% of SC done ... " << s << " spikes "<< endl;

			//! need to set phases and compute jacobians
	    
			net->set_simState_Nloc(stateTrajectory.at(permutation[s]));

			net->find_nextSpikes();
			
			net->evolve_spike(true);
	    
			net->multiply_JacobianONS(ons, NLE);

			ons->orthonormalizeONS(NLE);
		    
			for(int n=0; n<NLE; n++)
				if (myID == 0) out->randomLyapunovExponentsONS[i][n] += log(ons->get_normONS(n));

		}
		
	}
	
	//! postprocessing
	for( int i=0;i< in->randomizedLEspectra;i++)
		for(int n=0; n< NLE; n++)
			out->randomLyapunovExponentsONS[i][n] /= out->TC;
}

vector<vector<reell> > LE_analysis::histogram_uniform(vector<reell> data, int bins)
{
	//! Calculates the histogram with uniformly spaced bins on the x axis, the x-values return the center of the bins.

	//! Erase NANs from the distribution
	unsigned x = 0;
	do
	{
		if (gsl_isnan(data[x]))
			data.erase(data.begin() + x);
		else
			x++;
	}
	while (x < data.size());

	//! sort data in ascending order
	sort(data.begin(), data.end());

	//! set the uniformly distributed x values
	vector<vector<reell> > hist(bins, vector<reell> (2));

	reell dataMin = data[0];
	reell dataMax = data[data.size() - 1];

	reell dx = (dataMax - dataMin)/bins;
	reell dx2 = dx/2;

	for(int b=0; b<bins; b++)
	{
		hist[b][0] = dataMin + b*dx + dx2;					//centers of the bins
		hist[b][1] = 0;
	}

	//! fill the histogram
	vector<reell>::iterator it = data.begin();

	for(int b=0; b<bins; b++)
		while (*it <= hist[b][0] + dx2)
		{
			hist[b][1] = hist[b][1] + 1;

			if (it < data.end() - 1)
				it++;
			else
				break;			//shouldn't occur, since hist[b][bins-1] == maxData
		}

	//! normalize the histogram
	for(int b=0; b<bins; b++)
		hist[b][1] /= dx*data.size();				//all bins are equally space here (dx)!

	return hist;
}

/** possible extensions:
 *  logarithmic histograms
 */
