/*
 * LEneuron.cpp
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#include "LEneuron.h"

LE_neuron::LE_neuron(reell tau, vector<reell> reset, vector<reell> threshold)
{
	//! All time related quantities are expressed in terms of tauM internally.
	tauM = tau;
	
	//! Every neuron has some reset and threshold values.
	stateReset = reset;
	stateThreshold = threshold;

	Iext = 0;

	calcSpikeTime = false;
	
	spike_time = 0;
	prespike_time = 0;
	
	measure_post_pre = false;
}

LE_neuron::~LE_neuron()
{

}

reell LE_neuron::get_externalCurrent()
{
	return Iext;
}


reell LE_neuron::get_spikeTime()
{
	//! Returns the next spike time of the neuron.
	//! If it is not up to date, then it is calculated first.
	if (calcSpikeTime)
		calc_spikeTime();
	
	if (spikeTime < 0)
	{
		cout << "The neurons next spike time " << spikeTime << " was in the past!" << endl;
		throw(1);
	}
	return spikeTime*tauM;
}

vector<reell> LE_neuron::get_sim_state()
{
	return state;
}

void LE_neuron::set_sim_state(vector<reell>& newState)
{
	state = newState;
// 	if ( newState[0] > stateThreshold[0]) {cout << "state was set by SimRep to above threshold. It is reset to threshold\n";}
	state[0] = (state[0] < stateThreshold[0]) ? state[0] : stateThreshold[0];
	calcSpikeTime = true;
}



reell LE_neuron::update_ISIdecorr_pre_post()				// added by Alex (ISIdecorr), time from prespike to postspike
{
// 	cout << "spike - prespike time: " << spike_time - prespike_time << endl;
	
	ISI_tmp = spike_time - prespike_time;
	
	measure_post_pre = true;				// enable measuring again
	
	return ISI_tmp;
}


reell LE_neuron::update_ISIdecorr_post_pre()				// added by Alex (ISIdecorr), time from postspike to prespike
{
// 	cout << "prespike - spike time: " << prespike_time - spike_time << endl;
	
	ISI_tmp = prespike_time - spike_time;
	
	measure_post_pre = false;				// disable measuring - first ISI was found
	
	return ISI_tmp;
}


//! This function can be overloaded by the derived classes to store some dummy variables.
void LE_neuron::dummy_evolve(reell dt, vector<reell>* dummy)
{

}

void LE_neuron::dummy_PostUpdate()
{

}

void LE_neuron::get_diff_state(vector<reell>* dummy)
{

}
