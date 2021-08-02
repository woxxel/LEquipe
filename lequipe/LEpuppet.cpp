/*
* LElif.cpp
*
*  Created on: 04.05.2012
*      Author: ryanAir
*/

#include "LEpuppet.h"

LE_puppet::LE_puppet(reell tau, vector<reell> puppetTrain) : LE_phaseneuron(tau, vector<reell> (1, 0), vector<reell> (1, 1))
{
	spikeTrain = puppetTrain;
	spike_idx = 0;
	reset();
}

LE_puppet::~LE_puppet()
{

}

//! The external current has no impact at all!
inline void LE_puppet::set_externalCurrent(reell& poissonrate){}

inline void LE_puppet::set_sim_state(vector<reell>& newState)
{
	state = vector<reell> (1,0);
	spike_idx = 0;
	reset();
}

inline reell LE_puppet::PTC(reell c)
{
      		cout << "A puppet neuron should never receive a spike. Something about your network topology is wrong!" << endl;
		throw(1);
    return 0;


}

inline reell LE_puppet::PTCprime(reell c)
{
		cout << "A puppet neuron should never receive a spike. Something about your network topology is wrong!" << endl;
  		cout << "Calculation of Lyapunov spectrum impossible for network with puppet neuron!" << endl;
		throw(1);

    return 0;

}

vector<reell> LE_puppet::phase2voltage(vector<reell>& phState)
{
	cout << "phase/voltage not available for puppet neurons!" << endl;
	throw(1);
	
	return phState;
}

vector<reell> LE_puppet::voltage2phase(vector<reell>& VState)
{
	cout << "phase/voltage not available for puppet neurons!" << endl;
	throw(1);
	
	return VState;
}

inline void LE_puppet::calc_spikeTime()
{
	cout << "calc_spikeTime is never needed for puppet neurons!" << endl;
	throw(1);
}

//! Set the state variable to the reset value.
inline void LE_puppet::reset()
{
	spikeTime = spikeTrain[spike_idx]/(tauM);
	spike_idx++;
	calcSpikeTime = false;
}

void LE_puppet::evolve_dt(reell dt, vector<reell>* dummy)
{
	dt /= tauM;

	spikeTime -= dt; //! This accumulates precision loss. If that is a problem, it should be changed.
//! For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.
	calcSpikeTime = false;
}