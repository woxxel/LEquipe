/*
* LElif.cpp
*
*  Created on: 04.05.2012
*      Author: ryanAir
*/

#include "LEpoisson.h"

LE_poisson::LE_poisson(reell tau, reell p, vector< reell > init) : LE_phaseneuron(tau, vector<reell> (1, 0), vector<reell> (1, 1))
{
    /* reell poissonrate = poissonrate;*/
    //! initialize random number generator
    gsl_rng_env_setup();
    const gsl_rng_type *TYPE = gsl_rng_mt19937;
    rng = gsl_rng_alloc(TYPE);
//     gsl_rng_set(rng, 1 +init[0]);
    poissonrate = p;
    w=GSL_POSINF;
//     reset();
//     cout << "Neuron set to seed: " << init[0] << ", spikeTime: " << spikeTime << endl;
    
}

LE_poisson::~LE_poisson()
{

}

inline void LE_poisson::set_sim_state(vector<reell>& newState)
{
	gsl_rng_set(rng, 1 + newState[0]);
	reset();
// 	cout << "reinitiated poisson seed to " << newState[0] << ", new spike time: " << spikeTime << endl;
}

//! The external current sets the phase velocity.
void LE_poisson::set_externalCurrent(reell& poissonrate)
{
    this->poissonrate = poissonrate;
    calcSpikeTime = true;
  //!   Make sure, that Iext(poisson neurons) = Poisson rate, otherwise the initial poisson rate is Iext
}

inline reell LE_poisson::PTC(reell c)
{
      		cout << "A poisson neuron should never receive a spike. Something about your network topology is wrong!" << endl;
		throw(1);
    return 0;
}

inline reell LE_poisson::PTCprime(reell c)
{
        		cout << "A poisson neuron should never receive a spike. Something about your network topology is wrong!" << endl;
  		cout << "Calculation of Lyapunov spectrum impossible for network with poisson neuron!" << endl;
		throw(1);

    return 0;
}

vector<reell> LE_poisson::phase2voltage(vector<reell>& phState)
{
  		cout << "phase/voltage conversion not implemented!" << endl;
		throw(1);
		
    return phState;
}

vector<reell> LE_poisson::voltage2phase(vector<reell>& VState)
{
  		cout << "phase/voltage conversion not implemented!" << endl;
		throw(1);
		
    return VState;
}

inline void LE_poisson::calc_spikeTime()
{	
	spikeTime = -log(1-gsl_ran_flat(rng, 0,1)) /(poissonrate*tauM/1000); //!< the next spike time of the neuron
	calcSpikeTime = false;		//! There is only change here due to evolve_dt!
}

//! calculate a new spike time
inline void LE_poisson::reset()
{
  //!< The output of gsl_ran_flat is on the half-open interval [a,b). For an Output 0, we would get spiketime = inf. Thefore we take (1-random number)  
	spikeTime = -log(1-gsl_ran_flat(rng, 0,1)) /(poissonrate*tauM/1000); //!< the next spike time of the neuron 
	calcSpikeTime = false;
}

void LE_poisson::evolve_dt(reell dt, vector<reell>* dummy)
{
	dt /= tauM;
	
	spikeTime -= dt; //! This accumulates precision loss. If that is a problem, it should be changed.
//! For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.

	calcSpikeTime = false;
}

inline void LE_poisson::evolve_spike(reell c)
{
	cout << "poisson neurons should not receive spikes!" << endl;

	throw(1);
}