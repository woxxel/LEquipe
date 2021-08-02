/*
* LElif.cpp
*
*  Created on: 04.05.2012
*      Author: ryanAir
*/
#include "LEpoissonInhomo.h"

LE_poissonInhomo::LE_poissonInhomo(reell r, reell tau, reell p, vector< reell > init,reell *TC_) : LE_phaseneuron(tau, vector<reell> (1, 0), vector<reell> (1, 1)),
        TC(TC_)
{
    //! initialize random number generator
    gsl_rng_env_setup();
    const gsl_rng_type *TYPE = gsl_rng_mt19937;
    rng = gsl_rng_alloc(TYPE);
    gsl_rng_set(rng, 1 +init[0]);
    poissonrate = p;
    fIn = r; // here we use the variable rapidness to transmit the value of the frequency of the sine of the poissonrate.
    reset();
}

LE_poissonInhomo::~LE_poissonInhomo()
{
}

//! The external current sets the sine gain.
void LE_poissonInhomo::set_externalCurrent(reell& Iext)
{
    this->Iext = Iext;
    calcSpikeTime = true;
}
inline reell LE_poissonInhomo::PTC(reell c)
{
    cout << "A poisson neuron should never receive a spike. Something about your network topology is wrong!" << endl;
    throw(1);
    return 0;
}

inline reell LE_poissonInhomo::PTCprime(reell c)
{
    cout << "A poisson neuron should never receive a spike. Something about your network topology is wrong!" << endl;
    cout << "Calculation of Lyapunov spectrum impossible for network with poisson neuron!" << endl;
    throw(1);
    return 0;
}

vector<reell> LE_poissonInhomo::phase2voltage(vector<reell>& phState)
{
  		cout << "phase/voltage conversion not implemented!" << endl;
		throw(1);
		
    return phState;

}

vector<reell> LE_poissonInhomo::voltage2phase(vector<reell>& VState)
{
  		cout << "phase/voltage conversion not implemented!" << endl;
		throw(1);
		
    return VState;

}

inline void LE_poissonInhomo::calc_spikeTime()
{
    calcSpikeTime = false;
}



//! Set the state variable to the reset value.
inline void LE_poissonInhomo::reset()
{
    poissonrateInh = poissonrate *(1+ Iext*sin(fIn/1000*2*M_PI* (*TC)  )) ;
    //!< The output of gsl_ran_flat is on the half-open interval [a,b). For an Output 0, we would get spiketime = inf. Thefore we take (1-random number)
    spikeTime = -log(1-gsl_ran_flat(rng, 0,1)) /(poissonrateInh*tauM/1000); //!< the next spike time of the neuron
    calcSpikeTime = false;
}

void LE_poissonInhomo::evolve_dt(reell dt, vector<reell>* dummy)
{
    dt /= tauM;
    spikeTime -= dt; //! This accumulates precision loss. If that is a problem, it should be changed.
//! For a demonstration run a distance measurement
    //! of a zero vector with spikeTime -= dt and you will see how the distance increases.
    calcSpikeTime = false;
}