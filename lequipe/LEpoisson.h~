/*
 * LEpoisson.h
 *
 *  Created on: 04.05.2012
 *      Author: ryanAir
 */

#ifndef LEPOISSON_H_
#define LEPOISSON_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>
#include "gsl/gsl_rng.h"
#include <gsl/gsl_randist.h>
using namespace std;

//! The poisson class

class LE_poisson : public LE_phaseneuron
{
public:
    LE_poisson(reell tau, reell poissonrate, vector< reell > init);
    virtual ~LE_poisson();

    virtual void set_externalCurrent(reell&);

protected:
    virtual void calc_spikeTime();
    virtual void evolve_dt(reell, vector<reell>*);
    virtual void set_sim_state(vector<reell> &);
    inline virtual reell PTC(reell);							//!< phase transition curve
    inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve
           virtual vector<reell> voltage2phase(vector<reell>&);						
           virtual vector<reell> phase2voltage(vector<reell>&);
    virtual void reset(); 								//! < the new spiketime of the poisson neuron is generated in reset
    virtual void evolve_spike(reell);

    gsl_rng* rng;								//!< random number generator
    reell poissonrate;



};

#endif /* LEPOISSON_H_ */