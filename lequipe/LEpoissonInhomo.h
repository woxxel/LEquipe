/*
 * LEpoisson.h
 *
 *  Created on: 04.05.2012
 *      Author: ryanAir
 */

#ifndef LEPOISSONINHOMO_H_
#define LEPOISSONINHOMO_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>
#include "gsl/gsl_rng.h"
#include <gsl/gsl_randist.h>
using namespace std;

//! The poisson class

class LE_poissonInhomo : public LE_phaseneuron
{
public:
    LE_poissonInhomo(reell r, reell tau, reell poissonrate, vector< reell > init, reell *TC_);
    virtual ~LE_poissonInhomo();

    virtual void set_externalCurrent(reell&);

protected:
    virtual void calc_spikeTime();
    virtual void evolve_dt(reell, vector<reell>*);
    inline virtual reell PTC(reell);							//!< phase transition curve
    inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve
	virtual vector<reell> voltage2phase(vector<reell>&);						
    virtual vector<reell> phase2voltage(vector<reell>&);	    
	virtual void reset(); 								//! < the new spiketime of the poisson neuron is generated in reset
    gsl_rng* rng;								//!< random number generator
    reell poissonrate;
    reell poissonrateInh;
    reell *TC;
    reell fIn;												//!<  rap is elsewhere used for onset rapidness is in this neuron  used for phase of sine). Ugly, but new variables would imply a different Hash.


};

#endif /* LEPOISSONINHOMO_H_ */