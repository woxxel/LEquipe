/*
 * LEpuppet.h
 *
 *  Created on: 04.05.2012
 *      Author: Alex
 */

#ifndef LEPUPPET_H_
#define LEPUPPET_H_

#include "LEphaseneuron.h"
#include <iostream>
#include <vector>
// #include "gsl/gsl_rng.h"
// #include <gsl/gsl_randist.h>
using namespace std;

//! The puppet neuron class

class LE_puppet : public LE_phaseneuron
{
public:
    LE_puppet(reell tau, vector<reell> puppetTrain);
    virtual ~LE_puppet();

    virtual void set_externalCurrent(reell&);
    virtual void set_sim_state(vector<reell> &);
    


protected:
    virtual void calc_spikeTime();
    virtual void evolve_dt(reell, vector<reell>*);
    inline virtual reell PTC(reell);							//!< phase transition curve
    inline virtual reell PTCprime(reell);						//!< derivative of phase transition curve
           virtual vector<reell> voltage2phase(vector<reell>&);						
           virtual vector<reell> phase2voltage(vector<reell>&);
    virtual void reset(); 								//! < the new spiketime of the poisson neuron is generated in reset
    
    vector<reell> spikeTrain;  
    int spike_idx;

};

#endif /* LEPUPPET_H_ */