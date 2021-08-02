/*
 * LEphaseneuron.h
 *
 *  Created on: 02.01.2012
 *      Author: mik
 */

#ifndef LEPHASENEURON_H_
#define LEPHASENEURON_H_

#include "LEneuron.h"
#include <iostream>
#include <vector>

using namespace std;




//! The generic phase neuron class.

class LE_phaseneuron : public LE_neuron
{
	public:
		LE_phaseneuron(reell, vector<reell>, vector<reell>);
		virtual ~LE_phaseneuron();

		virtual int get_stateDim();
		virtual vector<reell> get_VRep_state();
		virtual vector<reell> get_PhRep_state();
		
		virtual void set_VRep_state(vector<reell>&);
		virtual void set_PhRep_state(vector<reell>&);
		
		virtual void set_externalCurrent(reell&) = 0;

		virtual void evolve_dt(reell, vector<reell>*);
		virtual void evolve_spike(reell);
		virtual void reset();

		virtual vector<vector<reell> > calc_JacElem_postsynaptic(vector<reell>*, reell);
		virtual vector<vector<reell> > calc_JacElem_self(vector<reell>*);
		virtual vector<vector<reell> > calc_JacElem_spiking(vector<reell>*);
		virtual void dummy_Jacobian_spiking(vector<reell>*);

	protected:
		        void calc_spikeTime();
		virtual reell PTC(reell) = 0;							//!< This method returns the new state value from the model specific phase transition curve.
		virtual reell PTCprime(reell) = 0;						//!< This method returns the derivative of the model specific phase transition curve
		virtual vector<reell> voltage2phase(vector<reell>&) = 0;						
		virtual vector<reell> phase2voltage(vector<reell>&) = 0;						

		reell w;												//!< phase velocity
};

#endif /* LEPHASENEURON_H_ */



