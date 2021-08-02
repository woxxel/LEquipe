/*
 * LElif.cpp
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */

#include "LElif.h"

LE_lif::LE_lif(reell tau) : LE_phaseneuron(tau, vector<reell> (1, 0), vector<reell> (1, 1))
{
	//! The rheobase current is always 1 in the LIF model.
	//! The reset value is 0 and the threshold is 1.
	Irheo = 1;
}

LE_lif::~LE_lif()
{

}

//! The external current sets the phase velocity.
void LE_lif::set_externalCurrent(reell& Iext)
{
	this->Iext = Iext;
	w = 1/log(1 + 1/Iext);						//!< Eq.~(5.7)
	calcSpikeTime = true;
}

vector<reell> LE_lif::phase2voltage(vector<reell>& phState)
{
      vector<reell> Vstatetmp(2, 0);
      Vstatetmp[0] = (1 + this->Iext)*(1 - exp(-phState[0]/w));

      return Vstatetmp;
}

vector<reell> LE_lif::voltage2phase(vector<reell>& Vstate)
{
      vector<reell> phStatetmp(2, 0);
      phStatetmp[0] = w*log( (1 + this->Iext)/(1 + this->Iext + Vstate[0]) );
      return phStatetmp;
}
inline reell LE_lif::PTC(reell c)
{
	return -w*log(exp(-state[0]/w) - c/(1 + Iext));
}

inline reell LE_lif::PTCprime(reell c)
{
	reell expTmp = exp(-state[0]/w);
	return expTmp/(expTmp - c/(1 + Iext));
}
