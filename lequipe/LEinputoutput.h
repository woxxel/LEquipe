/*
 * LEinputoutput.h
 *
 *  Created on: 11.01.2012
 *      Author: mik
 */


// #ifndef LEINPUTOUTPUT_H_
// #define LEINPUTOUTPUT_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <netcdfcpp.h>
#include "structures.h"
// #include "gsl/gsl_rng.h"
#include "gsl/gsl_math.h"
// #include <gsl/gsl_randist.h>
#include <limits>
#ifdef PAR
	#include <mpi.h>										//!< for parallel simulations
#endif

using namespace std;

class LE_input_output
{
	public:
		LE_input_output();
		virtual ~LE_input_output();

		void read_netcdf(string, string, string, string, string, st_in*);
		void write_netcdf(string, st_out*, string, string, string, bool);

		void display_in(st_in*);
		void display_out(st_out*);

		void gather_Nall(vector<reell> *);
		void gather_Nall(vector<int> *);

	protected:
		int myID, nP;
private:
private:
    reell IextScaling;
};



// #endif /* LEINPUTOUTPUT_H_ */
