/*
 * LEtwoDlinear.cpp
 *
 *  Created on: 14.02.2012
 *      Author: max & mik
 */

#include "LEtwoDlinear.h"

LE_twoDlinear::LE_twoDlinear(reell tauM, st_twoDlinear paras) : LE_neuron(tauM, vector<reell> (2, 0), vector<reell> (2, 1))
{
	//! we need the machine precision for the calculation of the spike time
	MACHINE_PRECISION = numeric_limits<reell>::epsilon();

	//! The state variable is two-dimensional in the twoDlinear model.
	state = vector<reell> (2);
	
	/**
	 * recall:	alpha	beta	gamma	delta	tauS	Cw
	 * cLIF		 1	0	0	1	0.5	0
	 * RF		-1	1	1	0	2	(0,1)
	 */


	if (paras.tauS == tauM)
	{
		cout << "tauS and tauM cannot be identical." << endl;
		throw(1);
	}

	tauS = paras.tauS/tauM;				//!< tauS in units of tauM
	alp    	=  paras.alpha;
	bet    	=  paras.beta;
	gam    	=  paras.gamma;  			//!correct scaling to map back to conductance models
	deltauS	=  paras.delta/tauS;
	Cw	=  paras.Cw;
	
	//! The rheobase current is always 1 in the LIF model.
	//! The reset value is 0 and the threshold is 1.
	//! for resonator(bet!=0), the rheobase current scales with beta 
	Irheo =1;//+bet*5/8;//1+bet/2;


	//! global variable for the calculation of the real and complex eigenvalues

	//real part
	realPart = -(1 + 1/tauS)/2;

	//complex part
	reell radicandOmega	= realPart*realPart - (1 - alp*bet)/tauS;
	omegaMOD = sqrt(fabs(radicandOmega));
	iscomplex = (radicandOmega < 0) ? true : false;

	//eigenvalues
	lambdaPlus = realPart + omegaMOD;
	lambdaMinus = realPart - omegaMOD;


	//!DIfferential Matrix:used in propspike and storing jacobian biulding blocks, store as vector A=(A11,A12,A21,A22)
	A[0] = -1; 						//!>dissipativve voltage
	A[1] = alp; 						//!>current coupling to voltage
	A[2] = bet/tauS;					//!?voltage couplin gto current
	A[3] = -1/tauS;						//!>disipative current


	//!initialize the 3 biulding blocks of jacobian:
	//!1)the 2x2 linear operator, S(deltat), will be define das dummy
	//!2)two 2x1 vectors, dzdt, 2a) one for spiking and 2b) one for post syn pop, but the latter without J since J is neuron specific
	dzdt[0] = stateThreshold[0]*A[0];// + A[1]*bet/tauS);
	dzdt[1] = stateThreshold[0]*A[2];// + A[3]*bet/tauS);
	dzdt[2] = -(A[0]*gam + A[1]*deltauS);//!withoutJ but with coefficient alpha, since its value depends on neuron...but will only use on post syn pop so doesn't matter...
	dzdt[3] = -(A[2]*gam + A[3]*deltauS);//!withoutJ,
	//!3) dtaudz, will be defined as dummy 
}

LE_twoDlinear::~LE_twoDlinear()
{

}

int LE_twoDlinear::get_stateDim()
{
	return 2;
}

vector<reell> LE_twoDlinear::get_PhRep_state()
{
	    vector<reell> stateTmp (2);
	    
	    reell calcSpikeTimeTmp = calcSpikeTime;
	    reell spikeTimeTmp=spikeTime;
	    
	    reell treset=calc_resetTime();			//! reset crossing time < 0 (current time)
	    calc_spikeTime();					//! threshold crossing time > 0 (current time)
	    
	    stateTmp[1]=1/(spikeTime-treset);			//!> where phstate[1] comes from in set_PhRep_state
	    stateTmp[0]=-treset*stateTmp[1];
	    cout<<"V:"<<state[0]<<" \t W:"<<state[1]<<endl;
	    cout<<"period:"<< stateTmp[1]<<" \t to spike:"<< spikeTime<< "\t from reset:"<< treset<< "\n\n"<< endl;
	    spikeTime=spikeTimeTmp;
	    calcSpikeTime=calcSpikeTimeTmp;
	    return stateTmp; 
}

void LE_twoDlinear::set_PhRep_state(vector<reell>& phstate)
{
	reell dt;
	reell St[4], Sr[4];
	dt=(1-phstate[0])/phstate[1];					//! time to threshold

	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		St[0] = expRdt*( cosOmega + sinOmega*(A[0]-realPart)/omegaMOD );
		St[1] = expRdt*sinOmega*A[1]/omegaMOD;
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		St[0] = coeff*A[0] + tempS;
		St[1] = coeff*A[1];
	}
	
	dt=-phstate[0]/phstate[1];					//! time back to reset
		
	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		Sr[0] = expRdt*( cosOmega + sinOmega*(A[0]-realPart)/omegaMOD );
		Sr[1] = expRdt*sinOmega*A[1]/omegaMOD;
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		Sr[0] = coeff*A[0] + tempS;
		Sr[1] = coeff*A[1];
	}
	
	//! solving equation to threshold for W and then plugging into the equation from reset and solving for V gives:
	state[0] = ( Sr[1] * (stateThreshold[0] - kappaV) + St[1]*kappaV )/(St[0]*Sr[1]-Sr[0]*St[1]);
	state[1] = (stateThreshold[0]-kappaV-St[0]*state[0])/St[1]+kappaW;
	state[0]+= kappaV;
	
	
	if ( state[0] > stateThreshold[0]) {cout << "state was set by PhRep to above threshold. It is reset to thresdold";}
	state[0] = (state[0] < stateThreshold[0]) ? state[0] : stateThreshold[0];
}

vector<reell> LE_twoDlinear::get_VRep_state()
{	    
	    return state; 
}

void LE_twoDlinear::get_diff_state(vector<reell>* dummy)
{
	(*dummy)[0]=A[0]*state[0]+A[1]*state[1]+Cv;
	(*dummy)[1]=A[2]*state[0]+A[3]*state[1]+Cw;
}
void LE_twoDlinear::set_VRep_state(vector<reell>& Vstate)
{	    
	state=Vstate;
	if ( state[0] > stateThreshold[0]) {cout << "state was set by VRep to above threshold. It is reset to thresdold";}
	
	state[0] = (state[0] < stateThreshold[0]) ? state[0] : stateThreshold[0];
}

void LE_twoDlinear::set_externalCurrent(reell& Iext)
{
	this->Iext = Iext;
	
	//!once Iext is defined we can calculate all remaining model parameters.
	
	Cv = Irheo + Iext;

	kappaV = (Cv + alp*Cw)/(1 - alp*bet);
	kappaW = (bet*Cv + Cw)/(1 - alp*bet);
	C3 = (kappaV - stateThreshold[0]);

// 	if C3 less than 0, outside of restrictive conditions for root finding with tguess
// 	if (C3 < 0)
// 	{
// 	    cout << "C3 = " << C3 << " < 0 ! Root finding will fail with kappaV= " << Iext+ Irheo <<endl;
// 	    throw(1);
// 	}

	if (iscomplex==false)
		C3 *= 2; //!just the way C3 is defined in the real case

	calcSpikeTime = true;
}

inline void LE_twoDlinear::calc_spikeTime()
{

	double tol = MACHINE_PRECISION;	//!>function value tolerance
	
	if (state[0] < 1)
	{
		reell C1,C2; 
		C1 = state[0] - kappaV;
		C2 = (alp*(state[1] - kappaW) - (1 + realPart)*C1)/omegaMOD;
		
		//!assign function parameters and define root function object
// 		rootparams.C1=C1;
// 		rootparams.C2=C2;
// 		rootparams.LE2DlinearOBJ=this;
		//gsl_function_fdf FDF;
		//FDF.params = (void*) &rootparams;
		
		fdfval.C1=C1;
		fdfval.C2=C2;
		fdfval.C3=C3;
		
// 			cout << "C1 = " << C1 << endl;
// 			cout << "C2 = " << C2 << endl;
// 			cout << "C3 = " << C3 << endl;
		//!The following assigns the correct (real/complex))root function, and 
		//!computes tguess to fall in interval in which points converge to correct root and run root finfing algorithm
		
		HomemadeRtFn fPoint;
		reell tguess,tlower,tupper;
		bool bisectFlag = false;
		
		if (iscomplex == true)
		{
			fPoint = &LE_twoDlinear::FComplex;
			
			//!for complex eigen values

			//!find extrema between -pi/2 and pi/2
			reell temp1 = omegaMOD*C2 + realPart*C1;
			reell temp2 = realPart*C2 - omegaMOD*C1;
			reell tExtrema = atan( -temp1/temp2 )/omegaMOD;
			
			//!if tExtrema negative then look at next extrema time (which is necessarily positive)
			tExtrema += (tExtrema<0)? M_PI/omegaMOD : 0;
			
			if ((C2 > 0) || (realPart*C1 > -omegaMOD*C2))	//tExtrema is a maximum, 
			{
				tlower = 0;
				tupper = tExtrema;
				
				//is first positive inflection point in this region?
				reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
				t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
				
				if ( (tlower < t_inflect) && (t_inflect < tupper) )
				{
					FComplex(fdfval,t_inflect);
					reell f_inflect = fdfval.fval;
					
					if (f_inflect < -100*tol) 	//place to the right of inflection point
					{
						tlower = t_inflect;
						tguess = tlower+10*tol;
// 										cout<<"Compmax1a: inflect withinbounds and neg"<< endl;

					}
					else if (f_inflect > 100*tol)	//place to the left of inflection point
					{
						tupper = t_inflect;
						tguess = tupper-10*tol;
// 										cout<<"Compmax1b: inflect withinbounds and pos"<< endl;

					}
					else //!root is near inflection point so derivative methods will fail, run bisection
					{
						bisectFlag = true;
						tlower = t_inflect-100*tol;
						tupper = t_inflect+100*tol;
						tguess = (tlower + tupper)/2;
// 						cout<<"Compmax1:root near inflection point. Start with bisection!!"<<endl;
					}
				}
				else //!no inflection point between 0 and tExtrema
				{
					tguess = tlower;
// 									cout<<"Compmax:no inflect within bounds"<< endl;

				}
			}
			else						//tExtream is a minimum
			{
			  
				//valid region is between tExtrema=tmin and tExtrema+pi/omega
				tlower = tExtrema; //min
				tupper = tExtrema + M_PI/omegaMOD;//max
				
				//there must be an inflection point in this region
				reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
				t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
				
				FComplex(fdfval,t_inflect);
				reell f_inflect = fdfval.fval;
				
				if (f_inflect < -100*tol)
				{
					tlower = t_inflect;
					tguess = tlower+10*tol;
// 					cout<<"Compmin1a: inflect neg"<< endl;

				}
				else if (f_inflect > 100*tol)
				{
					tupper = t_inflect;
					tguess = tupper-10*tol;
// 					cout<<"Compmin1b: inflect pos"<< endl;

				}
				else //!root is near inflection point so derivative methods will fail, run bisection
				{
					bisectFlag=true;
					tlower = t_inflect-100*tol;
					tupper = t_inflect+100*tol;
					tguess = (tlower + tupper)/2;
// 					cout<<"Compmin:2root near inflection point. Start with bisection!!"<<endl;
				}
			}
		}
		else //real eigen values
		{
		
			fPoint=&LE_twoDlinear::FReal;	
			
			if (C2<C1)		//!the extrema exists and is a min
			{
				//!place tguess between inflection point and root, where function is mono-convex/concave so derivative-based methods work
				reell temp 	= lambdaMinus/lambdaPlus;
				reell tmin 	= log(   temp   * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
				if (C2<C1-10*tol) //precision fails on t_inflect when C1,C2 are too similar
				{	
					reell t_inflect	= log(temp*temp * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
					FReal(fdfval, t_inflect);
					reell f_inflect = fdfval.fval;
	// 				reell f_inflect = rootRealf(t_inflect, &rootparams);
					
	// 				tlower=(tmin>0)? tmin:0;
	// 				tupper=GSL_POSINF;
	// 				tguess=(tmin>0)? tmin+1:1;
					
					if (f_inflect < -100*tol)
					{
						tlower = (t_inflect<0)? 0:t_inflect;
						tupper = GSL_POSINF;
						tguess = tlower;
// 						cout<<"Realmin1a: inflect neg"<< endl;
						
					}
					else if (f_inflect>100*tol)
					{
						tlower = (tmin<0)? 0 : tmin;
						tupper = t_inflect;
						tguess = tupper;//-10*tol;
	// 					tguess = (tlower+tupper)/2;

// 						cout<<"Realmin1b: inflect pos"<< endl;

					}
					else //!root is near inflection point so derivative methods will fail, run bisection
					{
						bisectFlag = true;
						tlower = t_inflect-100*tol;
						tupper = t_inflect+100*tol;
						tguess = (tlower + tupper)/2;
						cout<<"Realmin1b: root near inflection point. Start with bisection!!"<<endl;
					}
				}
				else  
				{	
					tlower = (tmin<0)? 0 : tmin; 
					tupper = GSL_POSINF;
					tguess = tlower;
					//cout<<"Realmin: C1~C2!"<<endl;
				}
			}
			else if (C2 > -C1) 	//! the extrema exists and is a max
			{
				
				//! max must occur at positive times since value at t=0 negative and asymptote is positive
				//! max must be positive, so region between t=0 and root is mono-concave

				tlower = 0;
				tupper = log( lambdaMinus/lambdaPlus*(C2 - C1)/(C1 + C2) ) / (2*omegaMOD);//tmax		
				tguess = tlower;
// 				cout<<"Realmax: root before max. Start with bisection!!"<<endl;	
// 				cout << "C1 = " << C1 << endl;
// 				cout << "C2 = " << C2 << endl;
// 				cout << "C3 = " << C3 << endl;

			}
			else  			//!there is no extrema, then function monotonic so any value will do
			{
				tlower = 0;
				tupper = GSL_POSINF;
				tguess = 10*tol;
// 				cout<<"Real noExt:"<<endl;		
			}
		}
		
		//!run  algorithm
		int iter = 0, max_iter = 100;
		reell told = 0;
		reell tnew = tguess;
		cout.precision(18);

		do
		{
			//!ensure consistent bracketing
			if (tupper<tlower)
			{
				cout<< iter<<": brackets inverted!"<<endl;
				cout<< C1 << " "<<C2<<" "<<C3<<endl;
				cout<< "up:"<< tupper << " low:" << tlower <<endl;
				throw(1);
			}
			if ( (tnew < tlower) || (tnew > tupper) )
			{
				cout<< iter<<": tnew out of bracket"<<endl;
				throw(1);
			}	
			  
			iter++;
				
			//update values
// 			toldold = told;
			told = tnew;			
			(this->*fPoint)(fdfval, told);

			//update brackets
			tlower=(fdfval.fval>=0)? tlower:told;
			tupper=(fdfval.fval>=0)? told:tupper;
			
			//!if root doesn't converge by 10 iterations (usually because oscillating around precision), switch to bisection
			if (iter>30)
				bisectFlag=true;
			
			if (bisectFlag==false)
			{

				tnew = told - (fdfval.fval)/(fdfval.dfval);
				
				if ( (tnew < tlower) || (tnew > tupper) )
				{
// 					cout<<"Out of bounds: bisect in newton"<<endl;
// 					if (tupper==GSL_POSINF)
// 					{
// 						cout<<"check"<<endl;
// 						tupper=toldold;
// 						do 
// 						{tupper+=10*fabs(told-toldold);}
// 						while (rootRealf(tupper, &rootparams)<0); 
// 					}
					tnew=(tupper+tlower)/2;	
				}
			}
			else
			{
// 				if (tupper==GSL_POSINF)
// 				{
// 					cout<<"fix upper bound"<<endl;
// 					tupper=toldold;
// 					do 
// 					{tupper+=10*fabs(told-toldold);}
// 					while (rootRealf(tupper, &rootparams)<0); 
// 				}
				tnew= (tupper+tlower)/2;
			}	
// 			printf ("%5d %10.18f %+10.18f %10.18f\n", iter, tupper-tlower, tnew-told,fdfval.fval);
		}
		while ( ( (fabs(tnew-told)> tol) && (fabs(tlower-tupper)>tol) ) && (iter <= max_iter) );	//! Can also include: && (fabs(fdfval.fval)> tol) 
		
		//!the above algorithm picks out the maximum voltage if the voltage never crosses threshold, thus flag this case:
		(this->*fPoint)(fdfval, tnew);
		if (fdfval.dfval<-2*tol)
		{
			//if voltage at estimated spike time is less than threshold then remove this neuron from the list ()
			spikeTime=GSL_POSINF;
			cout<<"SPIKETIME:Max V is not within 2*pred of threshold!:"<<fdfval.dfval<<endl;
		}
		else
			spikeTime=tnew;			  

		if (!(iter < max_iter))
		{
			cout.precision(16);
			cout << " Maxed out iterations in root finding algorithm" << endl;
			cout << "state[0]" << state[0] << endl;
			cout << "state[1]" << state[1] << endl;
			cout << "C1 = " << C1 << endl;
			cout << "C2 = " << C2 << endl;
			cout << "C3 = " << C3 << endl;
			cout << "maschine precision="<<MACHINE_PRECISION << endl;
			throw(1);
		}
	}
	else if (state[0] > 1+MACHINE_PRECISION)
	{
		cout << "phase greater than 1! V = " << state[0] << " I = " << state[1] << endl;
		cout<<"Iext"<< this->Iext;
		throw(1);
	}
	else
	{
		spikeTime = 0;
	}
	
	if (spikeTime < -MACHINE_PRECISION)
	{
	    cout << "The neurons next spike time from root finding " << spikeTime << " was in the past!" << endl;
	    cout << "state[0]" << state[0] << endl;
	    cout << "state[1]" << state[1] << endl;
	    cout.precision(18);
// 	    cout << "C1=" << C1 << endl;
// 	    cout << "C2=" << C2 << endl;
// 	    cout << "C3=" << C3 << endl;
	    throw(1);
	}
	calcSpikeTime = false;
}

inline reell LE_twoDlinear::calc_resetTime()
{
	
	cout<<"calc reset runs! Works in principle, but unused, so should be checked"<<endl; throw(1);
	
	reell resetTime;
	reell C1,C2; 
	if (state[0]==0)
		resetTime=0;
	else if ( (state[0]<1) || (state[0]==1) ) //!for non reset, non suprathreshold voltages
	{
		reell C1,C2; 
		double tol = MACHINE_PRECISION;	//!>function value tolerance
		reell tguess,tlower,tupper;
		bool bisectFlag=false;
		HomemadeRtFn fPoint;
		
		C1 = state[0] - kappaV;
		C2 = (alp*(state[1] - kappaW) - (1 + realPart)*C1)/omegaMOD;

		if (state[0]>0) //!if the current state is above reset:
		{
// 			cout<<"RESET:init state between reset and thresh: move down to reset"<<endl;
			
			fdfval.C1=C1;
			fdfval.C2=C2;

			
			//!The following assigns the correct (real/complex))root function, and 
			//!compute tguess to fall in interval in which points converge to correct root and run root finfing algorithm
			
			if (iscomplex == true)
			{
				fdfval.C3=C3+stateThreshold[0]; //sets the voltage to 0.
				fPoint=&LE_twoDlinear::FComplex;
				
				//!for complex eigen values

				//!find extrema between -pi/2 and pi/2
				reell temp1=omegaMOD*C2 + realPart*C1;
				reell temp2=realPart*C2 - omegaMOD*C1;
				reell tExtrema = atan( -temp1/temp2 )/omegaMOD;
				
				//!if tExtrema postive then look at next extrema time (which is necessarily negative)
				tExtrema -= (tExtrema>0)? M_PI/omegaMOD : 0;
				
				if ((C2 > 0) || (realPart*C1 > -omegaMOD*C2))		// tExtrema is a minimum,  shouldn't happen because in spiektimecase
				{
				  
					//inflection point
					reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
			
					(this->*fPoint)(fdfval, tExtrema);
					if (fdfval.fval<0)
					{
						//valid region is between tExtrema=tmin and tExtrema+pi/omega
						tlower= tExtrema; //min
						tupper= 0;//max
						
						//pick smallest negative one
						t_inflect -= (t_inflect > 0)? M_PI/omegaMOD : 0;
						//is first negative inflection point in this region?
					}
					else
					{
						tlower= tExtrema-2*M_PI/omegaMOD; //min
						tupper= tExtrema-M_PI/omegaMOD;//max

						t_inflect -= (t_inflect > 0)? M_PI/omegaMOD : 0;
						t_inflect -= M_PI/omegaMOD;
					}
					
					if ( (tlower < t_inflect) && (t_inflect < tupper) )
					{
						(this->*fPoint)(fdfval, t_inflect);
						reell f_inflect=fdfval.fval;
						
						if (f_inflect < -100*tol) 	//place to the right of inflection point
						{
							tlower = t_inflect;
							tguess = tlower+10*tol;
							cout<<"RESETCompmax1a: inflect withinbounds and neg"<< endl;
						}
						else if (f_inflect > 100*tol)	//place to the left of inflection point
						{
							tupper = t_inflect;
							tguess = tupper-10*tol;
							cout<<"RESETCompmax1b: inflect withinbounds and pos"<< endl;
						}
						else //root is near inflection point so derivative methods will fail, run bisection
						{
							bisectFlag = true;
							tlower = -100*tol;
							tupper = 100*tol;
							tguess = (tlower + tupper)/2;
							cout<<"RESETCompmax1:root near inflection point. Start with bisection!!"<<endl;
						}
					}
					else
						tguess = tupper;
				}
				else							//tExtream is a minimum
				{
					cout << "RESET:extrema is a minimum!" << endl;
					//valid region is between tExtrema=tmin and tExtrema+pi/omega
					tlower = tExtrema; //min is empirically always seen to be less than reset
					tupper = 0; //treset must be negative
					
					//there must be an inflection point in this region
					reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
					t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
					
					//is there must be an inflection point in this region?
					if ( (t_inflect > tlower) && (t_inflect < tupper))
					{
						//is it above or below reset votlage?
						(this->*fPoint)(fdfval, t_inflect);
						reell f_inflect=fdfval.fval;
					
						if (f_inflect < -100*tol)
						{
							tlower = t_inflect;
							tguess = tlower+10*tol;
							cout<<"RESET:Compmin1a: inflect withinbounds and neg"<< endl;
						}
						else if (f_inflect > 100*tol)
						{
							tupper = t_inflect;
							tguess = tupper-10*tol;
							cout<<"RESET:Compmin1a: inflect withinbounds and pos"<< endl;

						}
						else //root is near inflection point so derivative methods will fail, run bisection
						{
							bisectFlag = true;
							tlower = -100*tol;
							tupper = 100*tol;
							tguess = (tlower + tupper)/2;
							cout<<"RESET:Compmin1:root near inflection point. Start with bisection!!"<<endl;

						}
					}
					else
					{
						tguess = tlower+tol;
						cout<<"RESET:Compmax:no inflect within bounds"<< endl;
					}
				}
			}
			else //real eigen values
			{
			
				fPoint=&LE_twoDlinear::FReal;	
				fdfval.C3=C3+2*stateThreshold[0];
				reell temp1 = omegaMOD*C2 + realPart*C1;
				reell temp2 = realPart*C2 - omegaMOD*C1;			
				if (C2<C1) 		//!the extrema exists and is a min, min mustb e less than 0.
				{
					reell temp = lambdaMinus/lambdaPlus;
					reell tmin 	= log(   temp   * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
					if (C2<C1-10*tol) //precision fails on t_inflect when C1,C2 are too similar
					{
						//!place tguess between inflection point and root, where function is mono-convex/concave so derivative-based methods work					reell t_inflect	= log(temp*temp * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
						
						reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
						t_inflect -= (t_inflect > 0)? M_PI/omegaMOD : 0;

						(this->*fPoint)(fdfval, t_inflect);
						reell f_inflect = fdfval.fval;
						
						if (f_inflect < -100*tol)
						{
							tlower = t_inflect;//potentially creates problems when swithcing to bisection
							tupper = 0;
							tguess = tlower+10*tol;
							cout<<"RESET:Realmin1a: inflect neg"<< endl;
						}
						else if (f_inflect > 100*tol)
						{
							tupper = (t_inflect>0)? 0: t_inflect;
							tlower = tmin;
							tguess = tupper-10*tol;
							cout<<"RESET:Realmin1b: inflect pos"<< endl;
						}
						else //!root is near inflection point so derivative methods will fail, run bisection
						{
							bisectFlag = true;
							tlower = t_inflect-100*tol;
							tupper = t_inflect+100*tol;
							tguess = (tlower + tupper)/2;
							cout<<"RESET:Realmin1b: root near inflection point. Start with bisection!!"<<endl;
						}
					}
					else //for C1 ~ C2
					{
						tlower =  GSL_NEGINF;
						tupper = (tmin>0)? 0 : tmin;
						tguess = tupper;
						cout<<"RESET:Realmin: C1~C2!"<<endl;
					}
				}
				else if (C2 > -C1) 	//! the extrema exists and is a max
				{
					
					//! max must occur at positive times as in calcspiketime
					//! max must be positive, so region between t=0 and root is mono-concave

					tlower = GSL_NEGINF;
					tupper = 0;
					tguess = tupper-tol;
					cout<<"RESET:Realmax: root before max. Start with bisection!!"<<endl;		

				}
				else  			//!there is no extrema, then function monotonic so any value will do
				{
					tlower = GSL_NEGINF;
					tupper = 0;					
					tguess = tupper-1;
// 					cout<<"RESET:Real noExt:"<<endl;
				}
			}
			
			
		}
		else //if (state[0] <0) //pick init based on moving up to Vreset (ie. like orig spike time algorithm)
		{
			cout<<"RESET:init state below reset"<<endl;
			
			fdfval.C1=C1;
			fdfval.C2=C2;
			
			//!The following assigns the correct (real/complex))root function, and 
			//!computes tguess to fall in interval in which points converge to correct root and run root finfing algorithm
			
			if (iscomplex == true)
			{
				fPoint = &LE_twoDlinear::FComplex;
				fdfval.C3=C3+stateThreshold[0]; //sets teh voltage to 0.
				//!for complex eigen values

				//!find extrema between -pi/2 and pi/2
				reell temp1 = omegaMOD*C2 + realPart*C1;
				reell temp2 = realPart*C2 - omegaMOD*C1;
				reell tExtrema = atan( -temp1/temp2 )/omegaMOD;
				
				//!if tExtrema negative then look at next extrema time (which is necessarily positive)
				tExtrema += (tExtrema<0)? M_PI/omegaMOD : 0;
				
				if ((C2 > 0) || (realPart*C1 > -omegaMOD*C2))	//tExtrema is a maximum, 
				{
					tlower = 0;
					tupper = tExtrema;
					
					//is first positive inflection point in this region?
					reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
					t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
					
					if ( (tlower < t_inflect) && (t_inflect < tupper) )
					{
						FComplex(fdfval,t_inflect);
						reell f_inflect = fdfval.fval;
						
						if (f_inflect < -100*tol) 	//place to the right of inflection point
						{
							tlower = t_inflect;
							tguess = tlower+10*tol;
											cout<<"RESET:Compmax1a: inflect withinbounds and neg"<< endl;

						}
						else if (f_inflect > 100*tol)	//place to the left of inflection point
						{
							tupper = t_inflect;
							tguess = tupper-10*tol;
											cout<<"RESET:Compmax1b: inflect withinbounds and pos"<< endl;

						}
						else //!root is near inflection point so derivative methods will fail, run bisection
						{
							bisectFlag = true;
							tlower = t_inflect-100*tol;
							tupper = t_inflect+100*tol;
							tguess = (tlower + tupper)/2;
							cout<<"RESET:Compmax1:root near inflection point. Start with bisection!!"<<endl;
						}
					}
					else //no inflection point between 0 and tExtrema
					{
						tguess = tlower;
	// 									cout<<"Compmax:no inflect within bounds"<< endl;

					}
				}
				else						//tExtream is a minimum
				{
				  
					//valid region is between tExtrema=tmin and tExtrema+pi/omega
					tlower = tExtrema; //min
					tupper = tExtrema + M_PI/omegaMOD;//max
					
					//there must be an inflection point in this region
					reell t_inflect = atan( (realPart*temp1 + omegaMOD*temp2) / (omegaMOD*temp1 - realPart*temp2) )/omegaMOD;
					t_inflect += (t_inflect<0)? M_PI/omegaMOD : 0;
					
					FComplex(fdfval,t_inflect);
					reell f_inflect = fdfval.fval;
					
					if (f_inflect < -100*tol)
					{
						tlower = t_inflect;
						tguess = tlower+10*tol;
						cout<<"RESET:Compmin1a: inflect  neg"<< endl;

					}
					else if (f_inflect > 100*tol)
					{
						tupper = t_inflect;
						tguess = tupper-10*tol;
						cout<<"RESET:Compmin1b: inflect pos"<< endl;

					}
					else //!root is near inflection point so derivative methods will fail, run bisection
					{
						bisectFlag=true;
						tlower = t_inflect-100*tol;
						tupper = t_inflect+100*tol;
						tguess = (tlower + tupper)/2;
						cout<<"RESET:Compmin:2root near inflection point. Start with bisection!!"<<endl;
					}
				}
			}
			else //real eigen values
			{
			
				fPoint=&LE_twoDlinear::FReal;	
				fdfval.C3=C3+2*stateThreshold[0]; //sets teh voltage to 0.
				if (C2<C1)		//!the extrema exists and is a min
				{
					//!place tguess between inflection point and root, where function is mono-convex/concave so derivative-based methods work
					reell temp 	= lambdaMinus/lambdaPlus;
					reell tmin 	= log(   temp   * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
					if (C2<C1-10*tol) //precision fails on t_inflect when C1,C2 are too similar
					{	
						reell t_inflect	= log(temp*temp * (C2 - C1)/(C1 + C2) ) / (2*omegaMOD);
						FReal(fdfval, t_inflect);
						reell f_inflect = fdfval.fval;
		// 				reell f_inflect = rootRealf(t_inflect, &rootparams);
						
		// 				tlower=(tmin>0)? tmin:0;
		// 				tupper=GSL_POSINF;
		// 				tguess=(tmin>0)? tmin+1:1;
						
						if (f_inflect < -100*tol)
						{
							tlower = (t_inflect<0)? 0:t_inflect;
							tupper = GSL_POSINF;
							tguess = tlower;
							cout<<"RESET:Realmin1a: inflect neg"<< endl;
							
						}
						else if (f_inflect>100*tol)
						{
							tlower = (tmin<0)? 0 : tmin;
							tupper = t_inflect;
							tguess = tupper;//-10*tol;
		// 					tguess = (tlower+tupper)/2;

							cout<<"RESET:Realmin1b: inflect pos"<< endl;

						}
						else //!root is near inflection point so derivative methods will fail, run bisection
						{
							bisectFlag = true;
							tlower = t_inflect-100*tol;
							tupper = t_inflect+100*tol;
							tguess = (tlower + tupper)/2;
							cout<<"RESET:Realmin1b: root near inflection point. Start with bisection!!"<<endl;
						}
					}
					else  
					{	
						tlower = (tmin<0)? 0 : tmin; 
						tupper = GSL_POSINF;
						tguess = tlower;
						cout<<"RESET:Realmin: C1~C2!"<<endl;
					}
				}
				else if (C2 > -C1) 	//! the extrema exists and is a max
				{
					
					//! max must occur at positive times as in calcspiketime
					//! max must be positive, so region between t=0 and root is mono-concave
					tlower = 0;
					tupper = log( lambdaMinus/lambdaPlus*(C2 - C1)/(C1 + C2) ) / (2*omegaMOD);//tmax		
					tguess = tlower+tol;
					cout<<"RESET:Realmax: root before max. Start with bisection!!"<<endl;		

				}
				else  			//!there is no extrema, then function monotonic so any value will do
				{
					tlower = 0;
					tupper = GSL_POSINF;
					tguess = 10*tol;
					cout<<"RESET:Real noExt:"<<endl;		
				}
			}
		}
		
		//!run  algorithm
		int iter = 0, max_iter = 100;
		reell told = 0;
		reell tnew = tguess;
		//!run through algorithm
		cout.precision(18);

		do
		{

			//!ensure consistent bracketing
			if (tupper<tlower)
			{
				cout<< iter<<"RESET:: brackets inverted!"<<endl;
				cout<< C1 << " "<<C2<<" "<<C3<<endl;
				cout<< "up:"<< tupper << " low:" << tlower <<endl;
				throw(1);
			}
			if ( (tnew < tlower) || (tnew > tupper) )
			{
				cout<< iter<<"RESET:: tnew out of bracket"<<endl;
				throw(1);
			}	
			
			  
			iter++;
				
			told = tnew;
			(this->*fPoint)(fdfval, told);
			
			//update brackets
			tlower=(fdfval.fval>=0)? tlower:told;
			tupper=(fdfval.fval>=0)? told:tupper;
			
			//!if root doesn't converge by 10 iterations (usually because oscilalting aroudn precision), switch to bisection
			if (iter>30)
				bisectFlag=true;
			
			//!homemade newton:
			if (bisectFlag==false)
			{
				tnew = told - (fdfval.fval)/(fdfval.dfval);
				
				//! apply bisection if tnew goes out of bounds.
				if ( (tnew < tlower) || (tnew > tupper) )
				{
	// 					cout<<"bisect in newton"<<endl;
	// 					if (tupper==GSL_POSINF)
	// 					{
	// 						cout<<"check"<<endl;
	// 						tupper=toldold;
	// 						do 
	// 						{tupper+=10*fabs(told-toldold);}
	// 						while (rootRealf(tupper, &rootparams)<0); 
	// 					}
					tnew=(tupper+tlower)/2;	
				}
			}
			else
			{
	// 				cout<<"bisect after 10 iterations"<<endl;
	// 				if (tupper==GSL_POSINF)
	// 				{
	// // 					cout<<"fix upper bound"<<endl;
	// 					tupper=toldold;
	// 					do 
	// 					{tupper+=10*fabs(told-toldold);}
	// 					while (rootRealf(tupper, &rootparams)<0); 
	// 				}
					
				
				//! apply bisection
				tnew= (tupper+tlower)/2;
			}	
				//printf ("%5d %10.18f %+10.18f %10.18f\n", iter, tupper-tlower, tnew-told,fdfval.fval);
		}	
		while ( ( (fabs(tnew-told)> tol) && (fabs(tlower-tupper)>tol) ) && (iter <= max_iter) );//&& (oscflag==0) );	 fabs(tupper-tlower)>tol) && 						//! homemade conditions && (fabs(fdfval.fval)> tol) 										//! homemade estiamte
		
		//the above algorithm picks out the minimum/maximum voltage if the voltage never crosses reset, thus flag this case:
		(this->*fPoint)(fdfval, tnew);
		if (fabs(fdfval.dfval)<10*tol)
			resetTime=tnew;	
		else
		{	
			cout<<"RESET: Max V is not within 10*prec of threshold!:"<<fdfval.dfval<< " phase velocity not defined!"<<endl;
			cout<<"difft:"<<fabs(tnew-told)<<endl;
			resetTime=tnew;
		}
		
		if (!(iter < max_iter))
		{
			cout.precision(16);
			cout << " RESET:Maxed out iterations in root finding algorithm" << endl;
			cout << "state[0]" << state[0] << endl;
			cout << "state[1]" << state[1] << endl;
			cout << "C1 = " << C1 << endl;
			cout << "C2 = " << C2 << endl;
			cout << "C3 = " << C3 << endl;
			cout << "maschine precision="<<MACHINE_PRECISION << endl;
			throw(1);
		}
	}
	else
	{
		cout << "RESET:phase greater than 1! V = " << state[0] << " I = " << state[1] << endl;
		cout<<"Iext"<< this->Iext;
		cout << "C1 = " << C1 << endl;
		cout << "C2 = " << C2 << endl;
		cout << "C3 = " << C3 << endl;
		throw(1);	
	}
	
	return resetTime;
}

void LE_twoDlinear::dummy_evolve(reell dt, vector<reell>* dummy)
{
	//!the following calculates the propogation matrix,S, 
	//!it is only run for one (the 1st) neuron, but S is valid for all neurons (even with heterogeneous taus)
	//!The dummy variable is used in evolve_dt for all neurons.
	
	dt /= tauM;
	
	
	(*dummy).resize(4);

	reell S[4];

	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		S[0] = expRdt*( cosOmega + sinOmega*(A[0]-realPart)/omegaMOD );
		S[1] = expRdt*sinOmega*A[1]/omegaMOD;
		S[2] = expRdt*sinOmega*A[2]/omegaMOD;
		S[3] = expRdt*( cosOmega + sinOmega*(A[3]-realPart)/omegaMOD );
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		S[0] = coeff*A[0] + tempS;
		S[1] = coeff*A[1];
		S[2] = coeff*A[2];
		S[3] = coeff*A[3] + tempS;
	}

	//!assign to dummy
	(*dummy)[0] = S[0];
	(*dummy)[1] = S[1];
	(*dummy)[2] = S[2];
	(*dummy)[3] = S[3];

}

void LE_twoDlinear::evolve_dt(reell dt, vector<reell>* dummy)
{
	//! increments time between spikes in the network
	jacISI += dt;
	
	//!evolve neurons by dt. dummy is the propogration matrix for this time.
	reell v_temp=state[0];
	state[0] = (*dummy)[0]*(state[0] - kappaV) + (*dummy)[1]*(state[1]-kappaW)+kappaV;
	state[1]  = (*dummy)[2]*(v_temp - kappaV)  + (*dummy)[3]*(state[1]-kappaW)+kappaW;

	//!adjust the neuron's spike time
	calcSpikeTime = true;
	
	//! We cannot just use spikeTime -= dt because this acccumulates precision losses. For a demonstration run a distance measurement
	//! of a zero vector with spikeTime -= dt and you will see how the distance increases.
	

}

void LE_twoDlinear::evolve_spike(reell c)
{
	//! \a state[0] - the voltage and the current \a state[1] can change
	state[0] +=     gam*c;
	state[1] += deltauS*c;

	calcSpikeTime = true;
}

//! Set the state variable to the reset value.
void LE_twoDlinear::reset()
{
	//!voltage and current reset
	state[0]  = stateReset[0];
	state[1]  += stateReset[1];//reset value not clear, for resonator with slow currents, no rest is fine

// 	state[1] = -2*(Iext+Irheo);//reset value not clear, for resonator with slow currents, no rest is fine

	calcSpikeTime = true;
}

void LE_twoDlinear::dummy_PostUpdate()
{
	//! This is saved for each neuron, but only to be used in dummy_Jacobian_spiking
	initSpikingV = state[0];
	initSpikingW = state[1];
		
	//! reset at network spike time
	jacISI=0;  
}

void LE_twoDlinear::dummy_Jacobian_spiking(vector<reell>* dummy)
{
	//! the following calculates the propogation matrix and dtau/dz. It is run for the spikign neuron.
	(*dummy).resize(6);

	reell dt = jacISI/tauM; //convert to units of spiking neuron's 'membrane time constant	
	
	reell dS11dt;
	reell dS12dt;
	reell S[4];
	
	if (iscomplex==true) 
	{
		//!some temp variables
		reell argO	= omegaMOD*dt;
		reell expRdt	= exp(realPart*dt);
		reell cosOmega	= cos(argO);
		reell sinOmega  = sin(argO);
		
		//!the propogation matrix
		S[0] = expRdt*( cosOmega + sinOmega*(A[0] - realPart)/omegaMOD );
		S[1] = expRdt*sinOmega*A[1]/omegaMOD;
		S[2] = expRdt*sinOmega*A[2]/omegaMOD;
		S[3] = expRdt*( cosOmega + sinOmega*(A[3] - realPart)/omegaMOD );
		
		dS11dt = realPart*S[0] - expRdt*(omegaMOD*sinOmega - cosOmega*(A[0] - realPart));
		dS12dt = realPart*S[1] + expRdt*cosOmega*A[1];
	}
	else
	{	
		//!some temp variables
		reell expplus	= exp(lambdaPlus*dt);
		reell expminus  = exp(lambdaMinus*dt);
		reell coeff	= (expplus - expminus)/(2*omegaMOD);
		reell tempS	= expplus - coeff*lambdaPlus;
		
		//!the propogation matrix
		S[0] = coeff*A[0] + tempS;
		S[1] = coeff*A[1];
		S[2] = coeff*A[2];
		S[3] = coeff*A[3] + tempS;
		
		reell lambexpplus = lambdaPlus*expplus;
		reell lambexpminus = lambdaMinus*expminus;

		dS11dt = lambexpplus + (lambexpplus - lambexpminus)/(2*omegaMOD)*(A[0] - lambdaPlus);
		dS12dt = (lambexpplus - lambexpminus)*A[1]/(2*omegaMOD);
	}
	

// 	dS11dt /= tauM; //convert back to milliseconds as if it was the spiking neuron. assumes gamma same for all.
// 	dS12dt /= tauM;


	(*dummy)[0] = S[0];
	(*dummy)[1] = S[1];
	(*dummy)[2] = S[2];
	(*dummy)[3] = S[3];
	
	//dtaudz in temporal units of milliseconds
	(*dummy)[4] = S[0]/( -dS11dt/S[0]*(stateThreshold[0] - kappaV - S[1]*(initSpikingW - kappaW)) - dS12dt*(initSpikingW - kappaW) );
	(*dummy)[5] = S[1]/( -dS12dt/S[1]*(stateThreshold[0] - kappaV - S[0]*(initSpikingV - kappaV)) - dS11dt*(initSpikingV - kappaV) );
}

//! The diagonal elements of the single spike Jacobian for the postsynaptic neurons.
vector<vector<reell> > LE_twoDlinear::calc_JacElem_postsynaptic(vector<reell>* dummy, reell c)
{

	vector<vector<reell> > jacTmp(2, vector<reell>(4));				//!< the vector with Jacobian elements														
	//! 4 elements for the diagonal and the nondiagonal part, each.
	jacTmp[0][0] = (*dummy)[0];	
	jacTmp[0][1] = (*dummy)[1];
	jacTmp[0][2] = (*dummy)[2];		
	jacTmp[0][3] = (*dummy)[3];
	
// 	(*dummy)[4]/=tauM; //convert back to units of membrane time constant. This time of the post synaptic neuron.
// 	(*dummy)[5]/=tauM;
	
	jacTmp[1][0] = c*dzdt[2]*(*dummy)[4];															
	jacTmp[1][1] = c*dzdt[2]*(*dummy)[5];																						
	jacTmp[1][2] = c*dzdt[3]*(*dummy)[4];																		
	jacTmp[1][3] = c*dzdt[3]*(*dummy)[5];								
	return jacTmp;
}

vector<vector<reell> > LE_twoDlinear::calc_JacElem_self(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));
	//!the diagonal elements are just the propagation matrix, S
	jacTmp[0][0] = (*dummy)[0];											
	jacTmp[0][1] = (*dummy)[1];											
	jacTmp[0][2] = (*dummy)[2];											
	jacTmp[0][3] = (*dummy)[3];
	return jacTmp;
}

vector<vector<reell> > LE_twoDlinear::calc_JacElem_spiking(vector<reell>* dummy)
{
	vector<vector<reell> > jacTmp(1, vector<reell>(4));
	
// 	(*dummy)[4]/=tauM; //convert back to units of membrane time constant. This time of the spiking neuron.
// 	(*dummy)[5]/=tauM; 
	
	jacTmp[0][0] = (*dummy)[0] + dzdt[0]*(*dummy)[4];
	jacTmp[0][1] = (*dummy)[1] + dzdt[0]*(*dummy)[5];
	jacTmp[0][2] = (*dummy)[2] + dzdt[1]*(*dummy)[4];
	jacTmp[0][3] = (*dummy)[3] + dzdt[1]*(*dummy)[5];
	return jacTmp;
}



void LE_twoDlinear::FReal(LE_twoDlinear::st_ftuple &fdfval, reell tval)
{

  	//! Some temporary variables
	reell expPlus  = exp(lambdaPlus*tval);
	reell expMinus = exp(lambdaMinus*tval);
	
	//! The function value and the derivative.
	fdfval.fval  = fdfval.C1*(expPlus + expMinus) + fdfval.C2*(expPlus - expMinus) + fdfval.C3;
	fdfval.dfval  = fdfval.C1*(lambdaPlus*expPlus + lambdaMinus*expMinus) + fdfval.C2*(lambdaPlus*expPlus - lambdaMinus*expMinus);

}

void LE_twoDlinear::FComplex(LE_twoDlinear::st_ftuple &fdfval, reell tval)
{
	//! Some temporary variables
	reell argO 		    = omegaMOD*tval;
	reell cosOmegatemp 	= cos(argO);
	reell sinOmegatemp 	= sin(argO);
	reell expRdttemp 	= exp(realPart*tval);
	reell temp 		    = expRdttemp*(fdfval.C1*cosOmegatemp + fdfval.C2*sinOmegatemp);
	
	//! The function value and the derivative.
	fdfval.fval = temp + fdfval.C3;
	fdfval.dfval = realPart*temp - omegaMOD*expRdttemp*(fdfval.C1*sinOmegatemp - fdfval.C2*cosOmegatemp);
} 