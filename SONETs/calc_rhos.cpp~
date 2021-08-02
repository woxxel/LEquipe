#include<iostream>
#include<stdio.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_roots.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_integration.h>
#include "calc_rhos.hpp"

using namespace std;


// function and structure declarations
double integrate_gaussian(double rho0, void *parameters);
struct Qfunction_params { double xth, rho; };
struct int_gauss_params{ double p1, p2, sec_mom; };




////////////////////////////////////////////////////
// calc_rho_given_alpha
// need to calculate value of rho, the correlation
// between two Gaussian random variables
// so that the covariance of the resulting Bernoulli
// variables will be alpha*p1*p2
// solve this 1D equation for rho numerically
// (The variance of each Gaussian is determined by 
// the requirement that the Bernoulli obtained by
// thresholding at 1 has probability pj of being 1.)
/////////////////////////////////////////////////////

double calc_rho_given_alpha(double p1, double p2, double alpha,
			    int &status) {

  // if correlation or a probability is zero, return rho=0
  if(alpha==0 || p1==0 || p2==0) {
    status = 0;
    return 0.0;
  }
  if(alpha < -1) {
    cerr << "alpha < -1, cannot calc rho" << endl;
    status = -1;
    return 0.0;
  }
  
  // desired second moment among bernoulli random variables
  double b_sec_mom = p1*p2*(1+alpha);

  // set up gsl function FF to point to function
  // integrate_gaussian with
  // parameters determined by arguments p and bcorr
  struct int_gauss_params para;
  para.p1 = p1;
  para.p2 = p2;
  para.sec_mom = b_sec_mom;
  gsl_function FF;		//! defines FF(double x, void * parameters), that has members "function" and "params"
  FF.function = &integrate_gaussian;	//! hand function to FF
  FF.params = &para;			//! hand parameters to FF



  int iter=0, max_iter=1000;
  double rho=0;
  
  // we know rho has to be in interval [x_lo,x_hi]
  double x_lo = -1, x_hi=1;
  if (alpha > 0)
   x_lo=0;
  else x_hi=0;

  // Initialize solver and iterate until converge to solution
  const gsl_root_fsolver_type *T=gsl_root_fsolver_brent;	//! Brent = quick root finding algorithm
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &FF, x_lo, x_hi);	//! initialize "solver" s for function FF in search intervall [x_lo,x_hi]

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);	//! update state of solver
    if(status)
      break;

    rho = gsl_root_fsolver_root (s);		//! current (best) estimate of root
    x_lo = gsl_root_fsolver_x_lower (s);	//! current lower and upper boundary estimate
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, 0.00001);	//! test if solution converged (absolute-, relative error)
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);

  return rho;

}


// Q function, when integrated from xth to infinity,
// gives the second moment of the Bernoulli random variables
double Qfunction(double x, void *parameters){

  struct Qfunction_params *para = (struct Qfunction_params *)parameters;
  double rho = para->rho;
  double xth = para->xth;
  double Q=gsl_cdf_gaussian_Q(xth-rho*x, sqrt(1-rho*rho));	//! 1-cumulative gaussian with variance sqrt(1-rho^2), mu = 0
  double g= 1/sqrt(2*M_PI)*exp(-0.5*x*x)*Q;
  return g;

}


// Take the integral of the Gaussian that corresponds to the
// second moment of the Bernoulli random variables.
// Subtract off their required value so that function is 
// zero when found correct value of rho
double integrate_gaussian(double rho0, void *parameters){

  struct int_gauss_params * para = (struct int_gauss_params *) parameters;
  double p1=para->p1;
  double p2=para->p2;
  double sec_mom=para->sec_mom;
  double rho=rho0;
  double xth1=gsl_cdf_ugaussian_Qinv(p1);	//! inverse(1-(cumulative gaussian with variance 1, mu = 0)) phi^(-1)(p1)
  double xth2=gsl_cdf_ugaussian_Qinv(p2);

  struct Qfunction_params qparams;
  qparams.xth=xth1;
  qparams.rho=rho;
  gsl_function G;
  G.function = &Qfunction;
  G.params = &qparams;
  
  //! integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(1000);	//! workspace for storage of integration data (up to 1000 intervalls)
  double res, err;
  gsl_integration_qagiu(&G, xth2, 1e-7, 1e-7, 1000, ww, &res, &err);		//! integrate G from xth2 to infty
  gsl_integration_workspace_free (ww);

  res -= sec_mom;
  return res;

}
