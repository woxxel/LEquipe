#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include "secorder_rec_1p.hpp"
#include "calc_stats_1p.hpp"


using namespace std;



int main(int argc, char *argv[]) {
  
  // optional flags
  int deterministic_seed = 0; // if nonzero, use for random num gen seed

  ///////////////////////////////////////////////////////////////
  // Read fourteen input parameters 
  // N_nodes      number of nodes in the network
  // p            the probability of any one given connection
  // alpha_recip_{min/max/n}
  // alpha_conv_{min/max/n}
  // alpha_div_{min/max/n}
  // alpha_chain_{min/max/n}
  ///////////////////////////////////////////////////////////////


  if(argc !=15) {
    cerr << "Requires eight parameters: N_nodes p alpha_recip_{min/max/n} alpha_conv_{min/max/n} alpha_div_{min/max/n} alpha_chain_{min/max/n} \n";
    exit(-1);
  }

  int N_nodes = atoi(argv[1]);
  double p = atof(argv[2]);
  double alpha_recip_min=atof(argv[3]);
  double alpha_recip_max=atof(argv[4]);
  int alpha_recip_n=atoi(argv[5]);
  
  double alpha_conv_min=atof(argv[6]);
  double alpha_conv_max=atof(argv[7]);
  int alpha_conv_n=atoi(argv[8]);

  double alpha_div_min=atof(argv[9]);
  double alpha_div_max=atof(argv[10]);
  int alpha_div_n=atoi(argv[11]);

  double alpha_chain_min=atof(argv[12]);
  double alpha_chain_max=atof(argv[13]);
  int alpha_chain_n=atoi(argv[14]);

  

  cout << "Generating " << alpha_recip_n*alpha_conv_n*alpha_div_n*alpha_chain_n << " networks with\n"
       << "alpha_recip in [" << alpha_recip_min << ", " << alpha_recip_max << "]\n"
       << "alpha_conv in [" << alpha_conv_min << ", " << alpha_conv_max << "]\n"
       << "alpha_div in [" << alpha_div_min << ", " << alpha_div_max << "]\n"
       << "alpha_chain in [" << alpha_chain_min << ", " << alpha_chain_max << "]\n";

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

  // set the seed
  // if deterministic_seed, use that for seed
  if(deterministic_seed)
    gsl_rng_set(r,deterministic_seed); 
  // else use time in seconds for the seed
  else
    gsl_rng_set(r, time(NULL));


  double d_alpha_recip=(alpha_recip_max-alpha_recip_min)/alpha_recip_n;
  double d_alpha_conv=(alpha_conv_max-alpha_conv_min)/alpha_conv_n;
  double d_alpha_div=(alpha_div_max-alpha_div_min)/alpha_div_n;
  double d_alpha_chain=(alpha_chain_max-alpha_chain_min)/alpha_chain_n;

  for(int i_recip=0; i_recip < alpha_recip_n; i_recip++) 
    for(int i_conv=0; i_conv < alpha_conv_n; i_conv++) 
      for(int i_div=0; i_div < alpha_div_n; i_div++)
	for(int i_chain=0; i_chain < alpha_chain_n; i_chain++) {

	  // randomly generated alpha values Latin Hypercube intervals
	  cout << "\nLatin hypercube: [" 
	       << alpha_recip_min+i_recip*d_alpha_recip 
	       << ", " << alpha_recip_min+(i_recip+1)*d_alpha_recip
	       << "] x ["
	       << alpha_conv_min+i_conv*d_alpha_conv 
	       << ", " << alpha_conv_min+(i_conv+1)*d_alpha_conv
	       << "] x ["
	       << alpha_div_min+i_div*d_alpha_div 
	       << ", " << alpha_div_min+(i_div+1)*d_alpha_div
	       << "] x ["
	       << alpha_chain_min+i_chain*d_alpha_chain 
	       << ", " << alpha_chain_min+(i_chain+1)*d_alpha_chain
	       << "]\n";

	  // try multiple times to get a matrix in this hypercube
	  int N_tries = 100;

	  for(int j=0; j< N_tries; j++) {
	    double alpha_recip = gsl_ran_flat(r, alpha_recip_min+i_recip*d_alpha_recip,
					      alpha_recip_min+(i_recip+1)*d_alpha_recip);
	    double alpha_conv = gsl_ran_flat(r, alpha_conv_min+i_conv*d_alpha_conv,
					     alpha_conv_min+(i_conv+1)*d_alpha_conv);
	    double alpha_div = gsl_ran_flat(r, alpha_div_min+i_div*d_alpha_div,
					    alpha_div_min+(i_div+1)*d_alpha_div);
	    double alpha_chain = gsl_ran_flat(r, alpha_chain_min+i_chain*d_alpha_chain,
					     alpha_chain_min+(i_chain+1)*d_alpha_chain);

	    if(alpha_chain <= -3) {
	      cout << "Using min alpha_chain\n";
	    }
	    else if(alpha_chain <= -2) {
	      cout << "Using max alpha_chain\n";
	    }
	    
	    
	    gsl_matrix *W;

	    W=secorder_rec_1p(N_nodes,p, alpha_recip, alpha_conv, 
			      alpha_div, alpha_chain, r);
	    
	    if(W) {

	      // matrix files will be stored in data directory
	      // with filename determined by the six input parameters
	      mkdir("data",0755);  // make data directory, if it doesn't exist
	      char FNbase[200];
	      sprintf(FNbase,"_%i_%1.3f_%1.3f_%1.3f_%1.3f_%1.3f",N_nodes,p,alpha_recip, alpha_conv, alpha_div, alpha_chain);
	      
	      char FN[200];
	      FILE *fhnd;
	      strcpy(FN, "data/w");
	      strcat(FN, FNbase);
	      strcat(FN, ".dat");
	      fhnd = fopen(FN, "w");
	      if(fhnd==NULL) {
		cerr << "Couldn't open outfile file " << FN << "\n";
		exit(-1);
	      }
	      
	      for(int i=0; i<N_nodes;i++) {
		for(int j=0; j<N_nodes; j++) {
		  fprintf(fhnd, "%i ", (int) gsl_matrix_get(W,i,j));
		}
		fprintf(fhnd,"\n");
	      }
	      fclose(fhnd);
	      
	      cout << "done\n";
	      
	      
	      ////////////////////////////////////////////////////////////
	      // Calculate the covariance structure the adjacency matrix
	      // This should approximately agree with the input alphas
	      ////////////////////////////////////////////////////////////
	      
	      cout << "Testing statistics of W ...";
	      cout.flush();
	      double phat, alphahat_recip, alphahat_conv, alphahat_div,
		alphahat_chain, alphahat_other;
	      
	      calc_phat_alphahat_1p(W, N_nodes, phat, alphahat_recip, 
				    alphahat_conv, alphahat_div,
				    alphahat_chain, alphahat_other);
	      
	      strcpy(FN, "data/stats");
	      strcat(FN, FNbase);
	      strcat(FN, ".dat");
	      fhnd = fopen(FN, "w");
	      if(fhnd==NULL) {
		cerr << "Couldn't open outfile file " << FN << "\n";
		exit(-1);
	      }
	      fprintf(fhnd, "%e %e %e %e %e %e\n", phat, alphahat_recip, 
		      alphahat_conv, alphahat_div, alphahat_chain, alphahat_other);
	      fclose(fhnd);
	      
	      
	      cout << "done\n";
	      cout << "\nActual statistics of matrix:\n";
	      cout << "phat = " << phat << "\n";
	      cout << "alphahats:\n";
	      cout << "alpha_recip = " << alphahat_recip
		   << ", alpha_conv = " << alphahat_conv
		   << ", alpha_div = " << alphahat_div
		   << ", alpha_chain = " << alphahat_chain
		   << ", alpha_other = " << alphahat_other << "\n";
	      
	      

	      gsl_matrix_free(W);
	      
	      // since was successful, 
	      // break out of loop attempting for this latin hypercube
	      break;
	    }
	  }
	}

  
  return 0;

}
