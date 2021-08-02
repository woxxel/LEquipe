#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

gsl_matrix* secorder_rec_1p(int N_nodes, double p, 
			    double alpha_recip, double alpha_conv, 
			    double alpha_div, double alpha_chain,
			    gsl_rng *rng);
