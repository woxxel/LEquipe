#include <gsl/gsl_matrix.h>
#include <vector>

#include "structures.hpp"

using namespace std;

int calc_phat_alphahat_1p(gsl_matrix *W, st_out *out);

int calc_N_motifs(gsl_matrix *W,
		  double &N_edge, double &N_recip,
		  double &N_conv, double &N_div, double &N_chain,
		  double &N_other, st_out *out);
