#include <gsl/gsl_blas.h>

#include <iostream>
#include <vector>
#include "calc_stats_1p.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////
// calc_phat_alphahat_1p
//
// Calculate first and second order connectivity statistics from a 
// network adjacency matrix W with N_node nodes
////////////////////////////////////////////////////////////////////


int calc_phat_alphahat_1p(gsl_matrix *W, st_out *out){

  double N_edge, N_recip, N_conv, N_div, N_chain, N_other;
  
  calc_N_motifs(W, N_edge, N_recip, N_conv, N_div, N_chain, N_other, out);
// 		inDeg, outDeg, preSyn, postSyn);
  
  double c1 = out->N_nodes*(out->N_nodes-1.0);
  double c2 = out->N_nodes*(out->N_nodes-1.0)*(out->N_nodes-2.0)/2.0;
  double c3 = out->N_nodes*(out->N_nodes-1.0)*(out->N_nodes-2.0)*(out->N_nodes-3.0)/2.0;
  out->phat =N_edge/c1;
  out->alphahat_recip =N_recip/(c1/2.0*out->phat*out->phat)-1;
  out->alphahat_conv =  N_conv/(c2*out->phat*out->phat)-1;
  out->alphahat_div = N_div/(c2*out->phat*out->phat)-1;
  out->alphahat_chain = N_chain/(2*c2*out->phat*out->phat)-1;
  out->alphahat_other =  N_other/(c3*out->phat*out->phat)-1;

  return 0;

}

////////////////////////////////////////////////////////////////////
// calc_N_motifs
//
// Calculate number of edges and two edges motifs from a 
// network adjacency matrix W with N_node nodes
////////////////////////////////////////////////////////////////////

int calc_N_motifs(gsl_matrix *W,
	       double &N_edge, double &N_recip,
	       double &N_conv, double &N_div, double &N_chain,
	       double &N_other,
// 	       double *inDeg, double *outDeg,
// 	       vector<double> *preSyn, vector<double> *postSyn
	       st_out *out) {
  
  double *thedata = W->data;
  size_t tda = W->tda;
  
  double W_square_trace=0.0;
  N_edge=0;
  for(int i=0; i<out->N_nodes;i++)
    out->inDeg[i]=out->outDeg[i]=0.0;
  for(int i=0; i<out->N_nodes;i++)
  {
	  int i_tda = i*tda;
	  for(int j=0; j<out->N_nodes;j++)
	  {
		  int temp_in=thedata[i_tda+j];
		  int temp_out=thedata[j*tda+i];
		  if (temp_in == 1)	// save pre and postsynapses in corresponding arrays
		  {
			  out->preSyn.push_back(j);
			  out->inDeg[i]++;
		  }
		  if (temp_out == 1)
		  {
			  out->postSyn.push_back(j);
			  out->outDeg[i]++;
		  }
		    
		  W_square_trace += temp_in*thedata[j*tda+i];
		  N_edge += temp_in;
	  }
  }

  N_recip = W_square_trace/2;
  
  // N_chain, N_conv, and N_div
  N_chain = 0.0;
  N_conv = 0.0;
  N_div = 0.0;
  for (int i=0;i<out->N_nodes;i++){
    N_conv += out->inDeg[i]*out->inDeg[i];
    N_div += out->outDeg[i]*out->outDeg[i];
    N_chain += out->inDeg[i]*out->outDeg[i];
  }
  N_chain -= W_square_trace;
  N_conv -= N_edge;
  N_div -= N_edge;
  N_conv /= 2;
  N_div /= 2;
  
  N_other=N_edge*(N_edge-1.0)/2.0- (N_conv+N_div+N_chain+N_recip);
  
  return 0;

}

