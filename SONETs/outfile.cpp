#include <gsl/gsl_blas.h>

#include <iostream>
#include <vector>

#include "outfile.hpp"

void save_statistics(st_out *out)
{

	//!   save data in NcFile
	NcFile TopoSv(out->save_str.c_str(),NcFile::Replace);
	
	NcVar *phatP = TopoSv.add_var("phat", ncDouble);
	NcVar *alphahat_recipP = TopoSv.add_var("alphahat_recip",ncDouble);
	NcVar *alphahat_convP = TopoSv.add_var("alphahat_conv",ncDouble);
	NcVar *alphahat_divP = TopoSv.add_var("alphahat_div",ncDouble);
	NcVar *alphahat_chainP = TopoSv.add_var("alphahat_chain",ncDouble);
	
	phatP->put(&out->phat);
	alphahat_recipP->put(&out->alphahat_recip);
	alphahat_convP->put(&out->alphahat_conv);
	alphahat_divP->put(&out->alphahat_div);
	alphahat_chainP->put(&out->alphahat_chain);
	
	
	NcDim* N_nodesP = TopoSv.add_dim("N_nodes",out->N_nodes);
	NcVar *inDegP = TopoSv.add_var("inDegree",ncDouble,N_nodesP);
	NcVar *outDegP = TopoSv.add_var("outDegree",ncDouble,N_nodesP);
	
	inDegP->put(&out->inDeg[0],out->N_nodes);
	outDegP->put(&out->outDeg[0],out->N_nodes);
	
	int N_Syn = out->preSyn.size();
	NcDim* SynapsesP = TopoSv.add_dim("Synapses",N_Syn);
	NcVar *preSynP = TopoSv.add_var("preSyn",ncInt,SynapsesP);
	NcVar *postSynP = TopoSv.add_var("postSyn",ncInt,SynapsesP);
	
	preSynP->put(&out->preSyn.front(),N_Syn);
	postSynP->put(&out->postSyn.front(),N_Syn);
	
	
}