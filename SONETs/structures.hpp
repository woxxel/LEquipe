#include <vector>

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

using namespace std;

struct st_out
{
	string save_str;
	int N_nodes;
	
	double phat;
	double alphahat_recip;
	double alphahat_conv;
	double alphahat_div;
	double alphahat_chain;
	double alphahat_other;
	
	//assign size dynamically
	vector<double> inDeg;
	vector<double> outDeg;
	
	vector<double> preSyn;
	vector<double> postSyn;
};

#endif