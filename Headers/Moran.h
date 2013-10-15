#pragma once
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Hevolve.h>
#include<gsl/gsl_rng.h>
#include<Constants.h>
#include<gsl_randist.h>

using namespace std;

//In this header I can decide how many groups have to be cooperative!

const double num_coop=3; //This is the number of groups that start as cooperative! All the M_max - num_coop will start as defectors!

void initializeNandx(int M_max, double *x, double *N){ //Here I initialize all the cells to start with 40 agents. Ther first num_coop groups will be cooperative, the others defectors!
	int i;
	
	for(i=0;i < num_coop; i++){ //Initialized the cooperative groups!
		x[i]=1.;
		N[i]=40.;
	}
	for(i=num_coop; i<M_max;i++){ //Initializing all the rest of groups to defective
		x[i]=0.;
		N[i]=40.;
	}
	
	return;
}


