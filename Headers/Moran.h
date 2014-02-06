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

const double num_coop=6; //This is the number of groups that start as cooperative! All the M_max - num_coop will start as defectors!

void initializeNandx(double *x, double *Nc, double *Nd, Constants cons){ //Here I initialize all the cells to start with 40 agents. Ther first num_coop groups will be cooperative, the others defectors!
	int i;
	
	for(i=0;i < num_coop; i++){ //Initialized the cooperative groups!
		x[i]=1.;
		Nc[i]=40.;
		Nd[i]=0.;
	}
	for(i=num_coop; i<cons.M_max;i++){ //Initializing all the rest of groups to defective
		x[i]=0.;
		Nd[i]=40.;
		Nc[i]=0.;
	}
	
	return;
}

void initializeallgammas(double **G, double *Gamma,double *Nc, double *Nd, double *x, Constants cons){
	int i,j; // i is going to be the group index and j the reaction index 
	double average;
	
	initializeGamma(G,Gamma,Nc, Nd, x, cons); //I initialize the first gamma filling it up for i=0
	j=4; // j starts from 4 becauase Gamma[4] (i.e. the 5th reaction) is the first one I want to save!
	for(i=1;i<cons.M_max;i++){
		average=faverage(x[i],cons);
		G[i][0]=Nc[i]*g(x[i],cons)*fcoop(x[i],cons)/average;
		Gamma[j]=Gamma[j-1]+G[i][0];
		j++;
		G[i][1]=Nc[i]*d(Nc[i],Nd[i],cons); 
		Gamma[j]=Gamma[j-1]+G[i][1];
		j++;
		G[i][2]=g(x[i],cons)*Nd[i]*fdef(x[i],cons)/average;
		Gamma[j]=Gamma[j-1]+G[i][2];
		j++;
		G[i][3]=Nd[i]*d(Nc[i],Nd[i],cons); 
		Gamma[j]=Gamma[j-1]+G[i][3];
		j++;
	}	
	
	return ;
}

	
void initializeNandxdeterministic(double *x, double *Nc, double *Nd,double *N, Constants cons){ //Here I initialize all the cells to start with 40 agents. Ther first num_coop groups will be cooperative, the others defectors!
	int i;
	
	for(i=0;i < num_coop; i++){ //Initialized the cooperative groups!
		x[i]=1.;
		Nc[i]=40.;
		Nd[i]=0.;
		N[i]=40.;
	}
	for(i=num_coop; i<cons.M_max;i++){ //Initializing all the rest of groups to defective
		x[i]=0.;
		Nd[i]=40.;
		Nc[i]=0.;
		N[i]=40.;
	}
	
	return;
}
