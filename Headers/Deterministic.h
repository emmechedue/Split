#pragma once
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Constants.h>
#include<iomanip>

using namespace std;

double Nevolve(double Nc, double Nd, double x, double t, Constants cons){ //This function returns the approximate value of Nc+Nd after a time t
	double a;
	double N,Nold;
	
	if(cons.fitness==1){ //Computing a for the 2 different fitnesses
		a=(1.+cons.p*x)*(1.+cons.s*x*(cons.b-cons.c));
	}
	else{
		a=(1.+cons.p*x)*(1.-cons.s*x);
	}
	Nold=Nc+Nd;
	
	N=(a*cons.K*Nold*exp(a*t))/(a*cons.K+Nold*(exp(a*t)-1.));
	
	return N;
}

double xevolve(double Nc, double Nd, double x, double t, Constants cons){ //This function returns the approximate value of Nc+Nd after a time t
	double c;
	double newx;
	
	c=1./x+cons.p-1.;
	
	newx=1./(c*exp(cons.s*t)-cons.p+1.);
	
	return newx;
}

double inverseN(double Nc, double Nd, double x, Constants cons){ //This function returns the approximate value of t for wich N=N_max
	double a;
	double t;
	double N;
	
	if(cons.fitness==1){ //Computing a for the 2 different fitnesses
		a=(1.+cons.p*x)*(1.+cons.s*x*(cons.b-cons.c));
	}
	else{
		a=(1.+cons.p*x)*(1.-cons.s*x);
	}
	N=Nc+Nd;
	
	t=(1./a)*log((cons.N_max*(a*cons.K-N))/(N*(a*cons.K-cons.N_max)));
	
	return t;
}

int createcelldeterministic(int *M, int m,double *Nc, double *Nd, double *x, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. the function returns the index of the other cell that was created
	int n; //This is the index of one of the two new cells
	double rand;
	
	//********************Determine n**********
	if(*M>=cons.M_max){
		rand=gsl_rng_uniform(r)*(cons.M_max-1); //Generate a uniform random number between 0 and M_max-1
		n=ceil(rand)-1; //I subtract -1 because the full array goes from 0 to M_max-1
		if(n==m){
			n=n+1; //If I extract just the cell I am already splitting, I take the next one
		} 
	}
	else{
		n=*M; //not M+1 because of the index problem (it would be M+1-1)
		*M=*M+1;
	}
	
	//********creates the new cells ************
	fillcells(n,m,Nc,Nd,x,cons,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
	//cout<<"First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[n]+Nd[n]<<" bacteria"<<endl;
	
	return n;
}
	
