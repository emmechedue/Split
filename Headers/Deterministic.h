#pragma once
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<Constants.h>
#include<iomanip>

using namespace std;

void computeNcNd(double x, double N, double *Nc, double *Nd){ //Here I first doublecheck the limit cases and then I do the rest!
	int temp;
	double Ncdouble,Nddouble,Ntemp,fract;
	
	fract=modf(N,&Ntemp);
	if(fract>0.5){
		Ntemp=Ntemp+1;
	}
	
	if(x==0){
		*Nd=Ntemp;
		*Nc=0.;
	}
	else{
		if(x==1){
			*Nc=Ntemp;
			*Nd=0.;
		}
		else{
			
			if(x>0.5){
				temp=ceil(Ntemp*x);
			}
			else{
				temp=floor(Ntemp*x);
			}
			Ncdouble=(double) temp;
			fract=(double) floor(Ntemp);
			Nddouble=fract-Ncdouble;
			*Nc=Ncdouble;
			*Nd=Nddouble;
			
		}
	}
	return ;
}

double Nevolve(double Nold, double x, double t, Constants cons){ //This function returns the approximate value of Nc+Nd after a time t. Nold is Nc+Nd of the old cell
	double a;
	double N;
	
	if(cons.fitness==1){ //Computing a for the 2 different fitnesses
		a=(1.+cons.p*x)*(1.+cons.s*x*(cons.b-cons.c));
	}
	else if(cons.fitness==2){
		a=(1.+cons.p*x)*(1.-cons.s*x);
	}
	else{
		a=1.+cons.p*x;
	}
	
	N=(a*cons.K*Nold*exp(a*t))/(a*cons.K+Nold*(exp(a*t)-1.));
	return N;
}

double xevolve(double x, double t, Constants cons){ //This function returns the approximate value of Nc+Nd after a time t. If x is bigger than 0.99, it automatically returns 1 (this is done in order to have fixation in x=1, that is not given by the approximation of the equations! If x<0.01 it automatically returns 0 (just to avoid to divide by zero when I define fa
	double fa;
	double newx;
	
	if(x>=0.99){
		return 1.;
	}
	else{
		if(x<=0.01){
			return 0.;
		}
		else{
			fa=1./x+cons.p-1.;
			newx=1./(fa*exp(cons.s*t)-cons.p+1.);
		}
		return newx;
	}
}

double inverseN(double N, double x, Constants cons){ //This function returns the approximate value of t for wich N=N_max. N=Nc+Nd
	double a;
	double t;
	
	if(cons.fitness==1){ //Computing a for the 2 different fitnesses
		a=(1.+cons.p*x)*(1.+cons.s*x*(cons.b-cons.c));
	}
	else if(cons.fitness==2){
		a=(1.+cons.p*x)*(1.-cons.s*x);
	}
	else{
		a=1.+cons.p*x;
	}
	
	t=(1./a)*log((cons.N_max*(a*cons.K-N))/(N*(a*cons.K-cons.N_max)));
	
	return t;
}

int createcelldeterministic(int *M, int m,double *N, double *Nc, double *Nd, double *x, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. The function returns the index of the other cell that was created
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
	if(cons.choice==1){ //I have to define N[m] and N[n] according to the 2 different models
		N[m]=cons.N0;
		N[n]=cons.N0;
	}
	else{
		N[m]=Nc[m]+Nd[m];
		N[n]=Nc[n]+Nd[n];
	}
	return n;
}
	
