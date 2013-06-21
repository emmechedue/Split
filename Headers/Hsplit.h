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

/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
void fill1(int *C, int *D, double x, int N0, gsl_rng *r){ //This one gives C (the numbers of the cooperators in the cell) according to Bin(x,N0) and D as N0-C. So this is the propagule model
	int n;
	
	n=gsl_ran_binomial(r,x,N0);
	*C=n;
	*D=N0-n;
	return;
}

void fill2(int *C, int *D, double Nc, double Nd, gsl_rng *r){//This one gives C (the numbers of the cooperators in the cell) according to Bin(0.5,Nc) and D according to Bin(0.5,Nd), where Nc and Nd are the # of cooperators and defectors in the cell that just splitted. So this is the random splitting model
	int nd,nc;
	
	if(Nc<=0){
		nc=0;
	}
	else{
		nc=gsl_ran_binomial(r,0.5,Nc);
	}
	if (Nd<=0){
		nd=0;
	}
	else{
		nd=gsl_ran_binomial(r,0.5,Nd);
	}
	*C=nc;
	*D=nd;
	
	return;
}


void fillcells(int n, int m, double *Nc, double *Nd, double *x, Constants cons, gsl_rng *r){ //n is the cell that I want to fill, m is the cell from which I take the parameters, choice can take values 1 (for the first method) and 2 for the second method
	int C,D;
	
	//Here I choose which splitting I want to do
	if(cons.choice==1){ 
		fill1(&C,&D,x[m],cons.N0,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
		Nc[n]=C;
		Nd[n]=D;
		x[n]=Nc[n]/(Nc[n]+Nd[n]);
		fill1(&C,&D,x[m],cons.N0,r);
		Nc[m]=C;
		Nd[m]=D;
		x[m]=Nc[m]/(Nc[m]+Nd[m]);
		
	}
	else{ 
		do{ //Just to be sure that I am not creating an empty cell!!
			fill2(&C,&D,Nc[m],Nd[m],r);  // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
		}while((C==0)&&(D==0));
		Nc[n]=C;
		Nd[n]=D;
		x[n]=Nc[n]/(Nc[n]+Nd[n]);
		Nc[m]=Nc[m]-C; //The new Nc[m] is the old Nc[m] - the cooperators in the new cell => I'm conserving the total number of bacteria
		Nd[m]=Nd[m]-D;
		x[m]=Nc[m]/(Nc[m]+Nd[m]);
	}
	
	return;
}
	

int createcell(int M, int m,double *Nc, double *Nd, double *x, double *Gamma, double **G, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. the function returns the new value of M
	int n; //This is the index of one of the two new cells
	double rand;
	
	//********************Determine n**********
	if(M>=cons.M_max){
		rand=gsl_rng_uniform(r)*(cons.M_max-1); //Generate a uniform random number between 0 and M_max-1
		n=ceil(rand)-1; //I subtract -1 because the full array goes from 0 to M_max-1
		if(n==m){
			n=n+1; //If I extract just the cell I am already splitting, I take the next one
		} 
	}
	else{
		n=M; //not M+1 because of the index problem (it would be M+1-1)
		M=M+1;
	}
	
	//********creates the new cells and updates the Gamma and the G*****************
	fillcells(n,m,Nc,Nd,x,cons,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
	//cout<<"First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[n]+Nd[n]<<" bacteria"<<endl;
	if (n<m){ //  in the function I always have to have n<m
		updatebothG(G,Gamma,n,m,Nc,Nd,x,cons,4*M);
	}
	else{
		updatebothG(G,Gamma,m,n,Nc,Nd,x,cons,4*M);
	}
	return M;
}
	
		
