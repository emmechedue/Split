#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<Hevolve.h>
#include<Hsplit.h>
#include<Constants.h>

using namespace std;

/* In this program I evolve the intra-cell dynamics according to the mean field equations and the inter-cell dynamic with a probabilistic setting */

int main(){
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t; //The time
    double ts=0.001; //The timestep of the "integration" over N
    int i,iloop,j; 
    int TMAX; //Is the number of timesteps I have to do to arrive at T => T/ts
	int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
	ofstream filec, filet, fileN,filex;//Output files
	const char filenameN[]="ensambleN.txt"; //While ensambleN and x  will print out each of the <N> and <x>
    const char filenamex[]="ensamblex.txt";
    const char filenamet[]="time.txt";
    const char filenamec[]="parameters.txt"; //Just to print out all the parameters in the same folder at the end
    unsigned int seed; //Seed of the random number generator
    gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	double dummy;
	
	
	//******let's take the seed for the rng and initialize the rng******
	pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
	
	if(cons.choice==1){ //This is just to print the right model in the title!
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation with deterministic intra-cell dynamics reproducing the propagule with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation with deterministic intra-cell dynamics reproducing the propagule with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
	else{
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"Results for the simulation with deterministic intra-cell dynamics reproducing the random splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation with deterministic intra-cell dynamics reproducing the random splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
	
	//********************************Compute TMAX************
	TMAX=cons.T/ts;
	t=(int) TMAX; //Note that here I am just using t as dummy variable
	if(t!=TMAX){
		TMAX=floor(cons.T/ts)+1; //Here I am just adjusting for an extra step
	}
	//*********************************
		
	//*************Let's start the cycle********************
    for(iloop=0;iloop<cons.N_loop;iloop++){
    	
    
		//*********Let's initialize all**********
		t=0.;
		Nc[0]=cons.N0*cons.x0;
		x[0]=cons.x0;
		Nd[0]=cons.N0*(1.-cons.x0);
		M=1; //I start with one cell

		//*******end of initialization*********
		
		printiterens(Nc,Nd,1,fileN,filex); //Here I print for time==0
		
		//*****Start of the time evolution***********
		
		for(j=1;j<=TMAX;j++){ //This is the time loop! It is the equivalent of the do-while in the Main file! I stop before 
			
			t=t+ts; //Update the time to the new one!
			
			for(i=0;i<M;i++){
			
			//*******Start of the cell evolution***********
			
			
		}
		
	}
	
	return 0;
}


