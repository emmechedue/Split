#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<Hevolve.h>
#include<Hsplit.h>
#include<Constants.h>
#include<Deterministic.h>
#include<boost/numeric/odeint.hpp>
#include<vector>

using namespace std;
using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;


/* In this program I solve the mean field equations using boost. The program works pretty simlarly to Deterministic.cpp but I solve with boost instead of using my method.

Unfortunately, to simplify I will have to put all the parameters here, besides using the conf file. So MAKE SURE THAT THE DATA HERE ARE THE SAME AS IN THE CONF FILE!!!!! */

const int N0=4; //Initial number of bacteria in the cell
const double x0=0.5;//Initial fraction of cooperators in the cell
const double T=18.; //Time when the simulation stops
const double interval=0.05; //Time step for which I print my results in fast
const double s=0.05; //Selection's strenght
const double c=1.; //Cost in the original fitness
const double b=3.; //Benefit in the original fitness
const double p=10.; //Cooperators advantage
const int K=100.; //Carrying capacity
const int N_max=80; //The number of bacteria in the cell s.t. the cell splits
const int M_max=1000; //The maximum number of cells
const int N_loop=300; //The number of times I iterate
const int choice=2; //It's 1 if I want the propagule model and it is 2 if I want the random splitting model. It is 3 if I am choosing a deterministic splitting
const int fitness=2; //It's 1 if I am using the original fitnesses (the one from the paper), is 2 if I am using the approximated fitnesses
const double ts=0.0001; //The timestep in the semi-deterministic model. The timestep of the "integration" over N

/* The rhs of x' = f(x) 
Here is in the form of X(1)=x and X(2)=N
*/

double faveragenumeric(double x){
	double a;
	
	if (fitness==2){
		a=1-s*x;
	}
	else{
		a=1+s*x*(b-c);
	}

	return a;
}

void nx_evolver( const state_type &X , state_type &dxdt , const double /* t */)
{
    dxdt[0] = -s*(1+p*X[0])*X[0]*(1-X[0]);
    dxdt[1] = ((1+p*X[0])*faveragenumeric(X[0])-X[1]/K)*X[1];
}


int main(){
    Constants cons;
    double N[cons.M_max],Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t; //The time
    int i,iloop,j,n; 
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
	double dummy,tstar; //Dummy will generally be Nc+Nd. 
	int tempstep; //Tempstep is for printing purposes
	state_type X(2);
	state_type dxdt(2);
	size_t steps;
	
	
	//******let's take the seed for the rng and initialize the rng******
	pfile = fopen ("/dev/urandom", "r");
	TMAX=fread (&seed, sizeof (seed), 1, pfile); //I added the TMAX= ... just to not be bothered anymore by the warnings!
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
	
	if(cons.choice==1){ //This is just to print the right model in the title!
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation that checks with boost the dynamics reproducing the propagule with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation that checks with boost the dynamics reproducing the propagule with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
	else if(cons.choice==2){
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation that checks with boost the dynamics reproducing the random splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation that checks with boost the dynamics reproducing the random splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
	else{
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation that checks with boost the dynamics reproducing deterministic splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation that checks with boost the dynamics reproducing deterministic splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
	}
		
	
	//********************************Compute TMAX************
	TMAX=cons.T/cons.ts;
	dummy=(int) TMAX; 
	if(dummy!=TMAX){
		TMAX=floor(cons.T/cons.ts)+1; //Here I am just adjusting for an extra step
	}
	tempstep=(int) floor(cons.interval/cons.ts); //I need this for printing purposes. Note that in this way I might end up not printing with the timestep I originally wanted. So I always need interval to be a multiple of ts!!!
	//*********************************
	
	//Compute how many times I have to print (hopefully smaller or equal than TMAX) and check that interval>=ts
	if(cons.interval<cons.ts){
		cout<<"Print with bigger time intervals!!!"<<endl;
		exit(6);
	}
	//*************Let's start the cycle********************
    for(iloop=0;iloop<cons.N_loop;iloop++){
    	
    
		//*********Let's initialize all**********
		t=0.;
		for(i=1;i<cons.M_max;i++){
			Nc[i]=0;
			x[i]=-1;
			Nd[i]=0;
			N[i]=0;
		}
		N[0]=cons.N0;
		x[0]=cons.x0;
		computeNcNd(x[0], N[0], &Nc[0], &Nd[0]);
		M=1; //I start with one cell

		//*******end of initialization*********
		
		printiterens(Nc,Nd,1,fileN,filex); //Here I print for time==0
		
		//*****Start of the time evolution***********
		
		for(j=1;j<=TMAX;j++){ //This is the time loop! It is the equivalent of the do-while in the Main file! I stop for TMAX such that t=T
			
			t=t+cons.ts; //Update the time to the new one!
			
			for(i=0;i<M;i++){ //This is the cell loop
				//*******Start of the cell evolution***********
		
				
				X[0]=x[i];
				X[1]=N[i];
				steps = integrate( nx_evolver , X , 0.0 , ts , ts/2 );
				x[i]=X[0];
				N[i]=X[1];
				if(N[i]>=cons.N_max){ //Check if I need to do the splitting or not
					//********Here I perform the splitting**********
					//I don't correct the x for this tiny time difference, it should work in any case!!
					computeNcNd(x[i], cons.N_max, &Nc[i], &Nd[i]); //I compute the values of Nc and Nd at tstar
					n=createcelldeterministic(&M,i,N,Nc, Nd, x,cons, r); //Here I do the splitting and all the related things
					//So there is also no evolution after the split for the i-th cell
					if(i>n){ //So if i is bigger than n I finish the evolution, otherwise I will perform one entire step of evolution later
						X[0]=x[n];
						X[1]=N[n];
						steps = integrate( nx_evolver , X , 0.0 , ts , ts/2 );
						x[n]=X[0];
						N[n]=X[1];
						computeNcNd(x[n], N[n], &Nc[n], &Nd[n]);	
					}
				}
				else{ //It means that in this timestep there is no splitting for the i-th cell	
					computeNcNd(x[i], N[i], &Nc[i], &Nd[i]);
				}
			
				//*********End of the single cell evolution*************
	
			}
			//***********End of the cell loop at fixed time****************
		
			//Now to check if I have to print or not
			if(abs(t-cons.T)<cons.ts){ //If I am at the end of the file I print!
				printiterens(Nc,Nd,M,fileN,filex); //printing of the values in the row
				cout<<"The time is "<<t<<" and iloop is "<<iloop<<endl; //Just to check
			}
			else{ //If the time is a "multiple" of interval, then I print
				if((j%tempstep)==0){
					printiterens(Nc,Nd,M,fileN,filex); //printing of the values in the row
					cout<<"The time is "<<t<<" and iloop is "<<iloop<<endl; //Just to check
				}	
			}
			// End of the part inside the time loop
		}
		
		filex<<endl; //I print the \n in the 2 files!
	  	fileN<<endl;
	
	}
	
	//********************Here ends the main loop*******************
	
	//********Start printing the time
	t=(int)floor(cons.T/cons.interval); //I need it to print the time
	
	filet.open(filenamet,ios::out|ios::trunc); //Printing the time
    for(i=0;i<=t;i++){ 
    	filet<<i*cons.interval<<endl;
    	}
    //If T is not a multiple of interval I need the last step:
    tstar=cons.T/cons.interval;
	dummy= (int) tstar;
	if(tstar!=dummy){
		filet<<cons.T<<endl;
	}
    //************************************
    
    filet.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    filec.open(filenamec,ios::out|ios::trunc); //Now I just print all the parameters to a file!
    printparamloopdeterministic(filec,cons);
    filec.close();
	
	return 0;
}
