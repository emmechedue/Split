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

using namespace std;

//In this program I perform only a single loop but I save everything (meaning the behavior of any single cell)
//Note that in this program the w_s is set to 1!!!

//*****************************

int main(){
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t; //The time
    int i,j,n; 
    int TMAX; //Is the number of timesteps I have to do to arrive at T => T/ts
	int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
	ofstream filec, file, fileN,filex;//Output files
	const char filename[]="output.txt";
    const char filenameN[]="ensambleN.txt";
    const char filenamex[]="ensamblex.txt";
    const char filenamec[]="parameters.txt"; //Just to print out all the parameters in the same folder at the end
    unsigned int seed; //Seed of the random number generator
    gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	double dummy,tstar; //Dummy will generally be Nc+Nd. 
	int tempstep; //Tempstep is for printing purposes
    
    //*********Let's initialize all**********
	t=0.;
	for(i=1;i<cons.M_max;i++){
		Nc[i]=0;
		x[i]=-1;
		Nd[i]=0;
	}
	Nc[0]=cons.N0*cons.x0;
	x[0]=cons.x0;
	Nd[0]=cons.N0*(1.-cons.x0);
	M=1; //I start with one cell

	//*******end of initialization*********
    
    //******let's take the seed for the rng and initialize the rng******
	pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
    
    if(cons.choice==1){
		file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
		file<<"#Results for the simulation with deterministic intra-cell dynamics reproducing the propagule with"<<endl;
		file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		file<<"#Time  N   x    M"<<endl;
		myprint2(Nc,Nd,t,M,file);
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation with deterministic intra-cell dynamics reproducing the propagule with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		fileN<<"#In the form of N[t][m]"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation with deterministic intra-cell dynamics reproducing the propagule with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		filex<<"#In the form of x[t][m]"<<endl;
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
	}
	else{
		file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
		file<<"#Results for the simulation with deterministic intra-cell dynamics reproducing the random splitting with"<<endl;
		file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		file<<"#Time  N   x    M"<<endl;
		myprint2(Nc,Nd,t,M,file);
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation with deterministic intra-cell dynamics reproducing the random splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		fileN<<"#In the form of N[t][m]"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation with deterministic intra-cell dynamics reproducing the random splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		filex<<"#In the form of x[t][m]"<<endl;
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
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
	
    
    //*****Start of the evolution***********
     
   //*****Start of the time evolution***********
		
	for(j=1;/*j<=TMAX*/j<=30;j++){ //This is the time loop! It is the equivalent of the do-while in the Main file! I stop for TMAX such that t=T
		
		t=t+cons.ts; //Update the time to the new one!
		
		for(i=0;i<M;i++){ //This is the cell loop
			//*******Start of the cell evolution***********
		
			dummy=Nevolve(Nc[i]+Nd[i], x[i], cons.ts, cons); //Compute the approximate value of Nc+Nd after time ts
			if(dummy>=cons.N_max){ //Check if I need to do the splitting or not
				//********Here I perform the splitting**********
				tstar=inverseN(Nc[i]+Nd[i], x[i],cons); //Compute the time when N is roughly equal to N_max
				x[i]=xevolve(x[i], tstar, cons); //Compute the value of x after tstar
				Nc[i]=cons.N_max*x[i]; //I compute the values of Nc and Nd at tstar
				Nd[i]=cons.N_max-Nc[i];
				n=createcelldeterministic(&M,i,Nc, Nd, x,cons, r); //Here I do the splitting and all the related things
				//Now I have to finish the evolution for a time step ts-tstar for the i-th and n-th cell (maybe)
				dummy=Nevolve(Nc[i]+Nd[i], x[i], cons.ts-tstar, cons); //For the i-th cell
				x[i]=xevolve(x[i], cons.ts-tstar, cons);
				Nc[i]=dummy*x[i]; 
				Nd[i]=dummy-Nc[i];
				if(i>n){ //So if i is bigger than n I finish the evolution, otherwise I will perform one entire step of evolution later
					dummy=Nevolve(Nc[n]+Nd[n], x[n], cons.ts-tstar, cons); //The same for the n-th cell
					x[n]=xevolve(x[n], cons.ts-tstar, cons);
					Nc[n]=dummy*x[n]; 
					Nd[n]=dummy-Nc[n];						
				}
			}
			else{ //It means that in this timestep there is no splitting for the i-th cell
				x[i]=xevolve(x[i], cons.ts, cons);
				Nc[i]=dummy*x[i];
				Nd[i]=dummy-Nc[i];
			}
			
			//*********End of the single cell evolution*************
	
		}
		//***********End of the cell loop at fixed time****************
		
		//Now to check if I have to print or not
		if(abs(t-cons.T)<cons.ts){ //If I am at the end of the file I print!
			myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
			myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
			cout<<"The time is "<<t<<endl; //Just to check
		}
		else{ //If the time is a "multiple" of interval, then I print
			/*tstar=t/cons.interval;
			dummy=floor(tstar); */
			if((j%tempstep)==0){
				myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
				myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
				cout<<"The time is "<<t<<endl; //Just to check
			}	
		}
		myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
		// End of the part inside the time loop
       		
	}
    
    file.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    filec.open(filenamec,ios::out|ios::trunc); //Now I just print all the parameters to a file!
    printparamnoloop(filec,cons);
    filec.close();
        
    return 0;
}
