#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<Hevolve.h>
#include<Hsplit.h>
#include<Constants.h>
#include<Deterministic.h>

using namespace std;

/* In this program I evolve the intra-cell dynamics according to the mean field equations and the inter-cell dynamic with a probabilistic setting 
How do i do this: The random part is the same! The deterministic part is obtained with a two step evolution using approximated functions (check notebook)
*/

//Note that error 6 means that interval is smaller than ts

int main(){
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t; //The time
    int i,iloop,j,n; 
    int TMAX; //Is the number of timesteps I have to do to arrive at T => T/ts
    int numprint; //Is the number of times I want to print
	int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
	ofstream filec, filet, fileN,filex;//Output files
	const char filenameN[]="ensambleN.txt"; //While ensambleN and x  will print out each of the <N> and <x>
    const char filenamex[]="ensamblex.txt";
    const char filenamet[]="time.txt";
    const char filenamec[]="parameters.txt"; //Just to print out all the parameters in the same folder at the end
    unsigned int seed; //Seed of the random number generator
    gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	double dummy,tstar; //Dummy will generally be Nc+Nd
	
	
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
	TMAX=cons.T/cons.ts;
	dummy=(int) TMAX; 
	if(dummy!=TMAX){
		TMAX=floor(cons.T/cons.ts)+1; //Here I am just adjusting for an extra step
	}
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
		Nc[0]=cons.N0*cons.x0;
		x[0]=cons.x0;
		Nd[0]=cons.N0*(1.-cons.x0);
		M=1; //I start with one cell

		//*******end of initialization*********
		
		printiterens(Nc,Nd,1,fileN,filex); //Here I print for time==0
		
		//*****Start of the time evolution***********
		
		for(j=1;j<=TMAX;j++){ //This is the time loop! It is the equivalent of the do-while in the Main file! I stop for TMAX such that t=T
			
			t=t+cons.ts; //Update the time to the new one!
			
			for(i=0;i<M;i++){ //This is the cell loop
				//*******Start of the cell evolution***********
			
				dummy=Nevolve(Nc[i],Nd[i], x[i], cons.ts, cons); //Compute the approximate value of Nc+Nd after time ts
				if(dummy>=cons.N_max){ //Check if I need to do the splitting or not
					//********Here I perform the splitting**********
					tstar=inverseN(Nc[i],Nd[i], x[i],cons); //Compute the time when N is roughly equal to N_max
					x[i]=Nevolve(Nc[i],Nd[i], x[i], tstar, cons); //Compute the value of x after tstar
					Nc[i]=cons.N_max*x[i]; //I compute the values of Nc and Nd at tstar
					Nd[i]=cons.N_max-Nc[i];
					n=createcelldeterministic(&M,i,Nc, Nd, x,cons, r); //Here I do the splitting and all the related things
					//Now I have to finish the evolution for a time step ts-tstar for the i-th and n-th cell (maybe)
					dummy=Nevolve(Nc[i],Nd[i], x[i], cons.ts-tstar, cons); //For the i-th cell
					Nc[i]=dummy*x[i]; 
					Nd[i]=dummy-Nc[i];
					x[i]=Nevolve(Nc[i],Nd[i], x[i], cons.ts-tstar, cons);
					Nc[i]=dummy*x[i]; 
					Nd[i]=dummy-Nc[i];
					if(i>n){ //So if i is bigger than n I finish the evolution, otherwise I will perform one entire step of evolution later
						dummy=Nevolve(Nc[n],Nd[n], x[n], cons.ts-tstar, cons); //The same for the n-th cell
						Nc[n]=dummy*x[n]; 
						Nd[n]=dummy-Nc[n];
						x[n]=Nevolve(Nc[n],Nd[n], x[n], cons.ts-tstar, cons);
						Nc[n]=dummy*x[n]; 
						Nd[n]=dummy-Nc[n];						
					}
				}
				else{ //It means that in this timestep there is no splitting for the i-th cell
					Nc[i]=dummy*x[i]; 
					Nd[i]=dummy-Nc[i];
					x[i]=Nevolve(Nc[i],Nd[i], x[i], cons.ts, cons);
					Nc[i]=dummy*x[i]; 
					Nd[i]=dummy-Nc[i];
				}
				
				//*********End of the single cell evolution*************
	
			}
			//***********End of the cell loop at fixed time****************
		
			//Now to check if I have to print or not
			if(t>=cons.T){ //If I am at the end of the file I print!
				printiterens(Nc,Nd,M,fileN,filex); //printing of the values in the row
				cout<<"The time is "<<t<<" and iloop is "<<iloop<<endl; //Just to check
			}
			else{ //If the time is a "multiple" of interval, then I print
				tdummy=t/cons.interval;
				dummy= (int) tdummy;
				if(tdummy==dummy){
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
    for(i=0;i<=TI;i++){ 
    	filet<<i*cons.interval<<endl;
    	}
    //If T is not a multiple of interval I need the last step:
    tdummy=cons.T/cons.interval;
	dummy= (int) tdummy;
	if(tdummy!=dummy){
		filet<<cons.T<<endl;
	}
    //************************************
    
    filet.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    filec.open(filenamec,ios::out|ios::trunc); //Now I just print all the parameters to a file!
    printparamloop(filec,cons);
    filec.close();
	
	return 0;
}


