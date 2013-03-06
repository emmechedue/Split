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

//In this program I print the averages over the ensambles and the values of <N> and <x> for every single iteration
//Please note that in this program I will only use oldt and interval, so in this code: interval==intervalens
//Convention: average inside the iteration: < > ; ensamble average E_{ }
//Note that in this program the w_s is set to 1!!!

//*****************************

int main(){
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t ,oldt; //t is the time and oldt will be used to check whether or not print
    int i,l,m,emme=4*cons.M_max,iloop;
    double Gamma[emme]; //The array with all the partial sums
    double **G; //Matrix with all the gammas for all the cells in form of G[cell][reaction]
    double rand; 
    int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
    ofstream /*file, */filet, fileN,filex;//Output files
    //const char filename[]="output.txt"; //Here output.txt will output for each time step the ensamble average of <N> and <x>  and also the ensamble average of M
    const char filenameN[]="ensambleN.txt"; //While ensambleN and x  will print out each of the <N> and <x>
    const char filenamex[]="ensamblex.txt";
    const char filenamet[]="time.txt";
    unsigned int seed; //Seed of the random number generator
    /*double Nav[cons.N_loop]; //The arrays where I will store the values of <N> and <x>
    double xav[cons.N_loop];*/
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	double TI; //i need it to print the time!
	int dummy, enne; //Dummy is a dummy index needed for small loops, enne is taking care (in case) of  how many times is rand bigger than interval
      
      
     

   //******let's take the seed for the rng and initialize the rng******
	pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
	
	if(cons.choice==1){ //This is just to print the right model in the title!
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation reproducing the propagule with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation reproducing the propagule with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
	else{
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation reproducing the random splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		fileN<<"#In the form of E_{N[m][t]}"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation reproducing the random splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<"  N_loop="<<cons.N_loop<<endl;
		filex<<"#In the form of E_{x[m][t]}"<<endl;
		}
    
    //*****************************
    
    TI=(int)floor(cons.T/cons.interval);
    
    //*************Let's start the cycle********************
    for(iloop=0;iloop<cons.N_loop;iloop++){
    	
    
		//*********Let's initialize all**********
		t=0.;
		oldt=0.;
		Nc[0]=cons.N0*cons.x0;
		x[0]=cons.x0;
		Nd[0]=cons.N0*(1.-cons.x0);
		M=1; //I start with one cell

		
		G=new double* [cons.M_max]; //Create the Mx4 gamma matrix
		for(i=0; i<cons.M_max; i++){
		    G[i]=new double[4];
		}
		initializeGamma(G,Gamma,Nc,Nd,x,cons);
		//*******end of initialization*********
		
		printiterens(Nc,Nd,1,fileN,filex); //Here I print for time==0
		
		//*****Start of the evolution***********
		 
	   do{ 
		   rand=randlog(Gamma[4*M-1],r);//Samples the time at wich the next reaction happens;
		   t=t+rand; //Update the time
		   if(rand>cons.interval){ //Here is to check if I have to reprint the old situation before update the system!
		   		enne=floor(rand/cons.interval);
		   		for(dummy=0;dummy<enne;dummy ++){
		   			printiterens(Nc,Nd,M,fileN,filex);
		   			}
		   		rand=rand-cons.interval*enne;
		   }
		   oldt=oldt+rand; //Update oldt
		   //oldtensamble=oldtensamble+rand; //Update oldtensamble
		   rand=gsl_rng_uniform(r)*Gamma[4*M-1]; //Generates the random number to choose the reaction!
		   l=search(Gamma,4*M,rand); //Finds the reaction
		   m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
		   if(check(Nc, Nd, cons, m)==true){ //Of course I need to check if I have to split the cell or not
		   		M=createcell(M, m, Nc, Nd, x, Gamma, G, cons, r); 
				//cout<<endl<<endl<<"Now in the main: First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[1]+Nd[1]<<" bacteria"<<endl<<endl;
				//Here I do everything, I create the cell, I update the cells and then update the Gamma and G
				}
		   else{ //Of course if no cell splits, I just update the G and the Gamma, print and then sample for another reaction
				updateG(G,Gamma,m,Nc,Nd,x,cons,4*M); //Updates the G and the Gamma
				}
				
				
		  if(oldt>=cons.interval){ //Checks whether I have to print or not
					printiterens(Nc,Nd,M,fileN,filex); //printing of the values in the row
					oldt=oldt -cons.interval; //Subract by oldt the value of interval to start counting again 
					//count++;
					cout<<"The time is "<<t<<" and iloop is "<<iloop<<endl; //Just to check
				}
		
			/*   if(oldt>=cons.interval){ //Checks whether I have to print or not on ensamble.txt
					myprintensamble2(Nc,Nd,t,M,fileN,filex); //Printing the results on file ensamble; to create the movie
					oldt=oldt -cons.intervalens; //Subract by oldtensamble the value of intervalens to start counting again
					//cout<<"The time is "<<t<<endl; //Just to check
				}*/
				
	  }while(t<=cons.T);
	  //cout<<endl<<endl<<"gamma= "<<Gamma[4*M-1]<<endl<<endl;
		
	  filex<<endl; //I print the \n in the 2 files!
	  fileN<<endl;

    }
    
    
    //********************Here ends the loop
   
    filet.open(filenamet,ios::out|ios::trunc); 
    for(iloop=0;iloop<=TI;iloop++){
    	filet<<iloop*cons.interval<<endl;
    	}
    
    
    filet.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    
    
    
    return 0;
}
