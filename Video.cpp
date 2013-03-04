//The only difference between this program and Single_loop is that here I'm printing even if the group is empty
//Then in the python analysis I should do something like: if (M!=0) {print!}

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

//In this program I perform only a single loop but I save everything (meaning the behavior of any single cell)
//Note that in this program the w_s is set to 1!!!

//*****************************

int main(){
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t ,oldt; //t is the time and oldt will be used to check whether or not print
    int i,l,m,emme=4*cons.M_max;
    double Gamma[emme]; //The array with all the partial sums
    double **G; //Matrix with all the gammas for all the cells in form of G[cell][reaction]
    double rand;
    int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
    ofstream file, fileN,filex;//Output files
    const char filename[]="output.txt";
    const char filenameN[]="ensambleN.txt";
    const char filenamex[]="ensamblex.txt";
    unsigned int seed; //Seed of the random number generator
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	int dummy, enne; //Dummy is a dummy index needed for small loops, enne is taking care (in case) of  how many times is rand bigger than interval
    
    //*********Let's initialize all**********
    t=0.;
    oldt=0.;
    Nc[0]=cons.N0*cons.x0;
    x[0]=cons.x0;
    Nd[0]=cons.N0*(1.-cons.x0);
    M=1; //I start with one cell
    
    //Here I fill all the rest of the array with 0s
    for(i=1;i<cons.M_max;i++){
    	Nc[i]=0;
		x[i]=0;
		Nd[i]=0;
	}
    
    G=new double* [cons.M_max]; //Create the Mx4 gamma matrix
    for(i=0; i<cons.M_max; i++){
        G[i]=new double[4];
    }
    initializeGamma(G,Gamma,Nc,Nd,x,cons);
    //*******end of initialization*********
    
    //******let's take the seed for the rng and initialize the rng******
    pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
    
    file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
    file<<"#Results for the simulation reproducing the splitting with"<<endl;
    file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
    file<<"#Time  N   x    M"<<endl;
    myprint2(Nc,Nd,t,M,file);
    fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
    fileN<<"#Results for the simulation reproducing the splitting with"<<endl;
    fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
    fileN<<"#In the form of N[t][m]"<<endl;
    filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
    filex<<"##Results for the simulation reproducing the splitting with"<<endl;
    filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
    filex<<"#In the form of x[t][m]"<<endl;
    myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
    
    //*****Start of the evolution***********
     
   do{ 
        rand=randlog(Gamma[4*M-1],r);//Samples the time at wich the next reaction happens;
        t=t+rand; //Update the time
        if(rand>cons.interval){ //Here is to check if I have to reprint the old situation before update the system!
		  		enne=floor(rand/cons.interval);
		  		for(dummy=0;dummy<enne;dummy ++){
		  			myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
		  			myprint2(Nc,Nd,t,M,file);
		  			}
		  		rand=rand-cons.interval*enne;
		}
		oldt=oldt+rand; //Update oldt
        
		rand=gsl_rng_uniform(r)*Gamma[4*M-1]; //Generates the random number to choose the reaction!
		l=search(Gamma,4*M,rand); //Finds the reaction
        
		m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
        
		if(check(Nc, Nd, cons, m)==true){ //Of course I need to check if I have to split the cell or not
			M=createcell(M, m, Nc, Nd, x, Gamma, G, cons, r); //Here I do everything, I create the cell, I update the cells and then update the Gamma and G
			//cout<<endl<<endl<<"Now in the main: First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[1]+Nd[1]<<" bacteria"<<endl<<endl;
		}
		else{ //Of course if no cell splits, I just update the G and the Gamma, print and then sample for another reaction
		updateG(G,Gamma,m,Nc,Nd,x,cons,4*M); //Updates the G and the Gamma
		}
        
        
		if(oldt>=cons.interval){ //Checks whether I have to print or not
			myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
			myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
        	oldt=oldt -cons.interval; //Subract by oldt the value of interval to start counting again
        	cout<<"The time is "<<t<<endl; //Just to check
		}
		
        
	}while(t<=cons.T);
    
    file.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    
    
    return 0;
}
