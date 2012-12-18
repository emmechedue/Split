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

//Note that in this program the w_s is set to 1!!!

//*****************************

int main(){
    Constants consta;
    double Nc[consta.M_max], Nd[consta.M_max], x[consta.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t ,oldt; //t is the time and oldt will be used to check whether or not print
    int i,l,m,emme=4*consta.M_max;
    double Gamma[emme]; //The array with all the partial sums
    double **G; //Matrix with all the gammas for all the cells in form of G[cell][reaction]
    double rand;
    int M; //it can go from 1 to M_max and it's just to not waste time taking into account empty cells
    ofstream file,;//Output file 
    const char filename[]="output.txt";
    unsigned int seed; //Seed of the random number generator
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
    
    //*********Let's initialize all**********
    t=0.;
    oldt=0.;
    Nc[0]=consta.N0*consta.x0;
    x[0]=consta.x0;
    Nd[0]=consta.N0*(1.-consta.x0);
    M=1; //I start with one cell

    
    G=new double* [consta.M_max]; //Create the Mx4 gamma matrix
    for(i=0; i<consta.M_max; i++){
        G[i]=new double[4];
    }
    initializeGamma(G,Gamma,Nc,Nd,x,consta);
    //*******end of initialization*********
    
    //******let's take the seed for the rng and initialize the rng******
    pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
    
    file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
    file<<"#Results for the simulation reproducing the old results with"<<endl;
    file<<"# M_max="<<consta.M_max<<"  T="<<consta.T<<"  K="<<consta.K<<"  s="<<consta.s<<"  p="<<consta.p<<"  N0="<<consta.N0<<"  x0="<<consta.x0<<"  N_max="<<consta.N_max<<"  seed="<<seed<<endl;
    file<<"#Time  N   x"<<endl;
    myprint2(Nc,Nd,t,M,file);
    
    //*****Start of the evolution***********
     
   do{ 
        //cout<<endl<<endl<<"Gamma[emme-1] is"<<Gamma[emme-1]<<endl<<endl;
        rand=randlog(Gamma[4*M-1],r);//Samples the time at wich the next reaction happens;
        //cout<<endl<<"rand= "<<rand<<endl;
        t=t+rand; //Update the time
        oldt=oldt+rand; //Update oldt
        //cout<<endl<<"oldt= "<<oldt<<endl;
        rand=gsl_rng_uniform(r)*Gamma[4*M-1]; //Generates the random number to choose the reaction!
        //cout<<"check 1"<<endl;
        //cout<<"check 2"<<endl;
        l=search(Gamma,4*M,rand); //Finds the reaction
        //cout<<"check 3"<<endl;
        m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
        //cout<<"check 4"<<endl;
        if(check(Nc, Nd, consta, m)==true){//Here I have to do the splitting } 
        else{ //Of course if no cell splits, I just update the G and the Gamma, print and then sample for another reaction
        updateG(G,Gamma,m,consta,4*M); //Updates the G and the Gamma
        }
        //cout<<"check 5"<<endl;
        //myprint2(Nc,Nd,t,M,file); //Prints N average and x average at time t
        //cout<<"check 6"<<endl;
        if(oldt>=interval){ //Checks whether I have to print or not
        	myprint2(Nc,Nd,t,M,file_fast); //Printing the results on file fast. To create a picture
        	oldt=oldt -interval; //Subract by oldt the value of interval to start counting again
        	cout<<"The time is "<<t<<endl; //Just to check
        }
    }while(t<=T);
    
    file.close(); //Closing the files of output!
    file_fast.close();
    
    return 0;
}
