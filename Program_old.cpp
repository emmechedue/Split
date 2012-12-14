#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<headers.h>

using namespace std;

//**********Initializations**********
const int N0=12; //Initial number of bacteria in each cell
const int Nc0=6;//Initial number of cooperators in each cell
const double T=45.; //Time when the simulation stops
const double interval=0.001; //Time step for which I print my results in fast
const int M=2000; //Number of cells
const double b=3.;
const double c=1;
const double s=0.05; //Selection's strenght
const double p=10.; //Cooperators advantage
const double K=100.; //Carrying capacity

//Note that in this program the w_s is set to 1!!!

//*****************************

int main(){
    double Nc[M], Nd[M], x[M]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t ,oldt; //t is the time and oldt will be used to check whether or not print
    int i,l,m,emme=4*M;
    double Gamma[emme]; //The array with all the partial sums
    double **G; //Matrix with all the gammas for all the cells in form of G[cell][reaction]
    double rand;
    ofstream file,file_fast;//Output file and a file where I'm not going to print everything
    const char filename[]="output.txt";
    const char fname[]="fast.txt";
    unsigned int seed; //Seed of the random number generator
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
    
    //*********Let's initialize all**********
    t=0.;
    oldt=0.;
    for(i=0; i<M; i++){ //Initialize the matrices
        Nc[i]=Nc0;
        Nd[i]=N0-Nc0;
        x[i]=Nc[i]/(Nc[i]+Nd[i]);
    }
    G=new double* [M]; //Create the Mx4 gamma matrix
    for(i=0; i<M; i++){
        G[i]=new double[4];
    }
    initializeGamma(G,Gamma,M,Nc,Nd,x,p,s,K,b,c);
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
    file<<"# M="<<M<<"  T="<<T<<"  K="<<K<<"  s="<<s<<"  p="<<p<<"  N0="<<N0<<"  Nc0="<<Nc0<<"  seed="<<seed<<endl;
    file<<"#Time  N   x"<<endl;
    myprint2(Nc,Nd,t,M,file);
    file_fast.open(fname,ios::out|ios::trunc); //Open the output's fast_file and print the results for time=0
    file_fast<<"#Results for the simulation reproducing the old results with"<<endl;
    file_fast<<"# M="<<M<<"  T="<<T<<"  K="<<K<<"  s="<<s<<"  p="<<p<<"  N0="<<N0<<"  Nc0="<<Nc0<<"  seed="<<seed<<endl;
    file_fast<<"#Time  N   x"<<endl;
    myprint2(Nc,Nd,t,M,file_fast);
    
    //*****Start of the evolution***********
     
   do{ 
        //cout<<endl<<endl<<"Gamma[emme-1] is"<<Gamma[emme-1]<<endl<<endl;
        rand=randlog(Gamma[emme-1],r);//Samples the time at wich the next reaction happens;
        //cout<<endl<<"rand= "<<rand<<endl;
        t=t+rand; //Update the time
        oldt=oldt+rand; //Update oldt
        //cout<<endl<<"oldt= "<<oldt<<endl;
        rand=gsl_rng_uniform(r)*Gamma[emme-1]; //Generates the random number to choose the reaction!
        //cout<<"check 1"<<endl;
        //cout<<"check 2"<<endl;
        l=search(Gamma,emme,rand); //Finds the reaction
        //cout<<"check 3"<<endl;
        m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
        //cout<<"check 4"<<endl;
        updateG(G,Gamma,m,Nc,Nd,x,p,s,K,b,c,emme); //Updates the G and the Gamma
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
