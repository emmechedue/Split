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
    Constants cons;
    double Nc[cons.M_max], Nd[cons.M_max], x[cons.M_max]; //Coop. #, Def. # and fraction of coop. ****In form of N[cell]
    double t ,oldt, oldtensamble; //t is the time and oldt will be used to check whether or not print
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
    
    //*********Let's initialize all**********
    t=0.;
    oldt=0.;
    oldtensamble=0.;
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
    
    //******let's take the seed for the rng and initialize the rng******
    pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
    
    file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
    file<<"#Results for the simulation reproducing the old results with"<<endl;
    file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
    file<<"#Time  N   x"<<endl;
    myprint2(Nc,Nd,t,M,file);
    fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
    fileN<<"#Results for the simulation reproducing the old results with"<<endl;
    fileN<<"# M="<<M<<"  T="<<T<<"  K="<<K<<"  s="<<s<<"  p="<<p<<"  N0="<<N0<<"  Nc0="<<Nc0<<"  seed="<<seed<<endl;
    fileN<<"#In the form of N[t][m]"<<endl;
    filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
    filex<<"#Results for the simulation reproducing the old results with"<<endl;
    filex<<"# M="<<M<<"  T="<<T<<"  K="<<K<<"  s="<<s<<"  p="<<p<<"  N0="<<N0<<"  Nc0="<<Nc0<<"  seed="<<seed<<endl;
    filex<<"#In the form of x[t][m]"<<endl;
    myprintensamble2(Nc,Nd,t,M,fileN,filex);
    
    //*****Start of the evolution***********
     
   do{ 
        //cout<<endl<<endl<<"Gamma[emme-1] is"<<Gamma[emme-1]<<endl<<endl;
        rand=randlog(Gamma[4*M-1],r);//Samples the time at wich the next reaction happens;
        //cout<<endl<<"rand= "<<rand<<endl;
        t=t+rand; //Update the time
        oldt=oldt+rand; //Update oldt
        oldtensamble=oldtensamble+rand; //Update oldtensamble
        //cout<<endl<<"oldt= "<<oldt<<endl;
        rand=gsl_rng_uniform(r)*Gamma[4*M-1]; //Generates the random number to choose the reaction!
        //cout<<"check 1"<<endl;
        //cout<<"check 2"<<endl;
        l=search(Gamma,4*M,rand); //Finds the reaction
        //cout<<"check 3"<<endl;
        m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
        //cout<<"check 4"<<endl;
        if(check(Nc, Nd, cons, m)==true){ //Of course I need to check if I have to split the cell or not
        	M=createcell(M, m, Nc, Nd, x, Gamma, G, cons, r); //Here I do everything, I create the cell, I update the cells and then update the Gamma and G
        }
        else{ //Of course if no cell splits, I just update the G and the Gamma, print and then sample for another reaction
        updateG(G,Gamma,m,Nc,Nd,x,cons,4*M); //Updates the G and the Gamma
        }
        //cout<<"check 5"<<endl;
        //myprint2(Nc,Nd,t,M,file); //Prints N average and x average at time t
        //cout<<"check 6"<<endl;
        if(oldt>=interval){ //Checks whether I have to print or not
        	myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
        	oldt=oldt -interval; //Subract by oldt the value of interval to start counting again
        	cout<<"The time is "<<t<<endl; //Just to check
        }
        if(oldtensamble>=intervalens){ //Checks whether I have to print or not on ensamble.txt
        	myprintensamble2(Nc,Nd,t,M,fileN,filex); //Printing the results on file ensamble; to create the movie
        	oldtensamble=oldtensamble -intervalens; //Subract by oldtensamble the value of intervalens to start counting again
        	//cout<<"The time is "<<t<<endl; //Just to check
        }
    }while(t<=T);
    
    file.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    
    return 0;
}
