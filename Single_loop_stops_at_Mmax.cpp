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

/********************		IMPORTANT:                  ****************************************************

      THIS PROGRAM DOES THE SAME THING THAT SINGLE LOOP DOES BUT THE CONDITION TO STOP THE LOOP IS NOT THE TIME BUT I STOP AS SOON AS I REACH THE MAXIMUM NUMBER OF GROUPS!!!!!

****************************************************************************************************************************************************************************************/


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
    ofstream filec, file, fileN,filex;//Output files
    const char filename[]="output.txt";
    const char filenameN[]="ensambleN.txt";
    const char filenamex[]="ensamblex.txt";
    const char filenamec[]="parameters.txt"; //Just to print out all the parameters in the same folder at the end
    unsigned int seed; //Seed of the random number generator
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	int dummy, enne; //Dummy is a dummy index needed for small loops, enne is taking care (in case) of  how many times is rand bigger than interval
	bool checkiffixated=true; //This variable will check if a population fixated or not. It is set to true in the beginning just to get rid of the annoying warning.
	bool checkifMmax=false; //This variable will check if M is M_max or not. It is set to false in the beginning just to get rid of the annoying warning.
	int remainingsteps,anotherdummyindex; //I am going to use this to compute how many time steps I have left!  
    
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
	i=fread (&seed, sizeof (seed), 1, pfile); //I added the i= ... just to not be bothered anymore by the warnings!
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	//**********************************
    
    if(cons.choice==1){
		file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
		file<<"#Results for the simulation reproducing the propagule with"<<endl;
		file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		file<<left<<setw(12)<<"#Time"<<setw(12)<<"N"<<setw(15)<<"x"<<setw(12)<<"M"<<endl;
		myprint2(Nc,Nd,t,M,file);
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation reproducing the propagule with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		fileN<<"#In the form of N[t][m]"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation reproducing the propagule with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		filex<<"#In the form of x[t][m]"<<endl;
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
	}
	else if(cons.choice==2){
		file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
		file<<"#Results for the simulation reproducing the random splitting with"<<endl;
		file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		file<<left<<setw(12)<<"#Time"<<setw(12)<<"N"<<setw(15)<<"x"<<setw(12)<<"M"<<endl;
		myprint2(Nc,Nd,t,M,file);
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation reproducing the random splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		fileN<<"#In the form of N[t][m]"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation reproducing the random splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		filex<<"#In the form of x[t][m]"<<endl;
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
	}
	else{
		file.open(filename,ios::out|ios::trunc); //Open the output's file and print the results for time=0
		file<<"#Results for the simulation reproducing the deterministic splitting with"<<endl;
		file<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		file<<left<<setw(12)<<"#Time"<<setw(12)<<"N"<<setw(15)<<"x"<<setw(12)<<"M"<<endl;
		myprint2(Nc,Nd,t,M,file);
		fileN.open(filenameN,ios::out|ios::trunc); //Open the N's file 
		fileN<<"#Results for the simulation reproducing the deterministic splitting with"<<endl;
		fileN<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		fileN<<"#In the form of N[t][m]"<<endl;
		filex.open(filenamex,ios::out|ios::trunc); //Open the x's file and print the results for time=0
		filex<<"##Results for the simulation reproducing the deterministic splitting with"<<endl;
		filex<<"# M_max="<<cons.M_max<<"  T="<<cons.T<<"  K="<<cons.K<<"  s="<<cons.s<<"  p="<<cons.p<<"  N0="<<cons.N0<<"  x0="<<cons.x0<<"  N_max="<<cons.N_max<<"  seed="<<seed<<endl;
		filex<<"#In the form of x[t][m]"<<endl;
		myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex);
	}
    
    //*****Start of the evolution***********
     
   do{ 
        rand=randlog(Gamma[4*M-1],r);//Samples the time at wich the next reaction happens;
        if(t+rand>cons.T){ //If the new time would be bigger than T, I simply print out everything I have to print and then break the loop
		   		oldt=cons.T-t; //Note that here oldt has a different use, I'm simply using this variable since I don't need it anymore
		   		enne=floor(oldt/cons.interval);
		   		for(dummy=1;dummy<(enne+1);dummy ++){ //Note that here I have < of enne + 1 and not <=, since here I need the time explicitely! The enne+1 is because I need it to multply
		   			myprint2(Nc,Nd,t+dummy*cons.interval,M,file); 
		  			myprintensamble2(Nc,Nd,t+dummy*cons.interval,cons.M_max,fileN,filex);
		   			}
		   		myprint2(Nc,Nd,cons.T,M,file); //Printing the last time at T!!
		  		myprintensamble2(Nc,Nd,cons.T,cons.M_max,fileN,filex);
		   		break; //Exit from the do loop
		   }
        if(rand>cons.interval){ //Here is to check if I have to reprint the old situation before update the system!
		  		enne=floor(rand/cons.interval);
		  		for(dummy=1;dummy<=enne;dummy ++){  //<=because I start from 1
		  			myprint2(Nc,Nd,t+dummy*cons.interval,M,file); 
		  			myprintensamble2(Nc,Nd,t+dummy*cons.interval,cons.M_max,fileN,filex);
		  			cout<<"The time is "<<t+dummy*cons.interval<<endl; //Just to check
		  			}
		  		rand=rand-cons.interval*enne;
		  		t=t+enne*cons.interval;
		}
		t=t+rand; //Update the time. I have to do it after I do the check for the time
		oldt=oldt+rand; //Update oldt
        
		rand=gsl_rng_uniform(r)*Gamma[4*M-1]; //Generates the random number to choose the reaction!
		l=search(Gamma,4*M,rand); //Finds the reaction
        
		m=updateN(Nc, Nd,x,l); //Updates the variables at time i and returns the cell where the reaction happened
        
		if(check(Nc, Nd, cons, m)==true){ //Of course I need to check if I have to split the cell or not
			M=createcell(M, m, Nc, Nd, x, Gamma, G, cons, r); 	//Here I do everything, I create the cell, I update the cells and then update the Gamma and G
			//cout<<endl<<endl<<"Now in the main: First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[1]+Nd[1]<<" bacteria"<<endl<<endl;
		}
		else{ //Of course if no cell splits, I just update the G and the Gamma, print and then sample for another reaction
		//cout<<" la reazione e' avvenuta nella cellula m con Nc[m]:  "<<m<<"  "<<Nc[m]<<endl;
		updateG(G,Gamma,m,Nc,Nd,x,cons,4*M); //Updates the G and the Gamma
		}
        
        
		if(oldt>=cons.interval){ //Checks whether I have to print or not
			myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
			myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
        	oldt=oldt -cons.interval; //Subract by oldt the value of interval to start counting again
        	cout<<"The time is "<<t<<endl; //Just to check
		}
		
		checkiffixated=tochekciffixated(x,M,cons.M_max); //Here I check if I have to interrupt the do-while cycle or not!	
			
			if(M<cons.M_max){ //Here I check for the other condition of the do-while, cycle
				checkifMmax=false;
			}
			else{
				checkifMmax=true;
			}		
				
		
	}while((checkifMmax==false)&&(checkiffixated==false)); //I break the do while if M becomes M_max or all the groups in one population fixated!
	
	if((checkiffixated==true)||(checkifMmax==true)){ //If I stopped before TMAX I have to print for the rest of time steps the values of N and x
		remainingsteps=ceil((cons.T-t)/cons.interval);
		for(anotherdummyindex=0;anotherdummyindex<remainingsteps;anotherdummyindex++){
			t=t+cons.interval;
			myprint2(Nc,Nd,t,M,file); //Printing the results on file fast. To create a picture
			myprintensamble2(Nc,Nd,t,cons.M_max,fileN,filex); //Printing the results on file ensamble; to create the movie
		}
	}
	else{ //Otherwise I have to start again with a bigger T
		filex<<"Abort everything and start again with a bigger time!!!!!!!!!!!!!!!!";
		fileN<<"Abort everything and start again with a bigger time!!!!!!!!!!!!!!!!";
	}	
    
    file.close(); //Closing the files of output!
    filex.close();
    fileN.close();
    filec.open(filenamec,ios::out|ios::trunc); //Now I just print all the parameters to a file!
    printparamnoloop(filec,cons);
    filec.close();
    
    return 0;
}
