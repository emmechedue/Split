#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<Constants.h>

using namespace std;

/*
Constants::Constants(){
	N0=4; //Initial number of bacteria in the cell
	x0=0.5;//Initial fraction of cooperators in the cell
	T=2800.; //Time when the simulation stops
	interval=0.5; //Time step for which I print my results in fast
	s=0.38; //Selection's strenght
	c=1.; //Cost in the original fitness
	b=3.; //Benefit in the original fitness
	p=10.; //Cooperators advantage
	K=100.; //Carrying capacity
	N_max=80; //The number of bacteria in the cell s.t. the cell splits
	M_max=1000; //The maximum number of cells
	N_loop=200; //The number of times I iterate
	choice=2; //It's 1 if I want the propagule model and it is 2 if I want the random splitting model. It is 3 if I am choosing a deterministic splitting
	fitness=2; //It's 1 if I am using the original fitnesses (the one from the paper), is 2 if I am using the approximated fitnesses
	ts=0.0001; //The timestep in the semi-deterministic model. The timestep of the "integration" over N
}

Constants::~Constants(){}



*/
//This reads from a file:
Constants::Constants(){ //Note that name must be the entire path; i.e. "./config.conf"
	
	char line[256];
	int linenum=0;
	int count=0, M=15; //M is the amount of parameters I have to give, count will range from 0 to M-1
	double vector[M]; //will store the M parameters
	FILE *pfile;
	
	//pfile = fopen ("./config.conf" , "r");
	//pfile= fopen("/project/theorie/s/Stefano.Duca/Analysis/Prog/config.conf", "r"); //Here I have to put the folder where the config file will be!
	pfile= fopen("/project/theorie/s/Stefano.Duca/Analysis/I_need_this/config.conf", "r"); //Here I have to put the folder where the config file will be!
	
	while(fgets(line, 256, pfile) != NULL)
	{
		    

		    linenum++;
		    if(line[0] == '#') {continue;} //I'm going to the next line without reading and incrementing count

		    sscanf(line, "%*s %lf", &vector[count]);
			count ++;
	}

	N0=vector[0]; //Initial number of bacteria in the cell
	x0=vector[1];//Initial fraction of cooperators in the cell
	T=vector[2]; //Time when the simulation stops
	interval=vector[3]; //Time step for which I print my results in fast
	s=vector[4]; //Selection's strenght
	c=vector[5]; //Cost
	b=vector[6]; //Benefit
	p=vector[7]; //Cooperators advantage
	K=vector[8]; //Carrying capacity
	N_max=vector[9]; //The number of bacteria in the cell s.t. the cell splits
	M_max=vector[10]; //The maximum number of cells
	N_loop=vector[11]; //The number of times I iterate
	choice=vector[12]; //Is 1 if I want the propagule model and it is 2 if I want the random splitting model. It is 3 if I am choosing a deterministic splitting
	fitness=vector[13]; //Is 1 if I am using the original fitnesses (the one from the paper), is 2 if I am using the approximated fitnesses. Is 3 if I am imposing <f>=1, use ONLY when I have deterministic inner dynamics.
	ts=vector[14]; //The timestep in the semi-deterministic model. The timestep of the "integration" over N
}

Constants::~Constants(){}


