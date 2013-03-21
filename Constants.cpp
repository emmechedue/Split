#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<Constants.h>

using namespace std;
/*
Constants::Constants(){
	N0=4; //Initial number of bacteria in the cell
	x0=0.5;//Initial fraction of cooperators in the cell
	T=1; //Time when the simulation stops
	interval=0.0001; //Time step for which I print my results in fast
	s=0.05; //Selection's strenght
	p=10.; //Cooperators advantage
	K=100.; //Carrying capacity
	N_max=80; //The number of bacteria in the cell s.t. the cell splits
	M_max=100; //The maximum number of cells
	N_loop=100; //The number of times I iterate
	choice=1; //Is 1 if I want the propagule model and it is 2 if I want the random splitting model
}

Constants::~Constants(){}


*/

//This reads from a file:
Constants::Constants(){ //Note that name must be the entire path; i.e. "./config.conf"
	
	char line[256];
	int linenum=0;
	int count=0, M=11; //M is the amount of parameters I have to give, count will range from 0 to M-1
	double vector[M]; //will store the M parameters
	FILE *pfile;
	
	pfile = fopen ("./config.conf" , "r");
	//pfile= fopen("/project/theorie/s/Stefano.Duca/Analysis/Prog/config.conf", "r"); //Here I have to put the folder where the config file will be!
	
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
	p=vector[5]; //Cooperators advantage
	K=vector[6]; //Carrying capacity
	N_max=vector[7]; //The number of bacteria in the cell s.t. the cell splits
	M_max=vector[8]; //The maximum number of cells
	N_loop=vector[9]; //The number of times I iterate
	choice=vector[10]; //Is 1 if I want the propagule model and it is 2 if I want the random splitting model
}

Constants::~Constants(){}


