#pragma once
#include<cstdio>
#include<cstdlib>

class Constants{
public:
	Constants(); //Default constructor and destructor
	~Constants(); 
	
	int N0; //Initial number of bacteria in the cell
	double x0;//Initial fraction of cooperators in the cell
	double T; //Time when the simulation stops
	double interval; //Time step for which I print my results in fast
	double s; //Selection's strenght
	double c; //Cost in the original fitness
	double b; //Benefit in the original fitness
	double p; //Cooperators advantage
	double K; //Carrying capacity
	int N_max; //The number of bacteria in the cell s.t. the cell splits
	int M_max; //The maximum number of cells
	int N_loop; //The number of times I iterate
	int choice; //Is 1 if I want the propagule model and it is 2 if I want the random splitting model
	int fitness; //Is 1 if I am using the original fitnesses (the one from the paper), is 2 if I am using the approximated fitnesses
};


