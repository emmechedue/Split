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
	double b; //Parameter b as in the paper
	double c; //Parameter c as in the paper
	double s; //Selection's strenght
	double p; //Cooperators advantage
	double K; //Carrying capacity
	int N_max; //The number of bacteria in the cell s.t. the cell splits
	int M_max; //The maximum number of cells
};

