#include<cstdio>
#include<cstdlib>

class Constants{
public:
	Constants(); //Default constructor and destructor
	~Constants(); 
	
	int N0=12; //Initial number of bacteria in the cell
	double x0=0.5;//Initial fraction of cooperators in the cell
	double T=45.; //Time when the simulation stops
	double interval=0.001; //Time step for which I print my results in fast
	double b=3.; //Parameter b as in the paper
	double c=1.; //Parameter c as in the paper
	double s=0.05; //Selection's strenght
	double p=10.; //Cooperators advantage
	double K=100.; //Carrying capacity
	int N_max=80; //The number of bacteria in the cell s.t. the cell splits
	int M_max=1000; //The maximum number of cells
}

//In the future add private function for reading all the constants from a conf file!
