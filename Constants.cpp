#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<Constants.h>

using namespace std;

Constants::Constants(){
	N0=4; //Initial number of bacteria in the cell
	x0=0.5;//Initial fraction of cooperators in the cell
	T=3; //Time when the simulation stops
	interval=0.05; //Time step for which I print my results in fast
	intervalens=0.001; //Time step for which I print my results in ensamble.txt NOTE THAT: in Split I'm using only interval
	s=0.05; //Selection's strenght
	p=10.; //Cooperators advantage
	K=100.; //Carrying capacity
	N_max=80; //The number of bacteria in the cell s.t. the cell splits
	M_max=1000; //The maximum number of cells
	N_loop=100; //The number of times I iterate
}

Constants::~Constants(){}

//In the future, substitute this constructor with something that reads from a file
