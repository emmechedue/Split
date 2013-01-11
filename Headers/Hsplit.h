#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Hevolve.h>

using namespace std;

/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
    

int createcell(int M, int m,double *Nc, double *Nd double *x, double *Gamma, double **G, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. the function returns the new value of M
	int n; //This is the index of one of the two new cells
	double rand;
	
	if(M>=cons.M_max){
		rand=gsl_rng_uniform(r)*(cons.M_max-1); //Generate a uniform random number between 0 and M_max-1
		n=ceil(rand)-1; //I subtract -1 because the full array goes from 0 to M_max-1
		if(n==m){
			n=n+1;
		} //If I extract just the cell I am already splitting, I take the next one
	}
	else{
		n=M+1;
		M=M+1;
	}
	fillcell(n) 
	fillcell(m)
	updatenandm
	
		
