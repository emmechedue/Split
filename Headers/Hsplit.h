#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Hevolve.h>

using namespace std;

/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
    
void createcell(int *M, int m, double x, double *Gamma, double **G){ //M is the number of cells in total (it needs to go to M+1), m is the cell that is splitting, x is the fraction of cooperators in the cell m. Note that the function is automatically updating M to M+1!!!!
	
