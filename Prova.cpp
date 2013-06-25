#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<iostream>
#include<gsl/gsl_rng.h>
#include<Hevolve.h>
#include<Hsplit.h>
#include<Constants.h>
#include<Deterministic.h>

using namespace std;

int main(){
	Constants cons;
	double Nc[2],Nd[2],x[2];
	gsl_rng *r; //Pointer to the type of rng
	FILE *pfile; //file to read from /usr/urandom
	unsigned int seed; //Seed of the random number generator
	
	pfile = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (seed), 1, pfile);
	fclose(pfile);
	r = gsl_rng_alloc(gsl_rng_mt19937); //I'm using the "Mersenne Twister" generator!
	gsl_rng_set(r,seed); // Starting the generator
	
	Nc[0]=41.;
	Nd[0]=41.;
	Nc[1]=-1.;
	Nd[1]=-1.;
	fill3(1,0,x,Nc,Nd,r);
	cout<<Nc[0]<<"  "<<Nd[0]<<"  "<<Nc[1]<<"  "<<Nd[1]<<"  "<<x[0]<<"  "<<x[1]<<endl;
	
	return 0;
}
