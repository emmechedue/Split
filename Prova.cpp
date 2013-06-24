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
	double Nc, Nd, N,x,t;
	int i;
	int temp;
	
	/*N=cons.N0;
	x=cons.x0;
	Nc=x*N;
	Nd=N-Nc;
	t=0;
	cout<<setprecision(5)<<left<<setw(10)<<"N"<<setw(10)<<"x"<<setw(10)<<"Nc"<<setw(10)<<"Nd"<<setw(10)<<"t"<<endl;
	cout<<setprecision(5)<<left<<setw(10)<<N<<setw(10)<<x<<setw(10)<<Nc<<setw(10)<<Nd<<setw(10)<<t<<endl;
	for(i=0;i<=3;i++){
		N=Nevolve(N, x, cons.ts, cons);
		x=xevolve(x,cons.ts,cons);
		Nc=x*N;
		Nd=N-Nc;
		t=t+cons.ts;
		cout<<setprecision(5)<<left<<setw(10)<<N<<setw(10)<<x<<setw(10)<<Nc<<setw(10)<<Nd<<setw(10)<<t<<endl;
	}
	N=63.59;
	x=0.6;
	computeNcNd(x,N,&Nc,&Nd);
	cout<<setprecision(5)<<left<<setw(10)<<N<<setw(10)<<x<<setw(10)<<Nc<<setw(10)<<Nd<<setw(10)<<t<<endl;*/
	Nc=81.;
	temp=(int) Nc;
	Nd=temp/2;
	cout<<Nd<<endl;
	return 0;
}
