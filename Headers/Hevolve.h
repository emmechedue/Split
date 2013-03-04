#pragma once
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Constants.h>
#include<iomanip>

using namespace std;

/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
    
/*double fcoop(double x, Constants cons){
	double y;
	
	y=1.+cons.s*((cons.b-cons.c)*x-cons.c*(1-x));
	return y;
}

double fdef(double x, Constants cons){
	double y;
	
	y=1.+cons.s*cons.b*x;
	return y;
}*/ //Instead of using those one, I'm using fcoop=1-s, fdef= 1

double fcoop(double x, Constants cons){
	double y;
	
	y=1.-cons.s;
	return y;
}

double fdef(double x, Constants cons){
	
	return 1;
}

double g(double x, Constants cons){
	double y;
	
	y=1.+cons.p*x;
	return y;
}

double d(double Nc, double Nd, Constants cons){
	double y;
	
	y=(Nc+Nd)/cons.K;
	return y;
}    

double faverage(double x, Constants cons){ //computes the <f> as in the paper in order to use normalized fitness!
	double y;
	
	y=x*fcoop(x,cons)+(1.-x)*fdef(x,cons);
	return y;
}

void initializeGamma(double **G, double *Gamma,double *Nc, double *Nd, double *x, Constants cons){
    int j;
    double average;
    
    average=faverage(x[0],cons);
    j=0; //I have to create the first gamma by hand due to Gamma[0]
    G[0][0]=Nc[0]*g(x[0],cons)*fcoop(x[0],cons)/average;
    Gamma[0]=G[0][0];
    j++;
    G[0][1]=Nc[0]*d(Nc[0],Nd[0],cons); 
    Gamma[j]=Gamma[j-1]+G[0][1];
    j++;
    G[0][2]=g(x[0],cons)*Nd[0]*fdef(x[0],cons)/average;
    Gamma[j]=Gamma[j-1]+G[0][2];
    j++;
    G[0][3]=Nd[0]*d(Nc[0],Nd[0],cons); 
    Gamma[j]=Gamma[j-1]+G[0][3];
    return;
}


	

void myprint2(double *Nc,double *Nd,double t,int M, ofstream& file){ //Prints the time, N average and x average as defined in the very first paper (i.e. <x>=Sum(Nc)/Sum(N)) and also the number of cells M
    double Av,Ntot;
    int i; 
    
    file<<left<<setw(12)<<t; //Prints the time
    Ntot=0;
    for(i=0;i<M;i++){ //Computes the the average of N and prints it
        Ntot=Ntot+Nc[i]+Nd[i];
    }  
    Av=Ntot/M;
    file<<left<<setw(12)<<Av;
    Av=0;
    for(i=0;i<M;i++){ //Computes the the average of x as defined on the very first paper and prints it
        Av=Av+Nc[i];
        }
    Av=Av/Ntot;
    file<<left<<setw(15)<<Av<<M<<endl;
    return;
}

void myprintensamble2(double *Nc,double *Nd,double t,int M, ofstream& fileN, ofstream& filex){ //Prints all x in a file in the form x[t,m] and N in another file in the form N[t,m]
	double y;
	int i;

	//file<<t<<"    "; //Prints the time
	for(i=0; i<M; i++){ //Prints Nc+Nd and x
		y=Nc[i]+Nd[i];
		fileN<<left<<setw(7)<<y;
		y=Nc[i]/y;
		filex<<left<<setw(10)<<y;
	}
	fileN<<endl;
	filex<<endl;
	return ;
}	

void printiterens(double *Nc,double *Nd,int M, ofstream& fileN, ofstream& filex){ //Prints  <N> and <x> average as defined in the very first paper (i.e. <x>=Sum(Nc)/Sum(N)) 
    double Av,Ntot;
    int i; 
    
    Ntot=0;
    for(i=0;i<M;i++){ //Computes the the average of N and prints it
        Ntot=Ntot+Nc[i]+Nd[i];
    }  
    Av=Ntot/M;
    fileN<<left<<setw(12)<<Av;
    Av=0;
    for(i=0;i<M;i++){ //Computes the the average of x as defined on the very first paper and prints it
        Av=Av+Nc[i];
        }
    Av=Av/Ntot;
    filex<<left<<setw(15)<<Av;
    return;
}

int search(double *Gamma, int M, double x){ //Binary search
    int a,b,l,result;
    bool check;
    a=0;
    b=M-1;
    //cout<<endl<<"Random number is "<<x<<endl;
    do{
        l=(a+b)/2;
        if(x<=Gamma[l]){
            if((x>=Gamma[l-1])&&(l>0)){
            	result=l;
		   		check=true;}
            else{
		         if(l>0){
		            b=l;
		            check=false;}
		         else{
		         	result=0;
		         	check=true;
		         }
            }
        }
        else{
            if(x<=Gamma[l+1]){
                result=l+1;
                check=true;
            }
            else{
                a=l;
                check=false;}
        }
    }while(check==false);
   return result; 
}

int updateN(double *Nc, double *Nd,double *x, int l){  //Updates the N; l is the chosen reation from G (the one given by search), t is the time step (the new one, is just the index of the for!)
    int m,k; //k is the occured reation, m is the cell where the change occurred and is returned by the function
    
    m=l/4;
    k=l%4;
    switch(k){ //Update the appropriate number of bacteria
        case 0:
            Nc[m]=Nc[m]+1;
            break;
        case 1:
            Nc[m]=Nc[m]-1;
            break;
        case 2:
            Nd[m]=Nd[m]+1;
            break;
        case 3:
            Nd[m]=Nd[m]-1;
            break;
        default:
            cout<<"Error in updateN"<<endl;
            exit(1);
            break;
    }
    x[m]=Nc[m]/(Nc[m]+Nd[m]); //Update the x array
    return m;
}

bool check(double *Nc, double *Nd, Constants cons, int m){ //Check if the cell m (where the reaction occurred) has to split or not
	double N;
	bool asd=false;
	
	N=Nc[m]+Nd[m];
	if(N>cons.N_max){ //Checks if N>N_max
		asd=true;
	}
	return asd;
}

void updateG(double **G,double *Gamma, int m, double *Nc, double *Nd, double *x, Constants cons, int emme){ //Update the array G[m][] with the new Nc, Nd and x and the array Gamma
    double old[4];
    double average,sum;//Sum saves the difference of the old G[m][] with the new one;
    int i,a;
    
    for( i=0; i<4;i++){ //Save the changes of G[][]
        old[i]=G[m][i];
    }
    average=faverage(x[m],cons);
    G[m][0]=Nc[m]*g(x[m],cons)*fcoop(x[m],cons)/average; //Updates the G[][]
    G[m][1]=Nc[m]*d(Nc[m],Nd[m],cons); 
    G[m][2]=g(x[m],cons)*Nd[m]*fdef(x[m],cons)/average;
    G[m][3]=Nd[m]*d(Nc[m],Nd[m],cons);
    sum=0;
    for(i=0;i<4;i++){ //Compute the change
        sum=sum+G[m][i]-old[i];
    }
    a=4*m;
    if(m==0){ //Update Gamma[4*m]; I need to do in this way due to m=0
        Gamma[0]=G[0][0];
    }
    else{
        Gamma[a]=Gamma[a-1]+G[m][0];
    }
    for(i=1;i<4;i++) //Update the part of the Gamma[i] due to m
    {
        Gamma[i+a]=Gamma[i+a-1]+G[m][i];
    }
    for(i=a+4;i<emme;i++){ //I think this way is better because I have to make less calls (instead of Nd, Nc, x I just call sum)
        Gamma[i]=Gamma[i]+sum;
    }
    //cout<<"The gammas are: "<<G[0][0]<<"  "<<G[0][1]<<"  "<<G[0][2]<<"  "<<G[0][3]<<"and gamma j is "<<Gamma[emme-1]<<endl;
    return;
}

void updatebothG(double **G,double *Gamma, int n,int m, double *Nc, double *Nd, double *x, Constants cons, int emme){ //Update the arrays G[m][] and G[n][] with the new Nc, Nd and x and the array Gamma. Note that it requires that n<m!!!
    double old[4];
    double average, sum;//Sum saves the difference of the old G[m][] with the new one;
    int i,a;
    
    for( i=0; i<4;i++){ //Save the changes of G[][]
        old[i]=G[n][i];
    }
    average=faverage(x[n],cons);
    G[n][0]=Nc[n]*g(x[n],cons)*fcoop(x[n],cons)/average; //Updates the G[][]
    G[n][1]=Nc[n]*d(Nc[n],Nd[n],cons); 
    G[n][2]=g(x[n],cons)*Nd[n]*fdef(x[n],cons)/average;
    G[n][3]=Nd[n]*d(Nc[n],Nd[n],cons);
    sum=0;
    for(i=0;i<4;i++){ //Compute the change
        sum=sum+G[n][i]-old[i];
    }
    a=4*n;
    if(n==0){ //Update Gamma[4*n]; I need to do in this way due to n=0
        Gamma[0]=G[0][0];
    }
    else{
        Gamma[a]=Gamma[a-1]+G[n][0];
    }
    for(i=1;i<4;i++) //Update the part of the Gamma[i] due to n
    {
        Gamma[i+a]=Gamma[i+a-1]+G[n][i];
    }
    for(i=a+4;i<4*m;i++){ //Here I basically update all the Gamma from 4*n to 4*m (the one that is basically unchanged)
        Gamma[i]=Gamma[i]+sum;
    }
    //Now i do the same for m
    for( i=0; i<4;i++){ //Save the changes of G[][]
        old[i]=G[m][i];
    }
    average=faverage(x[m],cons);
    G[m][0]=Nc[m]*g(x[m],cons)*fcoop(x[m],cons)/average; //Updates the G[][]
    G[m][1]=Nc[m]*d(Nc[m],Nd[m],cons); 
    G[m][2]=g(x[m],cons)*Nd[m]*fdef(x[m],cons)/average;
    G[m][3]=Nd[m]*d(Nc[m],Nd[m],cons);
    sum=0;
    for(i=0;i<4;i++){ //Compute the change
        sum=sum+G[m][i]-old[i];
    }
    a=4*m;
    Gamma[a]=Gamma[a-1]+G[m][0]; //Here I don't do the check for m==0 because I'm sure that m!=0
    for(i=1;i<4;i++) //Update the part of the Gamma[i] due to m
    {
        Gamma[i+a]=Gamma[i+a-1]+G[m][i];
    }
    for(i=a+4;i<emme;i++){ //Here I update the rest, from 4*m to emme
        Gamma[i]=Gamma[i]+sum;
    }
    //cout<<"The gammas are: "<<G[0][0]<<"  "<<G[0][1]<<"  "<<G[0][2]<<"  "<<G[0][3]<<"and gamma j is "<<Gamma[emme-1]<<endl;
    return;
}

double randlog(double Gamma, gsl_rng *r){ //Generate the random number according to the distribution -ln(random[0,1]/Gamma)
    double x,y;
    
    y=gsl_rng_uniform(r);
    x=-log(y)/Gamma;
    return x;
}


