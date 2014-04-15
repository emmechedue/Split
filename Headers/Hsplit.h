#pragma once
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<Hevolve.h>
#include<gsl/gsl_rng.h>
#include<Constants.h>
#include<gsl_randist.h>

using namespace std;

// Note the function with the 2 at the end of the name are the ones without killing events, for the paper


/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
void fill1(int *C, int *D, double x, int N0, gsl_rng *r){ //This one gives C (the numbers of the cooperators in the cell) according to Bin(x,N0) and D as N0-C. So this is the propagule model
	int n;
	
	n=gsl_ran_binomial(r,x,N0);
	*C=n;
	*D=N0-n;
	return;
}

void fill2(int *C, int *D, double Nc, double Nd, gsl_rng *r){//This one gives C (the numbers of the cooperators in the cell) according to Bin(0.5,Nc) and D according to Bin(0.5,Nd), where Nc and Nd are the # of cooperators and defectors in the cell that just splitted. So this is the random splitting model
	int nd,nc;
	
	if(Nc<=0){
		nc=0;
	}
	else{
		nc=gsl_ran_binomial(r,0.5,Nc);
	}
	if (Nd<=0){
		nd=0;
	}
	else{
		nd=gsl_ran_binomial(r,0.5,Nd);
	}
	*C=nc;
	*D=nd;
	
	return;
}

void fill3(int n, int m, double *x, double *Nc, double *Nd, gsl_rng *r){//This one does the deterministic splitting, it just splits in two equal parts. It updates automatically the value of both cells. AGAIN REMEMBER THAT M IS THE OLD CELL AND N IS THE NEW CELL!!! If both Nc and Nd are odd, I puth the spare ones both in the same cell (to try to keep things balanced for small numbe rof agents)
	int temp1,temp2;
	int dice=-1; //Just to initialize it for something dummy, it will never be used unitialized
	
	
	
	temp1=(int) Nc[m];
	temp2=(int) Nd[m];
	if(((temp1%2)!=0)||((temp2%2)!=0)){ //I sample in which cell I want to put the spare agent/agents
		dice=gsl_rng_uniform_int(r,2);
	}
	Nc[n]=temp1/2;
	Nc[m]=temp1/2;
	if((temp1%2)!=0){
		if(dice==0){
			Nc[m]++;
		}
		else{
			Nc[n]++;
		}
	}
	
	Nd[n]=temp2/2;
	Nd[m]=temp2/2;
	if((temp2%2)!=0){
		if(dice==0){
			Nd[m]++;
		}
		else{
			Nd[n]++;
		}
	}
	x[n]=Nc[n]/(Nc[n]+Nd[n]);
	x[m]=Nc[m]/(Nc[m]+Nd[m]);
	
	return ;		
}
	
void fillcells(int n, int m, double *Nc, double *Nd, double *x, Constants cons, gsl_rng *r){ //n is the cell that I want to fill, m is the cell from which I take the parameters, choice can take values 1 (for the first method) and 2 for the second method
	int C,D;
	
	/*//Here I choose which splitting I want to do
	if(cons.choice==1){ 
		fill1(&C,&D,x[m],cons.N0,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
		Nc[n]=C;
		Nd[n]=D;
		x[n]=Nc[n]/(Nc[n]+Nd[n]);
		fill1(&C,&D,x[m],cons.N0,r);
		Nc[m]=C;
		Nd[m]=D;
		x[m]=Nc[m]/(Nc[m]+Nd[m]);
		
	}
	else{ 
		do{ //Just to be sure that I am not creating an empty cell!!
			fill2(&C,&D,Nc[m],Nd[m],r);  // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
		}while((C==0)&&(D==0));
		Nc[n]=C;
		Nd[n]=D;
		x[n]=Nc[n]/(Nc[n]+Nd[n]);
		Nc[m]=Nc[m]-C; //The new Nc[m] is the old Nc[m] - the cooperators in the new cell => I'm conserving the total number of bacteria
		Nd[m]=Nd[m]-D;
		x[m]=Nc[m]/(Nc[m]+Nd[m]);
	}*/
	switch(cons.choice){
		case 1:
			fill1(&C,&D,x[m],cons.N0,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
			Nc[n]=C;
			Nd[n]=D;
			x[n]=Nc[n]/(Nc[n]+Nd[n]);
			fill1(&C,&D,x[m],cons.N0,r);
			Nc[m]=C;
			Nd[m]=D;
			x[m]=Nc[m]/(Nc[m]+Nd[m]);
			break;
		case 2:
			do{ //Just to be sure that I am not creating an empty cell!!
				fill2(&C,&D,Nc[m],Nd[m],r);  // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
			}while((C==0)&&(D==0));
			Nc[n]=C;
			Nd[n]=D;
			x[n]=Nc[n]/(Nc[n]+Nd[n]);
			Nc[m]=Nc[m]-C; //The new Nc[m] is the old Nc[m] - the cooperators in the new cell => I'm conserving the total number of bacteria
			Nd[m]=Nd[m]-D;
			x[m]=Nc[m]/(Nc[m]+Nd[m]);
			break;
		case 3:
			fill3(n,m,x,Nc,Nd,r); //Here I do everything in one line!
			break;
		default:
			cout<<"ERROR IN THE CHOICE OF THE MODEL!!"<<endl;
			exit(101);
	}
	
	return;
}
	

int createcell(int M, int m,double *Nc, double *Nd, double *x, double *Gamma, double **G, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. the function returns the new value of M
	int n; //This is the index of one of the two new cells
	double rand;
	
	//********************Determine n**********
	if(M>=cons.M_max){
		rand=gsl_rng_uniform(r)*(cons.M_max-1); //Generate a uniform random number between 0 and M_max-1
		n=ceil(rand)-1; //I subtract -1 because the full array goes from 0 to M_max-1
		if(n==m){
			n=n+1; //If I extract just the cell I am already splitting, I take the next one
		} 
	}
	else{
		n=M; //not M+1 because of the index problem (it would be M+1-1)
		M=M+1;
	}
	
	//********creates the new cells and updates the Gamma and the G*****************
	fillcells(n,m,Nc,Nd,x,cons,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
	//cout<<"First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[n]+Nd[n]<<" bacteria"<<endl;
	if (n<m){ //  in the function I always have to have n<m
		updatebothG(G,Gamma,n,m,Nc,Nd,x,cons,4*M);
	}
	else{
		updatebothG(G,Gamma,m,n,Nc,Nd,x,cons,4*M);
	}
	return M;
}

int createcell2(int M, int m,double *Nc, double *Nd, double *x, double *Gamma, double **G, Constants cons, gsl_rng *r){ //M is the number of cells in total (it needs to go to M+1 or to stay M if M==M_max), m is the cell that is splitting, Nc,Nd and x are in the cell m. This functions splits the cell and checks wether I have to kill a cell or not. It accepts Nc, Nd and x because in this way I can use the same function to split in both ways. the function returns the new value of M
	int n; //This is the index of one of the two new cells
	double rand;
	
	//********************Determine n**********
	if(M>=cons.M_max){
		rand=gsl_rng_uniform(r)*(cons.M_max-1); //Generate a uniform random number between 0 and M_max-1
		n=ceil(rand)-1; //I subtract -1 because the full array goes from 0 to M_max-1
		if(n==m){
			n=n+1; //If I extract just the cell I am already splitting, I take the next one
		} 
	}
	else{
		n=M; //not M+1 because of the index problem (it would be M+1-1)
		M=M+1;
	}
	
	//********creates the new cells and updates the Gamma and the G*****************
	fillcells(n,m,Nc,Nd,x,cons,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
	//cout<<"First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[n]+Nd[n]<<" bacteria"<<endl;
	if (n<m){ //  in the function I always have to have n<m
		updatebothG2(G,Gamma,n,m,Nc,Nd,x,cons,4*M);
	}
	else{
		updatebothG2(G,Gamma,m,n,Nc,Nd,x,cons,4*M);
	}
	return M;
}
	
int createcellmoran(int M, int m,double *Nc, double *Nd, double *x, double *Gamma, double **G, Constants cons, gsl_rng *r){ /*Here I just do what I do in createcell but I also allow for the possibility to kill the cell that I just sampled.
I achieve this by cheating: when m comes up, in "createcell" I choose m+1 instead (it makes sense, due to the use of ceil), here if m comes up, I just procede normally but before I store the content of m+1 in an array. After fillcells and before updatebothG, I copy those values in the m-th place and then I am done! This relies on the fact that the cells are indistinguishable! 
At the end of the day I will end up with n and m switched! Basically the new cell should end up at the place of the old one but instead I switch m+1 and m (hence n)*/

	int n; //This is the index of one of the two new cells
	bool copycheck=false; //I am going to use this to check if I have to do the copy thing or not!
	double xcopy,Nccopy,Ndcopy;
	
	//********************Determine n**********
	if(M>=cons.M_max){
		n= gsl_rng_uniform_int(r,cons.M_max); //Generate a uniform random number between 0 and M_max-1 !! Here I generate all integers with equal probabilties
		if(n==m){ //If n==m, I copy the values and then I proceed with n
			copycheck=true;
			xcopy=x[m];
			Nccopy=Nc[m];
			Ndcopy=Nd[m];			
			n=n+1; //Now I take the next one!
		} 
	}
	else{
		n=M; //not M+1 because of the index problem (it would be M+1-1)
		M=M+1;
	}
	
	//********creates the new cells and updates the Gamma and the G*****************
	fillcells(n,m,Nc,Nd,x,cons,r); // It's important that I first create the n-cell and then the m one, because to create the cell I need the parameters of the m-th cell
	//cout<<"First cell now has "<<Nc[m]+Nd[m]<<" bacteria and second cell now has "<<Nc[n]+Nd[n]<<" bacteria"<<endl;
	if(copycheck==true){ //If nedeed I recopy the values of the (m+1)th cell in the mth cell!
		x[m]=xcopy; 
		Nc[m]=Nccopy;
		Nd[m]=Ndcopy;
	}
	
	if (n<m){ //  in the function I always have to have n<m
		updatebothG(G,Gamma,n,m,Nc,Nd,x,cons,4*M);
	}
	else{
		updatebothG(G,Gamma,m,n,Nc,Nd,x,cons,4*M);
	}
	return M;
}
	
