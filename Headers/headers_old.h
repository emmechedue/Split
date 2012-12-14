#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>

using namespace std;

/****************List of errors: ****************
    Error 1 is an error in the function upadateN
    */
    
    
double fcoop(double x, double b, double c, double s){
	double y;
	
	y=1.+s*((b-c)*x-c*(1-x));
	return y;
}

double fdef(double x, double b, double s){
	double y;
	
	y=1.+s*b*x;
	return y;
}

double g(double x, double p){
	double y;
	
	y=1+p*x;
	return y;
}

double d(double Nc, double Nd, double K){
	double y;
	
	y=(Nc+Nd)/K;
	return y;
}    

double faverage(double x, double b, double c, double s){ //computes the <f> as in the paper in order to use normalized fitness!
	double y;
	
	y=x*fcoop(x,b,c,s)+(1-x)*fdef(x,b,s);
	return y;
}

void initializeGamma(double **G, double *Gamma,int M, double *Nc, double *Nd, double *x, double p, double s, double K,double b, double c){
    int i=0,j;
    
    j=0; //I have to create the first gamma by hand due to Gamma[0]
    G[0][0]=Nc[i]*g(x[i],p)*fcoop(x[i],b,c,s)/faverage(x[i],b,c,s);
    Gamma[0]=G[0][0];
    j++;
    G[0][1]=Nc[i]*d(Nc[i],Nd[i],K); 
    Gamma[j]=Gamma[j-1]+G[0][1];
    j++;
    G[0][2]=g(x[i],p)*Nd[i]*fdef(x[i],b,s)/faverage(x[i],b,c,s);
    Gamma[j]=Gamma[j-1]+G[0][2];
    j++;
    G[0][3]=Nd[i]*d(Nc[i],Nd[i],K); 
    Gamma[j]=Gamma[j-1]+G[0][3];
    j++;
    //cout<<"The gammas are: "<<G[0][0]<<"  "<<G[0][1]<<"  "<<G[0][2]<<"  "<<G[0][3]<<"and gamma j is "<<Gamma[j-1]<<endl;
    for(i=1; i<M;i++){ //Create the Gammas; the order is G_(0->C), G_(c->0), G_(0->D), G_(D->0) and start back for the new cell
       G[i][0]=Nc[i]*g(x[i],p)*fcoop(x[i],b,c,s)/faverage(x[i],b,c,s);
       Gamma[j]=Gamma[j-1]+G[i][0];
       j++;
       G[i][1]=Nc[i]*d(Nc[i],Nd[i],K); 
       Gamma[j]=Gamma[j-1]+G[i][1];
       j++;
       G[i][2]=g(x[i],p)*Nd[i]*fdef(x[i],b,s)/faverage(x[i],b,c,s);
       Gamma[j]=Gamma[j-1]+G[i][2];
       j++;
       G[i][3]=Nd[i]*d(Nc[i],Nd[i],K);
       Gamma[j]=Gamma[j-1]+G[i][3];
       j++; 
       //cout<<"check to see if the foris doing something!!!!!!!!!";
    }//Here I also computed the overall Gamma
    return;
}


void myprint1(double *Nc,double *Nd,double *x,double t,int M, ofstream& file){ //Prints the time, N average and x average (i.e. <x>=Sum(x)/M)
    double Av;
    int i; 
    
    file<<t<<"   "; //Prints the time
    Av=0;
    for(i=0;i<M;i++){ //Computes the the average of N and prints it
        Av=Av+Nc[i]+Nd[i];
    }  
    Av=Av/M;
    file<<Av<<"   ";
    Av=0;
    for(i=0;i<M;i++){ //Computes the the average of x and prints it
        Av=Av+x[i];
    }  
    Av=Av/M;
    file<<Av<<endl;
    return;
}

void myprint2(double *Nc,double *Nd,double t,int M, ofstream& file){ //Prints the time, N average and x average as defined in the very first paper (i.e. <x>=Sum(Nc)/Sum(N))
    double Av,Ntot;
    int i; 
    
    file<<t<<"   "; //Prints the time
    Ntot=0;
    for(i=0;i<M;i++){ //Computes the the average of N and prints it
        Ntot=Ntot+Nc[i]+Nd[i];
    }  
    Av=Ntot/M;
    file<<Av<<"   ";
    Av=0;
    for(i=0;i<M;i++){ //Computes the the average of x as defined on the very first paper and prints it
        Av=Av+Nc[i];
        }
    Av=Av/Ntot;
    file<<Av<<endl;
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
            if(x>=Gamma[l-1]){result=l;
            check=true;}
            else{
                b=l;
                check=false;}
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
    int m,k; //m is the cell where the change occurred, l is the occured reation and is returned by the function
    
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

void updateG(double **G,double *Gamma, int m, double *Nc, double *Nd, double *x, double p, double s, double K,double b, double c, int emme){ //Update the array G[m][] with the new Nc, Nd and x and the array Gamma
    double old[4];
    double sum;//Sum saves the difference of the old G[m][] with the new one;
    int i,a;
    
    for( i=0; i<4;i++){ //Save the changes of G[][]
        old[i]=G[m][i];
    }
    G[m][0]=Nc[m]*g(x[m],p)*fcoop(x[m],b,c,s)/faverage(x[m],b,c,s); //Updates the G[][]
    G[m][1]=Nc[m]*d(Nc[m],Nd[m],K); 
    G[m][2]=g(x[m],p)*Nd[m]*fdef(x[m],b,s)/faverage(x[m],b,c,s);
    G[m][3]=Nd[m]*d(Nc[m],Nd[m],K);
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

double randlog(double Gamma, gsl_rng *r){ //Generate the random number according to the distribution -ln(random[0,1]/Gamma)
    double x,y;
    
    y=gsl_rng_uniform(r);
    x=-log(y)/Gamma;
    return x;
}


