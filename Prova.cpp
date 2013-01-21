#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>

using namespace std;

/*void print(Constants consta){
	double y;
	
	cout<<consta.T<<endl;
	
	y=consta.T/consta.b;
	cout<<y<<endl;
}
	

int main(){
	Constants consta;
	
	print(consta);
	cout<<consta.T<<endl;
	
	return 0;
}*/

void function(int *N,int *M){
	*N=54;
	*M=102;
}

int main(){
	
	int N,M;
	
	N=0;
	M=0;
	function(&N,&M);
	cout<<N<<endl<<M<<endl;
	return 0;
}
