#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<iostream>
#include<Constants.h>

using namespace std;

void print(Constants consta){
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
}
