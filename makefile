compile: Main.cpp Program.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Program.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall #-O2
		
clean: 
	rm *~
	rm ./Headers/*~
	rm output.txt
	
