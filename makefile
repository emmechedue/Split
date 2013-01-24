compile: Main.cpp Prova.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
		

single_loop: Single_loop.cpp Prova.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall #-O2

		
clean: 
	rm *~
	rm ./Headers/*~
	rm output.txt
	
