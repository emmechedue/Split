compile: Main.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
		

single_loop: Single_loop.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Single_loop.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

clean: 
	rm *.txt
	rm *~
	rm ./Headers/*~
	
deterministic: Deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
deterministic_single_loop: Single_loop_deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Single_loop_deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
