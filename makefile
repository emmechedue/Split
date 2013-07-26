compile: Main.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

single_loop: Single_loop.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Single_loop.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

clean: 
	rm *~
	rm ./Headers/*~
	rm *.txt
	
deterministic: Deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
deterministic_single_loop: Single_loop_deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Single_loop_deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
prova: Prova.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Prova.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
noselection: Noselection.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers Noselection.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
numeric_check: Numeric_check.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Numeric_check.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
stops_at_m: Stops_at_M.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Stops_at_M.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
