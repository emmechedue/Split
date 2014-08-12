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
	
cutted: Main_cutted.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main_cutted.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
mmax: Main_stop_at_M_max.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h
	g++ -I /usr/include/gsl -I ./Headers Main_stop_at_M_max.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
moran: Moran.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Moran.h
	g++ -I /usr/include/gsl -I ./Headers Moran.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
moran_deterministic: Moran_deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Moran.h
	g++ -I /usr/include/gsl -I ./Headers Moran_deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
single_loop_stops_at_m: Single_loop_stops_at_Mmax.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Single_loop_stops_at_Mmax.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
single_loop_fully_deterministic: Single_loop_fully_deterministic.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Single_loop_fully_deterministic.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
	
no_death: No_death.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries No_death.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

mutation: Mutation.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Mutation.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

mutation_nokill: Main_stop_at_M_max_mutation.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Main_stop_at_M_max_mutation.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2

mutation_singleloop: Single_loop_mutation.cpp Constants.cpp ./Headers/Hsplit.h ./Headers/Hevolve.h ./Headers/Constants.h ./Headers/Deterministic.h
	g++ -I /usr/include/gsl -I ./Headers -I /project/theorie/s/Stefano.Duca/Libraries Single_loop_mutation.cpp Constants.cpp -lgsl -lgslcblas -lm -Wall -O2
