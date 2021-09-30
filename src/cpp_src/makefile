

CC = g++
CFLAGS = -Wall -g
#If you want to profile with gprof you need to link with that library 

main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o
	$(CC) $(CFLAGS) -o main test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o -lgsl -lcblas -lgslcblas -lm

run: main 
	./main

test.o: test.cpp pot_nn_mwpc.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c test.cpp

pot_nn_mwpc.o: pot_nn_mwpc.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c pot_nn_mwpc.cpp

Term.o: Term.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c Term.cpp

LS_Solver.o: LS_Solver.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c LS_Solver.cpp

quantum_states.o: quantum_states.h
	$(CC) $(CFLAGS) -c quantum_states.cpp

clean:
	rm main pot_nn_mwpc.o test.o Term.o LS_Solver.o quantum_states.o
