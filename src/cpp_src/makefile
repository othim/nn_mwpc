

CC = g++
CFLAGS2 = -Wall -g
CFLAGS  = -Wall -g -I/Users/toliver/wigxjpf-1.11/inc/ -I/usr/local/Cellar/gsl/2.7/include/gsl
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/

#CFLAGS = -O3
#If you want to profile with gprof you need to link with that library 

main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm  

run: main 
	./main

test.o: test.cpp pot_nn_mwpc.h
	$(CC) $(CFLAGS) -c test.cpp

pot_nn_mwpc.o: pot_nn_mwpc.h
	$(CC) $(CFLAGS) -c pot_nn_mwpc.cpp

Term.o: Term.h
	$(CC) $(CFLAGS) -c Term.cpp

LS_Solver.o: LS_Solver.h
	$(CC) $(CFLAGS) -c LS_Solver.cpp

quantum_states.o: quantum_states.h
	$(CC) $(CFLAGS) -c quantum_states.cpp

clean:
	rm main pot_nn_mwpc.o test.o Term.o LS_Solver.o quantum_states.o
