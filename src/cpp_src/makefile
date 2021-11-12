

CC = g++
CFLAGS2 = -Wall -ggdb -std=c++11
CFLAGS  = -Wall -ggdb -std=c++11 -I/Users/toliver/wigxjpf-1.11/inc/ -I/usr/local/Cellar/gsl/2.7/include/gsl
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/

#CFLAGS = -O3
#If you want to profile with gprof you need to link with that library 

obs: compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm  

main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm  

run: main
	./main

compute_obs.o: compute_obs.cpp
	$(CC) $(CFLAGS) -c compute_obs.cpp

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

scattering.o: scattering.h
	$(CC) $(CFLAGS) -c scattering.cpp

physics_helpers.o: physics_helpers.h 
	$(CC) $(CFLAGS) -c physics_helpers.cpp	

clean:
	rm main pot_nn_mwpc.o test.o Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o
