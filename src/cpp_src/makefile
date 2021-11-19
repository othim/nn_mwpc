

CC = g++
CFLAGS2 = -Wall -ggdb -std=c++11
CFLAGS  = -Wall -ggdb -std=c++11 -I/Users/toliver/wigxjpf-1.11/inc/ -I/usr/local/Cellar/gsl/2.7/include/gsl
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/ 

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
NN_MWPC_FILES = compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

#CFLAGS = -O3
#If you want to profile with gprof you need to link with that library 

obs: compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) $(NN_MWPC_FILES) $(NIJM_FILES) -lgsl -lcblas -lgslcblas -lwigxjpf -lm -lgfortran 

main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm -lgfortran

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

pot_ext.o: pot_ext.h
	$(CC) $(CFLAGS) -c pot_ext.cpp


clean:
	rm main pot_nn_mwpc.o test.o Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o pot_ext.o
