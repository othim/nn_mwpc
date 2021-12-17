

CC = g++


# For mac
CFLAGS2 = -Wall -std=c++11 -O3 -stdlib=libc++
#CFLAGS  = -Wall -std=c++11 -O3 -I/Users/toliver/wigxjpf-1.11/inc/ -I/usr/local/Cellar/gsl/2.7/include/gsl
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/  

# For SO
CFLAGS = -Wall -fPIC -O3 -std=c++11 -stdlib=libc++ -I/opt/anaconda3/envs/nn-mwpc-env/include -I/opt/anaconda3/include/python3.8 -I/Users/toliver/wigxjpf-1.11/inc/ -I/usr/local/Cellar/gsl/2.7/include/gsl


# For Linux
#CFLAGS2 = -Wall -g -pg -std=c++11
#CFLAGS  = -Wall -g -pg -std=c++11 -I/net/home/toliver/wigxjpf-1.11/inc/ -I/net/home/toliver/phd/nn_mwpc/src/cpp_src/gsl-2.7/gsl 
#LDFLAGS =  -L/net/home/toliver/wigxjpf-1.11/lib/ -L/net/home/toliver/gsl/lib # -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/ 

# For SO
#CFLAGS = -Wall -fPIC -std=c++11 -I/net/home/toliver/miniconda3/include/python3.8 -I/net/home/toliver/miniconda3/envs/nn-mwpc-env/include  -I/net/home/toliver/wigxjpf-1.11/inc/ -I/net/home/toliver/phd/nn_mwpc/src/cpp_src/gsl-2.7/gsl

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
NN_MWPC_FILES = compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

#CFLAGS = -O3
#If you want to profile with gprof you need to link with that library 

obs: compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) $(NN_MWPC_FILES) $(NIJM_FILES) -lgsl -lgslcblas -lm  -lwigxjpf -lgfortran -g  

main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm -lgfortran

run: main
	./main

FILES_SO = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared -undefined dynamic_lookup  $(LDFLAGS)  $(NIJM_FILES) $(FILES_SO) -o nn_mwpc.so -lgsl -lgslcblas -lm  -lwigxjpf -lgfortran

pybind_interface.o: pybind_interface.h pybind_interface.cpp
	$(CC) $(CFLAGS) -c pybind_interface.cpp

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
	rm main nn_mwpc.so pybind11_interface.o compute_obs.o pot_nn_mwpc.o test.o Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o pot_ext.o compute_obs.o
