
# Choose sytem
# ------------
SYSTEM = MAC
#SYSTEM = LINUX
# ------------
# Note that yo might need to remove the -fopenmp flag on mac in the so: target.

# Compiler
CC = g++

# For mac
ifeq ($(SYSTEM), MAC)
	
CFLAGS2 = -Wall -ggdb -p -std=c++11 -O3 -stdlib=libc++
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/  \
		   -L/usr/local/opt/llvm/lib

# For SO
CFLAGS = -Wall -fPIC -ggdb -p -O3 -std=c++11 -stdlib=libc++ \
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl -I/usr/local/opt/llvm/include

# This flag is needed on mac...
UNDEF = -undefined dynamic_lookup
endif

# For Linux
ifeq ($(SYSTEM), LINUX)

CFLAGS2 = -Wall -g -O3 -pg -std=c++11
LDFLAGS =  -L/net/home/toliver/wigxjpf-1.11/lib/ -L/net/home/toliver/gsl/lib 

# For SO
CFLAGS = -Wall -fopenmp -O3 -fPIC -std=c++11 \
		 -I/net/home/toliver/miniconda3/include/python3.8 \
		 -I/net/home/toliver/miniconda3/envs/nn-mwpc-env/include \
		 -I/net/home/toliver/wigxjpf-1.11/inc \
		 -I/net/home/toliver/phd/nn_mwpc/src/cpp_src/gsl-2.7/gsl

UNDEF = 
endif

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
NN_MWPC_FILES = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
				quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

# If you want to profile with gprof you need to link with that library 
obs: $(NN_MWPC_FILES)
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) $(NN_MWPC_FILES) $(NIJM_FILES) -lgsl -lgslcblas -lm  -lwigxjpf -lgfortran -fopenmp 

# OLD
main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm -lgfortran

run: main
	./main

FILES_SO = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
		   quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared $(UNDEF) $(LDFLAGS)  $(NIJM_FILES) $(FILES_SO) -o nn_mwpc.so -lgsl -lgslcblas -lm  -lwigxjpf -lgfortran

pybind_interface.o: pybind_interface.h 
	$(CC) $(CFLAGS) $(PY_BIND) -c pybind_interface.cpp

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
	rm main nn_mwpc.so pybind_interface.o compute_obs.o pot_nn_mwpc.o test.o \
		Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o \
	   	pot_ext.o compute_obs.o
