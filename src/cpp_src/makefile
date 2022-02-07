
# Choose sytem
# ------------
#SYSTEM = MAC
#SYSTEM = LINUX
SYSTEM = TETRALITH

# ------------
# Note that yo might need to remove the -fopenmp flag on mac in the so: target.

# Compiler
#CC = g++
CC = icpc

# For mac
ifeq ($(SYSTEM), MAC)
	
CFLAGS2 = -Wall -ggdb -std=c++11 -O3 -stdlib=libc++
LDFLAGS =  -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/  \
		   -L/usr/local/opt/llvm/lib

# For SO
CFLAGS = -Wall -fPIC -ggdb -O3 -std=c++11 -stdlib=libc++ \
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl -I/usr/local/opt/llvm/include

# This flag is needed on mac...
UNDEF = -undefined dynamic_lookup
OPENMP =
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
OPENMP = -fopenmp
endif

ifeq ($(SYSTEM),TETRALITH)

WIGXJPF=/home/x_olith/wigxjpf-1.11

CFLAGS2 = -Wall -O3 -std=c++11
#LDFLAGS = -L$(WIGXJPF)/lib -L$(MKLROOT)/lib/intel64 -L/home/x_olith/gsl/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl
LDFLAGS = -L$(WIGXJPF)/lib -L/home/x_olith/gsl/lib 

#LDFLAGS +=  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#LDFLAGS +=  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
LDFLAGS += -Nmkl -mkl=parallel

#LDFLAGS = -L$(WIGXJPF)/lib -L/home/x_olith/gsl/lib

CFLAGS =-Wall -fopenmp -O3 -fPIC -std=c++11 -I$(MKLROOT)/include -I/home/x_olith/gsl/include/gsl -I$(WIGXJPF)/inc
#CFLAGS += -I/software/sse/easybuild/prefix/software/Python/3.8.2-GCCcore-6.4.0-nsc1/include/python3.8 -I/software/sse/easybuild/prefix/software/pybind11/2.4.3-GCCcore-6.4.0-Python-3.8.2-nsc1/include/python3.8 -I/home/x_olith/.local/include/python3.8 
CFLAGS += -I/home/x_olith/.conda/envs/nn-mwpc-env/include
CFLAGS += -I/home/x_olith/.conda/envs/nn-mwpc-env/include/python3.8
CFLAGS += -DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include

OPENMP=-fopenmp
endif

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
NN_MWPC_FILES = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
				quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

# If you want to profile with gprof you need to link with that library 
obs: $(NN_MWPC_FILES)
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) $(NN_MWPC_FILES) $(NIJM_FILES) -lgsl -lm  -lwigxjpf -lgfortran $(OPENMP)

# OLD
main: test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o
	$(CC) $(CFLAGS2) -o main $(LDFLAGS) test.o pot_nn_mwpc.o Term.o quantum_states.o LS_Solver.o scattering.o physics_helpers.o -lgsl -lcblas -lgslcblas -lwigxjpf -lm -lgfortran

run: main
	./main

FILES_SO = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
		   quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared $(UNDEF) $(LDFLAGS)  $(NIJM_FILES) $(FILES_SO) -o nn_mwpc.so -lgsl -lm  -lwigxjpf -lgfortran $(OPENMP)

pybind_interface.o: pybind_interface.cpp
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
