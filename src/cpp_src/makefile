
# Choose sytem
# ------------
#SYSTEM = MAC
SYSTEM = SUBATOM
#SYSTEM=TETRALITH
# ------------

# Choose blas library
BLAS_LIB = MKL
#BLAS_LIB = GSL

# Compiler
CC = g++
#CC = icpc

# Optimization
OPT = -O3
#OPT = -O3 # -Ofast can have some optimazation that deviate from c standard
# For Mac
# -------
ifeq ($(SYSTEM), MAC)

CFLAGS2 = -Wall -std=c++11 -stdlib=libc++ $(OPT)
LDFLAGS =  -undefined dynamic_lookup -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.3.0/lib/gcc/11/  \
		   -L/usr/local/opt/llvm/lib 

ifeq ($(BLAS_LIB), MKL)
LDFLAGS += -L/opt/intel/oneapi/mkl/2021.4.0/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

CFLAGS = -Wall -fPIC -std=c++11 -stdlib=libc++ $(OPT) \
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl -I/usr/local/opt/llvm/include
endif 

# For SUBATOM
# -----------
ifeq ($(SYSTEM), SUBATOM)

CFLAGS2 = -Wall -std=c++11 $(OPT)
LDFLAGS =  -L/net/home/toliver/wigxjpf-1.11/lib/ -L/net/home/toliver/gsl-2.7/lib -L/usr/local/lib -fopenmp

ifeq ($(BLAS_LIB), MKL)

LDFLAGS  += -L/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64 \
			-Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

CFLAGS = -Wall -fopenmp -fPIC -std=c++11 $(OPT) \
		 -I/net/home/toliver/miniconda3/include/python3.8 \
		 -I/net/home/toliver/miniconda3/envs/nn-mwpc-env/include \
		 -I/net/home/toliver/wigxjpf-1.11/inc \
		 -I/net/home/toliver/gsl-2.7/gsl
endif 


# Tetralith
# ---------
ifeq ($(SYSTEM),TETRALITH)

WIGXJPF=/home/x_olith/wigxjpf-1.11

CFLAGS2 = -Wall -std=c++11 $(OPT)
#LDFLAGS = -L$(WIGXJPF)/lib -L$(MKLROOT)/lib/intel64 -L/home/x_olith/gsl/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl
LDFLAGS = -L$(WIGXJPF)/lib -L/home/x_olith/gsl/lib -fopenmp

#LDFLAGS +=  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#LDFLAGS +=  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
ifeq ($(BLAS_LIB), MKL)
LDFLAGS += -Nmkl -qmkl=parallel
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif
#LDFLAGS = -L$(WIGXJPF)/lib -L/home/x_olith/gsl/lib

CFLAGS =-Wall -fopenmp $(OPT) -fPIC -std=c++11 -I$(MKLROOT)/include -I/home/x_olith/gsl/include/gsl -I$(WIGXJPF)/inc
#CFLAGS += -I/software/sse/easybuild/prefix/software/Python/3.8.2-GCCcore-6.4.0-nsc1/include/python3.8 -I/software/sse/easybuild/prefix/software/pybind11/2.4.3-GCCcore-6.4.0-Python-3.8.2-nsc1/include/python3.8 -I/home/x_olith/.local/include/python3.8 
CFLAGS += -I/home/x_olith/.conda/envs/nn-mwpc-env/include
CFLAGS += -I/home/x_olith/.conda/envs/nn-mwpc-env/include/python3.8
CFLAGS += -DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include

endif

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
CDBONN_FILES = cdbonn_interface.o cdbonn.o # Compiled separately from Fortran code
NN_MWPC_FILES = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
				quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

obs: $(NN_MWPC_FILES)
	$(CC) $(CFLAGS2) -o obs $(LDFLAGS) $(NN_MWPC_FILES) $(NIJM_FILES) $(CDBONN_FILES) -lgsl -lm  -lwigxjpf -lgfortran

FILES_SO = pybind_interface.o compute_obs.o pot_nn_mwpc.o Term.o \
		   quantum_states.o LS_Solver.o scattering.o physics_helpers.o pot_ext.o

so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared $(LDFLAGS) $(CDBONN_FILES) $(NIJM_FILES) $(FILES_SO) -o nn_mwpc.so -lgsl -lm  -lwigxjpf -lgfortran 

pybind_interface.o: pybind_interface.cpp
	$(CC) $(CFLAGS) $(PY_BIND) $(FLAGS) -c pybind_interface.cpp

compute_obs.o: compute_obs.cpp
	$(CC) $(CFLAGS) $(FLAGS) -c compute_obs.cpp

pot_nn_mwpc.o: pot_nn_mwpc.h
	$(CC) $(CFLAGS) $(FLAGS) -c pot_nn_mwpc.cpp

Term.o: Term.h
	$(CC) $(CFLAGS) $(FLAGS) -c Term.cpp

LS_Solver.o: LS_Solver.h
	$(CC) $(CFLAGS) $(FLAGS) -c LS_Solver.cpp

quantum_states.o: quantum_states.h
	$(CC) $(CFLAGS) $(FLAGS) -c quantum_states.cpp

scattering.o: scattering.h
	$(CC) $(CFLAGS) $(FLAGS) -c scattering.cpp

physics_helpers.o: physics_helpers.h 
	$(CC) $(CFLAGS) $(FLAGS) -c physics_helpers.cpp	

pot_ext.o: pot_ext.h
	$(CC) $(CFLAGS) $(FLAGS) -c pot_ext.cpp

clean:
	rm obs nn_mwpc.so pybind_interface.o compute_obs.o pot_nn_mwpc.o \
		Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o \
	   	pot_ext.o
