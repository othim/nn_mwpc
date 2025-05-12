# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----------------- Choose system, BLAS library and compiler ------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Choose sytem
# ------------
# The code can be compiled on a number of systems, which requires different 
# includes and flags. Here you can add your own system.

SYSTEM = MAC
#SYSTEM = SUBATOM
#SYSTEM = SUBFY
#SYSTEM = TETRALITH
# ------------

# Choose blas library (GSL is the easy choice, but slightly slower)
# -------------------
#BLAS_LIB = MKL
BLAS_LIB = GSL
# -------------------

# Compiler (when using MKL, icpc is preferable)
# --------
CC = g++
#CC = icpc
# --------

# Optimization
# ------------
OPT = -O3
# ------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# For Mac
# -----------------------------------------------------------------------------
ifeq ($(SYSTEM), MAC)

GSL_PATH = /usr/local/Cellar/gsl
WIG_PATH = /Users/toliver/wigxjpf-1.11
ANACONDA_PATH = /opt/anaconda3

CFLAGS = -Wall -fPIC -std=c++11 -stdlib=libc++ $(OPT) \
		 -I$(ANACONDA_PATH)/envs/nn-mwpc-env/include \
		 -I$(ANACONDA_PATH)/include/python3.8 -I$(WIG_PATH)/inc \
		 -I$(GSL_PATH)/2.7/include/gsl -I/usr/local/opt/llvm/include

LDFLAGS =  -undefined dynamic_lookup -L$(WIG_PATH)/lib  -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/12  \
		   -L/usr/local/opt/llvm/lib 

ifeq ($(BLAS_LIB), MKL)
LDFLAGS += -L/opt/intel/oneapi/mkl/2021.4.0/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

endif 

# -----------------------------------------------------------------------------
# For SUBATOM (linux system)
# -----------------------------------------------------------------------------
ifeq ($(SYSTEM), SUBATOM)

GSL_PATH = /net/home/toliver/gsl-2.7
WIG_PATH = /net/home/toliver/wigxjpf-1.11


CFLAGS = -Wall -fopenmp -fPIC -std=c++11 $(OPT) \
		 -I$(WIG_PATH)/inc \
		 -I$(GSL_PATH)/gsl

LDFLAGS =  -L$(WIG_PATH)/lib -L$(GSL_PATH)/lib -L/usr/local/lib -fopenmp

ifeq ($(BLAS_LIB), MKL)
# You can use this flag for the icpc compiler
#LDFLAGS += -qmkl
LDFLAGS += -Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group 
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

endif 

# -----------------------------------------------------------------------------
# For SUBFY (linux system)
# -----------------------------------------------------------------------------
ifeq ($(SYSTEM), SUBFY)

GSL_PATH = /chalmers/users/toliver/gsl
WIG_PATH = /chalmers/users/toliver/wigxjpf-1.13
ANACONDA_PATH = /chalmers/users/toliver/miniconda3

CFLAGS = -Wall -fopenmp -fPIC -std=c++11 \
		 -I$(ANACONDA_PATH)/envs/env_test/include/python3.8 -I$(ANACONDA_PATH)/envs/env_test/include/python3.8 \
		 -I$(ANACONDA_PATH)/envs/env_test/include/python3.8 -I$(ANACONDA_PATH)/envs/env_test/lib/python3.8/site-packages/pybind11/include \
		 -I$(WIG_PATH)/inc \
		 -I$(GSL_PATH)/include/gsl

# Linker flags for the wigxjpf and gsl libraries
LDFLAGS =  -L$(WIG_PATH)/lib -L$(GSL_PATH)/lib \
		   -L/usr/local/lib -fopenmp

ifeq ($(BLAS_LIB), MKL) # OneAPI not installed on this system
LDFLAGS += -Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group 
endif

ifeq ($(BLAS_LIB), GSL)
#LDFLAGS += -lgslcblas
endif

endif 


# -----------------------------------------------------------------------------
# Tetralith (cluster - not up to date)
# -----------------------------------------------------------------------------
ifeq ($(SYSTEM),TETRALITH)

WIGXJPF = /home/x_olith/wigxjpf-1.11

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
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Define all the files
# -----------------------------------------------------------------------------

NIJM_FILES = nijmegen_interface.o pnijm.o # Compiled from fortran code
CDBONN_FILES = cdbonn_interface.o cdbonn.o # Compiled separately from Fortran code
IDAHO_N3LO_FILES = idaho-n3lo-interface.o idaho-n3lo.o # Compiled separately from Fortran code

# nn_mwpc files
FILES_SO = pybind_interface.o Term.o potential.o\
		   quantum_states.o LS_Solver.o scattering.o physics_helpers.o potential_ext.o \
		   born_approx.o pybind_interface_dwb.o potential_mwpc.o \
		   pre_def_potentials.o
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Compile the code to a shared library.
# This the .so file can then be compiled into a python module
# or linked to another C++ program. NOTE that the order of the 
# linking flags might matter...
so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared  $(LDFLAGS) $(CDBONN_FILES) $(NIJM_FILES) $(IDAHO_N3LO_FILES) $(FILES_SO) -o libnn_mwpc.so -lgsl -lgslcblas -lm -lwigxjpf -lgfortran 

# PY_BIND is a flag that is set in nn_mwpc/install.sh to make sure
# that the correct includes are included for building the 
# python module.
pybind_interface.o: pybind_interface.cpp
	$(CC) $(CFLAGS) $(PY_BIND) -c pybind_interface.cpp

potential.o: potential.cpp
	$(CC) $(CFLAGS) $(PY_BIND) -c potential.cpp

compute_obs.o: compute_obs.cpp
	$(CC) $(CFLAGS) -c compute_obs.cpp

test_mwpc_dwb_interface.o: test_mwpc_dwb_interface.cpp
	$(CC) $(CFLAGS) -c test_mwpc_dwb_interface.cpp

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

potential_ext.o: potential_ext.h
	$(CC) $(CFLAGS) -c potential_ext.cpp

born_approx.o: born_approx.h
	$(CC) $(CFLAGS) -c born_approx.cpp

pybind_interface_dwb.o: pybind_interface_dwb.h
	$(CC) $(CFLAGS) -c pybind_interface_dwb.cpp

potential_mwpc.o: potential_mwpc.h
	$(CC) $(CFLAGS) -c potential_mwpc.cpp

pre_def_potentials.o: pre_def_potentials.h
	$(CC) $(CFLAGS) -c pre_def_potentials.cpp

clean:
	rm  libnn_mwpc.so pybind_interface.o compute_obs.o potential.o \
		Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o \
	   	potential_ext.o born_approx.o potential_mwpc.o pybind_interface_dwb.o *.tmp \
		pre_def_potentials.o
