# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----------------- Choose system, BLAS library and compiler ------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Choose sytem
# ------------
#SYSTEM = MAC
#SYSTEM = SUBATOM
SYSTEM = SUBFY
#SYSTEM=TETRALITH
# ------------

# Choose blas library
# -------------------
#BLAS_LIB = MKL
BLAS_LIB = GSL
# -------------------

# Compiler
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

# For Mac
# -------
ifeq ($(SYSTEM), MAC)

CFLAGS2 = -Wall -std=c++11 -stdlib=libc++ $(OPT)
LDFLAGS =  -undefined dynamic_lookup -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/12  \
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

# Linker flags for the wigxjpf and gsl libraries
LDFLAGS =  -L/net/home/toliver/wigxjpf-1.11/lib/ -L/net/home/toliver/gsl-2.7/lib -L/usr/local/lib -fopenmp

ifeq ($(BLAS_LIB), MKL)

# Using these linkerflags you need to do export the following paths:
#export LD_LIBRARY_PATH=/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/
#export LD_PRELOAD=/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so.2:/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/net/opt/intel/lib/intel64/libiomp5.so
#LDFLAGS  += -L/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64 \
			-Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

# Using this flag you can do when using the icpc compiler
#
#LDFLAGS += -qmkl

LDFLAGS += -Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group \

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

# For SUBFY
# -----------
ifeq ($(SYSTEM), SUBFY)

CFLAGS2 = -Wall -std=c++11 $(OPT)

# Linker flags for the wigxjpf and gsl libraries
LDFLAGS =  -L/chalmers/users/toliver/wigxjpf-1.13/lib -L/chalmers/users/toliver/gsl/lib -L/usr/local/lib -fopenmp

ifeq ($(BLAS_LIB), MKL)


LDFLAGS += -Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group \

endif

ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

CFLAGS = -Wall -fopenmp -fPIC -std=c++11 $(OPT) \
		 -I$(HOME)/miniconda3/include/python3.8 \
		 -I/net/home/toliver/miniconda3/envs/nn-mwpc-env/include \
		 -I$(HOME)/wigxjpf-1.13/inc \
		 -I$(HOME)/gsl/include/gsl
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
IDAHO_N3LO_FILES = idaho-n3lo-interface.o idaho-n3lo.o # Compiled separately from Fortran code

FILES_SO = pybind_interface.o Term.o potential.o\
		   quantum_states.o LS_Solver.o scattering.o physics_helpers.o potential_ext.o \
		   born_approx.o pybind_interface_dwb.o potential_mwpc.o \
		   pre_def_potentials.o

so: $(FILES_SO)
	$(CC) -Wall -std=c++11 -shared $(LDFLAGS) $(CDBONN_FILES) $(NIJM_FILES) $(IDAHO_N3LO_FILES) $(FILES_SO) -o libnn_mwpc.so -lgsl -lm  -lwigxjpf -lgfortran 

pybind_interface.o: pybind_interface.cpp
	$(CC) $(CFLAGS) $(PY_BIND) $(FLAGS) -c pybind_interface.cpp

potential.o: potential.cpp
	$(CC) $(CFLAGS) $(PY_BIND) $(FLAGS) -c potential.cpp

compute_obs.o: compute_obs.cpp
	$(CC) $(CFLAGS) $(FLAGS) -c compute_obs.cpp

test_mwpc_dwb_interface.o: test_mwpc_dwb_interface.cpp
	$(CC) $(CFLAGS) $(FLAGS) -c test_mwpc_dwb_interface.cpp

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

potential_ext.o: potential_ext.h
	$(CC) $(CFLAGS) $(FLAGS) -c potential_ext.cpp

born_approx.o: born_approx.h
	$(CC) $(CFLAGS) $(FLAGS) -c born_approx.cpp

pybind_interface_dwb.o: pybind_interface_dwb.h
	$(CC) $(CFLAGS) $(FLAGS) -c pybind_interface_dwb.cpp

potential_mwpc.o: potential_mwpc.h
	$(CC) $(CFLAGS) $(FLAGS) -c potential_mwpc.cpp

pre_def_potentials.o: pre_def_potentials.h
	$(CC) $(CFLAGS) $(FLAGS) -c pre_def_potentials.cpp

clean:
	rm  libnn_mwpc.so pybind_interface.o compute_obs.o potential.o \
		Term.o LS_Solver.o quantum_states.o scattering.o physics_helpers.o \
	   	potential_ext.o born_approx.o potential_mwpc.o pybind_interface_dwb.o *.tmp \
		pre_def_potentials.o
