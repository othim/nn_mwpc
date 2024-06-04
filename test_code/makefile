# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----------------- Choose system, BLAS library and compiler ------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Choose sytem
# ------------
#SYSTEM = MAC
SYSTEM = SUBATOM
#SYSTEM=TETRALITH
# ------------

# Choose blas library
# -------------------
BLAS_LIB = MKL
#BLAS_LIB = GSL
# -------------------

# Compiler
# --------
#CC = g++
CC = icpc
# --------

# Optimization
# ------------
OPT = -O3
# ------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

ifeq ($(SYSTEM), MAC)
CFLAGS = -I../src/cpp_src/ \
		 -Wall -fPIC -std=c++11 -stdlib=libc++ $(OPT)\
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 \
		 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl \
		 -I/usr/local/opt/llvm/include

LDPATH  = -L../src/cpp_src/
LDFLAGS = -lnn_mwpc -lgsl -lm -lwigxjpf -lgfortran -lgslcblas

LDPATH +=  -L/Users/toliver/wigxjpf-1.11/lib/ \
		   -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/12  \
		   -L/usr/local/opt/llvm/lib
endif


ifeq ($(SYSTEM), SUBATOM)

CFLAGS = -Wall -fopenmp -fPIC -std=c++11 \
		 -I../src/cpp_src/\
		 -I/net/home/toliver/miniconda3/include/python3.8 \
		 -I/net/home/toliver/miniconda3/envs/nn-mwpc-env/include \
		 -I/net/home/toliver/wigxjpf-1.11/inc \
		 -I/net/home/toliver/gsl-2.7/gsl

#LDFLAGS = -lnn_mwpc -lgsl -lm -lwigxjpf -lgfortran -lgslcblas -fopenmp

#LDPATH = -L../src/cpp_src/\
		 -L/net/home/toliver/wigxjpf-1.11/lib/\
		 -L/net/home/toliver/gsl-2.7/lib\
		 -L/usr/local/lib\
		 -L./

LDFLAGS = -lnn_mwpc -lm  -fopenmp
LDPATH = -L/usr/local/lib -L./
endif

NPHASE: test_nijm_phase.o
	$(CC) -Wall -std=c++11 -o NPHASE $(LDPATH) test_nijm_phase.o $(LDFLAGS)

OBS: test_obs_nijmegen1.o
	$(CC) -Wall -std=c++11 -o OBS $(LDPATH) test_obs_nijmegen1.o $(LDFLAGS)

WPC: test_LO_WPC.o
	$(CC) -Wall -std=c++11 -o WPC $(LDPATH) test_LO_WPC.o $(LDFLAGS)

POT: test_potentials.o
	$(CC) -Wall -std=c++11 -o POT $(LDPATH) test_potentials.o $(LDFLAGS)

DWB: test_dwb.o
	$(CC) -Wall -std=c++11 -o DWB $(LDPATH) test_dwb.o $(LDFLAGS)

test_nijm_phase.o: test_nijm_phase.cpp
	$(CC) $(CFLAGS) -c test_nijm_phase.cpp

test_obs_nijmegen1.o: test_obs_nijmegen1.cpp
	$(CC) $(CFLAGS) -c test_obs_nijmegen1.cpp

test_LO_WPC.o: test_LO_WPC.cpp
	$(CC) $(CFLAGS) -c test_LO_WPC.cpp

test_potentials.o: test_potentials.cpp
	$(CC) $(CFLAGS) -c test_potentials.cpp

test_dwb.o: test_dwb.cpp
	$(CC) $(CFLAGS) -c test_dwb.cpp

clean:
	rm test_nijm_phase.o test_obs_nijmegen1.o test_potentials.o test_LO_WPC.o\
		test_dwb.o
