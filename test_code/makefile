# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----------------- Choose system, BLAS library and compiler ------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Choose sytem
# ------------
SYSTEM = MAC
#SYSTEM = SUBATOM
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

CFLAGS = -I../src/cpp_src/ \
		 -Wall -fPIC -std=c++11 -stdlib=libc++ \
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 \
		 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl \
		 -I/usr/local/opt/llvm/include

LDPATH  = -L../src/cpp_src/
LDFLAGS = -lnn_mwpc
LDFLAGS = -lnn_mwpc -lgsl -lm -lwigxjpf -lgfortran -lgslcblas

LDPATH +=  -L/Users/toliver/wigxjpf-1.11/lib/ \
		   -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/12  \
		   -L/usr/local/opt/llvm/lib

test: test.o
	$(CC) -Wall -std=c++11 -o test $(LDPATH) test.o $(LDFLAGS)

obs: test_obs_nijmegen1.o
	$(CC) -Wall -std=c++11 -o obs $(LDPATH) test_obs_nijmegen1.o $(LDFLAGS)

WPC: test_LO_WPC.o
	$(CC) -Wall -std=c++11 -o WPC $(LDPATH) test_LO_WPC.o $(LDFLAGS)

test.o: test.cpp
	$(CC) $(CFLAGS) -c test.cpp

test_obs_nijmegen1.o: test_obs_nijmegen1.cpp
	$(CC) $(CFLAGS) -c test_obs_nijmegen1.cpp

test_LO_WPC.o: test_LO_WPC.cpp
	$(CC) $(CFLAGS) -c test_LO_WPC.cpp

clean:
	rm test.o test_obs_nijmegen1.o
