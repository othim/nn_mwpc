-----------
Description
-----------
nn_mwpc is a code for computing nucleon-nucleon observables from different 
potential models. The main purpose why the code was created was to make a 
potential class that can include arbitrary potential terms, which makes it easy
to check different power-countings. The code will also be able to perturbatively
include higher order potential contributions in the Born approximation.

The flow of the code is as follows.

- Build your potential by specifying what terms it should include

- Create a list of the possible quantum scattering channels

- Pass this information to the Lippmann Schwinger solver and get phase shifts.

- Pass the phase shifts to the code that computes the scattering observable of
  your choice.

The code is a C++ code. If you want to use all of its features, you have to
write this code in C++. However there is a python interface that is built for
easy access to certain key features that are pre-specified. To make more
features available from python you have to write the code in C++ and then
bind the calls to python via the interface. 

If you have any questions or suggestions of improvement you are very welcome to
open an issue on git or contact me 
(Oliver Thim, email: oliver.thim@chalmers.se).

-----------
DEFINITIONS
-----------

A COUPLED channel in np-scattering is a channel where there are four possible
states that can mix. E.g. the 3P0 channel is considered uncoupled since it 
cannot mix with any other L-value. For J>0 the the S=3 channel can always have
both L=J-1 and L=J+1.

-----------
CONVENTIONS
-----------

POTENTIAL

The potential is computed in a partial wave basis with the normalization 
<p',p> = \delta(p'-p)/p^2. There potential can precompute and save all matrix
elements that are not dependent on the on-shell energy.

LS-EQUATION

The LS-equatoin is written in the same partial wave basis as the potential.
This menas that there is no \pi/2 in front of the integral which appear in 
another widely used convention. This will affect the \rho_T parameter 
(see. eg Quantum Mechanics II - A second Course in Quantum Theory by Rubin H.
Landau) for a discussion. NOTE that this book uses the OTHER convention of
normalization of the parital wave states which produces a \pi/2 in the partial
wave LS eqiuation.

The LS-equation is solved using discretization of the momentum states
\int dp p^2 -> \sum_i w_i p_i^2. (see Rubin H. Landau Ch. 18)

OBSERVABLES

Observables are calculated from phase shifts by first computing the M-matrix 
and then computing the corresponding Saclay amplitudes which are then converted 
to observables. For the moment the trace computation directly from the M-matrix
is not implemented. 

------------------
COMPILING THE CODE
------------------

EXTERNAL LIBRARIES NEEDED
- gsl
- pybind11
- intel MKL (oneapi)

COMPILING

1. The first thing you need to do is to compile the external fortran code for the
nijmegen potential. This is done by running 'make cleanall' && 'make' in the 
nijm-fort-working-class directory. The object files 'pnijm.o' and 
'nijmegen_interface.o' then needs to be copied to 
the src/cpp_src directory to be accessed when compiling the C++ code. You can
see in the make file in src/cpp_src that it assumes that the fortran object
files are in that directory.

2. The C++ code is compiled by running the appropriate make command as defined 
in the make file in src/cpp_src depending on what target you want to have. 
e.g. 'make obs'.


!! UPDATE THIS !!
If you link to MKL instead of GSL you need some additional tricks. This is 
worth it because you can expect a speedup in the solution of LS-equation by a
factor 2-3 depending on the system. If you compile with the MKL flag in the 
makefile you will likely get errors that there are some missing shared library
files. Sometimes even worse, you will just get segfaults... The solution is a 
combination of what follows depending on the system you use.

    1. You need to locate the directory where MKL is installed. Sometimes the 
    enviroment variable $MKLROOT is set automatically, and sometimes not. In
    that directory under env/ or bin/ there is a shell script named
    vars.sh or setvars.sh. You need to run it by 'source setvars.sh'. This 
    will set some required enviroment variables. 

    2. If you still get errors you might need to export some enviroment variables
    manually. On the subatom computers you need to run:
    
    - 'export LD_LIBRARY_PATH=/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/'
    - 'export LD_PRELOAD=/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2: \
    /net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so: \
    /net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so.2: \
    /net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/net/opt/intel/lib/intel64/libiomp5.so'

    don't ask me why.....

    3. You must tell OpenMP and MKL how many threads they should use by running.
    The number of MKL threads should be one and the number of openMP threads
    depends on the machine and what you run. More is not always better!
    - 'export MKL_NUM_THREADS=1'
    - 'export OMP_NUM_THREADS=16'

    4. By doing all these steps you should get a factor 2-3 speedup on the 
    SUBATOM computers. Note that the export and setting of enviroment variables 
    needs to be done each time you go into a new terminal session.

    5. On Tetralith you just need to run the setvars.sh script, set the thread
    variables and run:
    'export LD_LIBRARY_PATH=<conda_env_dir>/.conda/envs/nn-mwpc-env:<gsl_dir>/gsl/lib'

!!
Note that the linking arguments might need a change since they are specific to
where your libraris are installed on your system. There is also different arguments
if you compile on Mac or Linux (check the makefile). If you compile the code on
a subatom computer you could probabily access my (Oliver Thim) libraries and 
the changes to the makefile should be minimal.

TESTING

AUTOMATED:

There is some autometed tests that are currently implemented. These can be run 
by running 'bash run_tests.sh' in this directory. The bash file will call 
./obs with different arguments to test different things. Note that it is likely
that some tests at the phase shift level will be marked as not passed. This 
can be convention dependent. There are currently come coupled channels that will
not pass, but if you run './obs phase not_test' and look at the ouput you see
that the code is indeed correct. 

NON AUTOMATED:

In the file compute_observable.cpp there are some different functionalities
acessed by giving different command line arguments. If you have compiled with
'make obs' you can run the ./obs excecutable. You have to give TWO arguments.
You can test phase shifts by runnning ./obs phase X, or some differnent observables
by running ./obs DSG X, or ./obs AYY X etc. Here if X = 'test' the code will 
run in test mode and provide less output. If X \neq 'test' this will not happen.

This will run the code with the
nijmegen1 potential and check the results agains data from nn-on-line saved in 
the data/ directory. When the observable data is generated by the ./obs program
it can be plotted with the python script plot_obs.py  by running
'python3 plot_obs.py <OBS_NAME>' in the src/python_src directory. Note that
the exact values of the observables and phase shifts are dependent on the
constants in the Constants.h file.

----------------------------
Compiling to a python module
----------------------------

If the makefile is set up correctly you just run 'bash install.sh' in the
nn_mwpc directory. This will build the code and install the python library
with pip. Then it is just to import it by 'import nn_mwpc' in python. Note that 
the library is sensitive to the python version used so make sure you are in the
conda environment nn-mwpc-env.

There is not so much documentation on the python interface. I recomend that 
one opens the interface files 'src/cpp_src/pybind_interfac.h and 
src/cpp_src/pybind_interface.cpp' and read the comments. There is also a file 
in this directory: 'module_test.py' that tests the library that you could run
after 'bash install.sh' to make sure that everything loads and runs fine. 
The code and comments in 'module_test.py' should give some guidance how to use 
the code by providing some examples.

-----------------------------
RUNNING THE CODE AND EXAMPLES
-----------------------------

C++


PYTHON

