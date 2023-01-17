# Description

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

# Definitions

A **coupled** channel in np-scattering is a channel where there are four possible
states that can mix. E.g. the $^3P_0$ channel is considered uncoupled since it 
cannot mix with any other L-value. For J>0 the the S=3 channel can always have
both L=J-1 and L=J+1.


# Conventions

## Potential and LS-equation

The potential is computed in a partial wave basis with the normalization 
$\langle p'|p \rangle = \delta(p'-p)/p^2$. The LS-equation for the potential
thus reads

```math
T^{js}_{l'l}(p', p) = V^{js}_{l'l}(p',p) \ + \sum_{l''}\int_0^\infty dk \ k^2 \ V^{js}_{l'l''}(p',k) \frac{m_N}{p^2-k^2+i\epsilon} T^{js}_{l''l}(k,p)$$,
```
i.e. without a factor of $2/\pi$ in front of the integral, which is another
common convention.

The partial wave decomposition is done according to the formulas in 
"MOMENTUM SPACE CALCULATIONS AND HELICITY FORMALISM IN NUCLEAR PHYSICS" by Erkelenz et al. (1971).
NOTE: that there is an error in eq. (4.22) which is pointed out in 
"Chiral effective field theory and nuclear forces" by Machleidt and Entem (2011).

We use the convention to add an additional minus sign w.r.t. the above calculation
on the potential elements that are off-diagonal in $l$. This added minus sign
is a basis convention that will be compensated in the calculation of the 
spin-scattering matrix.

The LS-equation can be solved either for the $T$-matrix or the $R$ matrix. 

To summarize: The potential that is given e.g. by the `pot_nn_mwpc.get_matrix()` is 
a potential that goes into the LS-equation above where the off-diagonal elements 
have an extra minus sign such that the M-matrix equation will have a $i^{l'-l}$
factor.

## Observables

Observables are calculated from phase shifts by first computing the M-matrix 
and then computing the corresponding Saclay amplitudes which are then converted 
to observables. With the convention to add the extra minus sign to the off-diagonal
elements give the following relation

```math
    \begin{align}
    &M^{s}_{m'_s m_s}(p,\theta_{cm},\phi) = \frac{\sqrt{4\pi}}{2ip} \sum_{j,l,l'} i^{l-l'} (2j+1)\sqrt{2l+1} \cdot \nonumber\\ &Y^{l'}_{m_s -m'_s}(\theta_{cm},\phi)
    \cdot \begin{pmatrix} l' & s & j  \\ m_s-m'_s & m'_s & -m_s \end{pmatrix} \cdot \begin{pmatrix} l & s & j \\ 0 & m_s& -m_s\end{pmatrix} \\
    &\cdot \left(S^{js}_{l'l}(p,p)-\delta_{l'l}\right)
    \end{align}
```

# COMPILING THE CODE


## EXTERNAL LIBRARIES NEEDED
- gsl (https://www.gnu.org/software/gsl/)
- pybind11 (https://pypi.org/project/pybind11/)
- intel MKL, (oneAPI)
- wigxjpf (http://fy.chalmers.se/subatom/wigxjpf/)

## COMPILING C++ CODE

1. The first thing you need to do is to compile the external fortran code for the
nijmegen potential and cdbonn potential. This is done by running
`$ bash compile_nijmegen.sh` and `$ bash compile_cdbonn.sh`. If this step fails 
you can try to perform the same tasks as the scripts manually.

2. Make sure that you go to the makefile in `src/cpp_src` and choose the correct
compiler, BLAS library and make sure that the paths to the external libraries
are changes to where you have installed them. Note that if you use MKL, you 
might need to source a setvars.sh file: `$ source intel/oneapi/setvars.sh`.

3. The C++ code is compiled by running the appropriate make command as defined 
in the make file in `src/cpp_src` depending on what target you want to have, e.g.
`make obs` or `make so`.

4. To set the number of allowed threads for MKL and openMP you set the 
the environment variables e.g.
```
$ export MKL_NUM_THREADS=1
$ export OMP_NUM_THREADS=16
```

Note that the linking arguments might need a change since they are specific to
where your libraris are installed on your system. There is also different arguments
if you compile on Mac or Linux (check the makefile). If you compile the code on
a subatom-computer you could probabily access my (Oliver Thim) libraries and 
the changes to the makefile should be minimal.

## COMPILING AND INSTALLING PYTHON MODULE

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



5. (Only on Tetralith) you just need to run the setvars.sh script, set the thread
variables and run:
```
$ export LD_LIBRARY_PATH=<conda_env_dir>/.conda/envs/nn-mwpc-env:<gsl_dir>/gsl/lib
```


# TESTING

## AUTOMATED

There is some autometed tests that are currently implemented. These can be run 
by running 'bash run_tests.sh' in this directory. The bash file will call 
./obs with different arguments to test different things. Note that it is likely
that some tests at the phase shift level will be marked as not passed. This 
can be convention dependent. There are currently come coupled channels that will
not pass, but if you run './obs phase not_test' and look at the ouput you see
that the code is indeed correct. 

## NON AUTOMATED

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


# RUNNING THE CODE AND EXAMPLES

C++


PYTHON

