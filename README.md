# Table of Contents

1. [Description](#description)
2. [Definitions](#definitions)
3. [Conventions](#conventions)
4. [Compiling the Code](#compiling-the-code)
5. [Testing the Code](#testing-the-code)
6. [Using the Python Module](#using-the-python-module)
6. [Appendix](#appendix)



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
write this code in C++. However, there is a python interface that is built for
easy access to certain key features that are pre-specified. To make more
features available from python you have to write the code in C++ and then
bind the calls to python via the interface. 

If you have any questions or suggestions for improvement you are very welcome to
open an issue on git or contact me 

(Oliver Thim, email: oliver.thim@chalmers.se).

# Definitions

- A **coupled** channel in np-scattering is a channel where there are four possible
states that can mix. e.g. the $^3P_0$ channel is considered uncoupled since it 
cannot mix with any other L-value. For J>0 the S=3 channel can always have
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
"Momentum Space Calculations and Helicity Formalism in Nuclear Physics" by Erkelenz et al. (1971).
NOTE: There is an error in eq. (4.22) which is pointed out in 
"Chiral effective field theory and nuclear forces" by Machleidt and Entem (2011).

We use the convention to add an additional minus sign w.r.t. the above calculation
on the potential elements that are off-diagonal in $l$. This added minus sign
is a basis convention that will be compensated in the calculation of the 
spin-scattering matrix.

The LS-equation can be solved either for the $T$-matrix or the $R$ matrix. The 
relation between the partial wave $S$-matrix and $T$-matrix in this basis reads

```math
     S^{js}_{l'l}(p',p) = \delta_{l'l} - 2\pi i \mu p \ T^{js}_{l'l}(p',p).
```

To summarize: The potential that is given e.g. by the `pot_nn_mwpc.get_matrix()` is 
a potential that goes into the LS-equation above where the off-diagonal elements 
have an extra minus sign such that the M-matrix equation will have a $i^{l-l'}$
factor.

## Observables

Observables are calculated from phase shifts by first computing the $M$-matrix 
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
The relation between the $M$-matrix elements and the Saclay amplitudes and observables
are found in "Formalism of nucleon-nucleon elastic scattering experiments" by Bystricky
et al. (1978).

The method `get_M_matrix_T(...)` in `src/cpp_src/scattering.h` will compute 
elements of the above-defined $M$-matrix. This is then used in the other 
methods in `scattering.h` to further compute observables.

# Compiling the Code


## External Libraries Needed
- gsl (https://www.gnu.org/software/gsl/)
- pybind11 (https://pypi.org/project/pybind11/)
- intel MKL, (oneAPI)
- wigxjpf (http://fy.chalmers.se/subatom/wigxjpf/)


## Setup
1. The first thing you need to do is to compile the external Fortran code for the
nijmegen1 and cdbonn potential. This is done by running
`$ bash compile_nijmegen.sh` and `$ bash compile_cdbonn.sh`. 

- If this step fails you can try to perform the same tasks as the scripts manually.

2. Make sure that you go to the makefile in `src/cpp_src` and choose the correct
compiler and BLAS library. Make sure that the paths to the external libraries
are changed to where you have installed them. Note that if you use MKL, you 
might need to source a setvars.sh file: `$ source /net/home/toliver/intel/oneapi/setvars.sh`.

3. If everything is setup correctly you should be able to go to `src/cpp_src` and run:
`$ make clean` and then `$ make` without errors.


### Compiling C++ Code (If you want to install and use the python module you can skip this section)

1. The C++ code is compiled by running the appropriate make command as defined 
in the makefile in `src/cpp_src`. 

- Always start with running `$ make clean` to remove any old object files.
- If you want to compile the code in order to make some tests from the 
`src/cpp_src/compute_observable.cpp` file you run `$ make obs` that will create
an executable named `obs`.
- There are two flags that can be used in order to compile the code
with non-default physics constants defined in `src/cpp_src/Constants.h`. This is
in order to match the benchmarks exactly. Not including this flag will set the 
default constant that is taken from PDG 2022. Examples:
    - `$ make obs FLAGS=-DANDREAS_CONST`
    - `$ make obs FLAGS=-DNIJM_CONST`
    - Defalut: `$ make obs`

2. The program is parallelized using OpenMP. 
The number of allowed threads for MKL and OpenMP are set prior to runtime using
the environment variables below.
```
$ export MKL_NUM_THREADS=1
$ export OMP_NUM_THREADS=16
```
- It is advised to keep the number of threads to MKL low and to OMP high for 
optimal performance. 
- The code is parallelized over the channels when solving 
for phase shifts, so if you just solve in one channel it might be beneficial 
to increase the number of threads given to MKL

3. You can test the program by running `$ make -j && ./obs phase nijm`, which 
will compute a bunch of phase shifts. You will most probably see that the test
does not pass. This is because the benchmarks are generated with slightly different
physics constants. Instead try `$ make -j FLAGS=-DNIJM_CONST && ./obs phase nijm`.
This should give that all tests pass, otherwise something is wrong.

### Compiling and Installing Python Module

1. Install and activate the conda environment that is defined in the file `environment.yml`
by running the commands 
- `conda env-create -f environment.yml`
- `conda activate nn-mwpc-env`

2. Make sure that you made steps 1-3 under Setup so that the makefile
is setup correctly.
    - Run `$ bash install.sh`, which is a shell script that will compile the 
    C++ code to a `.so` file and then making it into a python module with
    the use of pybind11.
    - If you run `$ conda list` you should find a module named `nn-mwpc`.

3. The program is parallelized using OpenMP. 
The number of allowed threads for MKL and OpenMP are set prior to runtime using
the environment variables below.
```
$ export MKL_NUM_THREADS=1
$ export OMP_NUM_THREADS=16
```

# Testing the code
If everything is correct, all tests should pass with no exceptions. If you are
using the python module make sure that you run both the test of the C++ code
and the python module.

## C++ code

- To run test on the C++ code, run `$ bash run_test.sh`. This will do a bunch of 
tests both with the Nijmegen potential and LO WPC potential. It will take a 
few minutes for these tests to run.

## Python module

- Make sure that you have installed the python module

- Run `$ python3 python_module/run_tests_on_module.py`


## (Non-automated tests - this can be skipped)

In the file compute_observable.cpp there are some different functionalities
accessed by giving different command line arguments. If you have compiled with
'make obs' you can run the ./obs executable. You have to give TWO arguments.
You can test phase shifts by running ./obs phase X, or some different observables
by running ./obs DSG X, or ./obs AYY X etc. Here if X = 'test' the code will 
run in test mode and provide less output. If X \neq 'test' this will not happen.

This will run the code with the
nijmegen1 potential and check the results against data from nn-on-line saved in 
the data/ directory. When the observable data is generated by the ./obs program
it can be plotted with the python script plot_obs.py  by running
'python3 plot_obs.py <OBS_NAME>' in the src/python_src directory. Note that
the exact values of the observables and phase shifts are dependent on the
constants in the Constants.h file.


# Using the Python Module

- The python module is made from the C++ class defined in `src/cpp_src/pybind_interface.h`.
The class definition contains a lot of comments that can be used to understand 
the functionalities of the class.

- If the module is installed you can import it by `import nn_mwpc` in python. 
Note that the library is sensitive to the python version used so make sure you are in the
conda environment `nn-mwpc-env`.

- There is a file `python_module/module_test.py` that contains examples of how to
compute different things using the module. This file also contains some useful 
comments. If the module is installed correctly you can run 
`$ python3 python_module/module_test.py`

# Appendix
## Misc notes
- (Only if you run on Tetralith) You need to run:
```
$ export LD_LIBRARY_PATH=<conda_env_dir>/.conda/envs/nn-mwpc-env:<gsl_dir>/gsl/lib
```
