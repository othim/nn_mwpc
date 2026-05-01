# Table of Contents

1. [Description](#description)
2. [Definitions](#definitions)
3. [Conventions](#conventions)
4. [Compiling the Code](#compiling-the-code)
5. [Using the Python Module](#using-the-python-module)
6. [Appendix](#appendix)

# Description

nn_mwpc is a code for computing nucleon-nucleon observables from different
potential models. The main purpose why the code was created was to make a
potential class that can include arbitrary potential terms, which makes it easy
to check different power countings. The code is also able to perturbatively
include higher order potential contributions when solving for the NN T-matrix.

The flow of the code is as follows.

- Build your potential by specifying what terms it should include

- Create a list of the possible quantum scattering channels

- Pass this information to the Lippmann Schwinger solver to get T-amplitudes,
  perturbatively or non-perturbatively.

- The amplitudes can then be used to compute: NN phase shifts and
  NN scattering observables.

The code is a C++ code. If you want to use all of its features, you have to
write this code in C++. However, there is a python interface that is built for
easy access to certain key features that are pre-specified. To make more
features available from python you have to write the code in C++ and then
bind the calls to python via the interface.

If you have any questions or suggestions for improvement you are very welcome to
open an issue on git or contact me

## Written by
Oliver Thim

Copyright (C) 2026 Oliver Thim

oliver.thim97@gmail.com

## License
GNU GENERAL PUBLIC LICENCE Version 3, 29 June 2007, see LICENSE.txt.
This program comes with ABSOLUTELY NO WARRANTY. It is the users responsibility
to verify that the the code is installed and used correctly.

- If you use the code in a research project, appropriate credit should be given
  to the author. This can be done by citing the relevant publications
  where this code is presented:

O. Thim, A. Ekström, and C. Forssén, Phys. Rev. C 109, 064001  (2024),

O. Thim, Chiral Effective Field Theory with Partly Perturbative Pions
Applied to the Few-Nucleon Sector, PhD Thesis, Chalmers, (2026).

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
on the potential elements that are off-diagonal in $l$. This added minus sign ($i^{l-l'}$)
is a basis convention that will be compensated in the calculation of the
spin-scattering matrix, (that is why the $i^{l-l'}$ factor appears again here).
This convention is sometimes known as the 'Machleidt
convention', and is the same convention as is used in the potential codes
for the nijmegen1 and cdbonn potentials that is included in the repo.

The LS-equation can be solved either for the $T$-matrix or the $R$ matrix. The
relation between the partial wave $S$-matrix and $T$-matrix in this basis reads

```math
     S^{js}_{l'l}(p',p) = \delta_{l'l} - 2\pi i \mu p \ T^{js}_{l'l}(p',p).
```

The nucleon mass $m_N = 2*\mu$.

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
- intel MKL, (optional)
- wigxjpf (http://fy.chalmers.se/subatom/wigxjpf/)

## Setup
1. The first thing you need to do is to compile the external Fortran code for the
nijmegen1 and cdbonn potential. This is done by running
`$ bash compile_nijmegen.sh` and `$ bash compile_cdbonn.sh`.

- If this step fails you can try to perform the same tasks as the scripts manually.

2.

- Make a basic conda environment from the `conda_environments/env_test.yml` file.
   Install pybind11 (in the same conda environment, it can be installed via both pip
   and conda I think.) Other python packages in the `conda_environments/envirnoment.yml`
   file are needed, but it is best to install them manually from your minimal
   environment env_test. I think that the python version needs to be 3.8.

- Install gsl and wigxjpf, and optionally intel MKL.

3. Go to the makefile in `src/cpp_src` change the include and link paths to
match the install paths of the external libraries on your system. Note that if you use MKL, you
might need to source a `setvars.sh` file e.g.: `$ source <install-path>/intel/oneapi/setvars.sh`.

4. If everything is setup correctly you should be able to go to `src/cpp_src` and run:
`$ make clean` and then `$ make` without errors. If you get errors, it is most likely
related to the various included packages. You might need to add additional flags
on your system etc. Sometimes you need to provide paths external libraries in the
`LD_PRELOAD` and `LD_LIBRARY_PATH` variables.

5. The program is parallelized using OpenMP.
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

6. Testing the code can be done by running `$ test_code/run_tests.sh`. This will compiles the
code in the `src/cpp_src` directory into an `.so` file and link it to various codes
in the directory `test_code`. (see the shell script). You will get various output
files that you can inspect.

## Compiling and Installing Python Module

1. Install and activate the conda environment that you have created, e.g.
- `$ conda activate env_test` (or make a better name)

2.
- Run `bash install.sh`, which is a shell script that will compile the
  C++ code to a `.so` file and then make it into a python module with
  the use of pybind11.

- If you run `$ conda list` you should find a module named `nn-mwpc`.

3. The program is parallelized using OpenMP.
The number of allowed threads for MKL and OpenMP are set prior to runtime using
the environment variables below.
```
$ export MKL_NUM_THREADS=1
$ export OMP_NUM_THREADS=16
```

4. To test that the module is installed correctly, you should be able to run
 ```
 $ python3 python_module/run_tests_on_module.py
 ```

# Using the Python Module

- If the module is installed you can import it by `import nn_mwpc` in python.

- The python module is made from the C++ classes defined in
`src/cpp_src/pybind_interface.h` and `src/cpp_src/pybind_interface_dwb.h`
The class definitions contains some comments that can be used to understand
the functionalities of the class. I also provide some example scripts.

- The file `python_module/module_test.py` contains examples of how to
compute different things using the nn_mwpc_interface module. This is used only
for non-perturbative computations. This file also contains some useful
comments.

- Example file for using the nn_mwpc_dwb_interface for perturbative computations
  is provided separately.

# Appendix
## Misc notes
- (Only if you run on Tetralith) You need to run:
```
$ export LD_LIBRARY_PATH=<conda_env_dir>/.conda/envs/nn-mwpc-env:<gsl_dir>/gsl/lib
```

- For some reason, the nijmegen and cdbonn potentials can not be used from
  the python module when `OMP_NUM_THREADS>1`. For using these potentials from
  the python interface, make sure that you run `$ export OMP_NUM_THREADS=1`
  before running.

- If you encounter MKL-related errors when running the python module you might
  try the following.
  - Use intel compilers when using MKL.
  - Check that the python package versions are correct in the conda environement.
  - You might try to update some package, e.g. numpy.
  - If there are specific library files that not are loaded correctly you can
    try to add them manually to `LD_LIBRARY_PATH` and `LD_PRELOAD`.
  - It is easiest to get the code running with GSL, so that is advised.
