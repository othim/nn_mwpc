
-----------
Description
-----------
nn_mwpc is a code for computing nucleon-nucleon observables from different 
potential models. The main purpose why the code was created was to make a 
potential class that can include arbitrary potential terms, which makes it easy
to check different power-coutings. The code will also be able to perturbatively
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

-----------
DEFINITIONS
-----------

A COUPLED channel in np-scattering is a channel where there are four possible
states that can mix.


-----------
CONVENTIONS
-----------

POTENTIAL

LS-EQUATION

OBSERVABLES


------------------
COMPILING THE CODE
------------------

EXTERNAL LIBRARIES



----------------
RUNNING THE CODE
---------------

C++


PYTHON

--------------------------
Compiling to python module
--------------------------


