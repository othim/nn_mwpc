/*
   LS_Solver.h
   This file is part of the nn_mwpc project.

   The LS_Solver class is a solver for the Lippmann-Schwinger equation
   "T = V + VGT". The main thing that needs to be specified is the
   potential, V, in the form of a matrix. The solver solves the LS-equation
   for a nucleon-nucleon system in a partial wave basis, and thus the output is
   the T-matrix in a given partial wave channel.

   The code also have functionality to solve the LS equation for a range of
   parameters (most importantly maximum angular momentum j_max) which allows for
   calculation of different spin-observebles, differential and total crossection.

   Oliver Thim 2021-09 --
   Department of Physics, Chalmers
*/

class LS_Solver
{
private:

public:

};
