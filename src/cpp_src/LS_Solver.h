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
#include "pot_nn_mwpc.h"
#include "gsl_matrix.h"
#include <cmath>

struct Phase_shifts_chn {double delta_p; double delta_m; double epsilon;};

class LS_Solver
{
private:
   Potential_mwpc* pot_V_;
   std::vector<qs:quantum_channel> channels_;

   double* p_grid_;
   double* w_grid_;
   std::size_t mom_grid_size_;
   double cutoff_Lambda_;

   unsigned int J_max_;

   double cutoff_Lambda_;
   bool cutoff_enabled_;

   bool relcorr_enabled_;

   double* setup_D_vector(double q_on_shell);
   gsl_matrix* setup_F_matrix(double q_on_shell, gsl_matrix* V_mtx);

public:

   LS_Solver();
   ~LS_Solver();

   // Returns an array of phase shifts in the convention: ...
   Phase_shifts_chn solve_in_chn(solve_in_chn(double q_on_shell, qs::quantum_channel chn, bool rel_correction, bool cutoff_on);

};
