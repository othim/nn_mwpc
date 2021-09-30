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

#ifndef LS_SOLVER_H
#define LS_SOLVER_H

#include "pot_nn_mwpc.h"
#include "gsl_matrix.h"
#include "gsl_sf_trig.h"
#include "gsl_blas.h"
#include "gsl_vector.h"
#include "gsl_linalg.h"
#include "quantum_states.h"
#include "gsl_permutation.h"
#include "gsl_integration.h" // GL integration
#include <cmath>


class LS_Solver
{
private:
   Potential_mwpc* pot_V_;

   // TODO REMOVE NOT USED
   std::vector<qs::quantum_channel> channels_;

   double* p_grid_;
   double* w_grid_;
   std::size_t mom_grid_size_;

   unsigned int J_max_;

   double cutoff_Lambda_;
   bool cutoff_enabled_;

   bool relcorr_enabled_;

   gsl_vector*  setup_D_vector(double q_on_shell, bool coupled, double mu);
   gsl_matrix* setup_F_matrix(bool coupled, gsl_vector* D_vector, gsl_matrix* V_mtx);
  

public:

   LS_Solver(std::vector<qs::quantum_channel> channels, Potential_mwpc* pot_V, unsigned int mom_grid_size=100,
      double mom_grid_scale=100.0, bool cutoff_enabled_ = true, double cutoff_Lambda_ = 450.0, bool relcorr_enabled = true);

   ~LS_Solver();
   void gauss_legendre_inf_mesh(unsigned int Numper_of_points, double scale,double** p,double** w);

   // Returns an array of phase shifts in the convention: ...
   Phase_shifts_chn solve_in_chn(double T_lab, qs::quantum_channel chn, bool rel_correction);

};

#endif