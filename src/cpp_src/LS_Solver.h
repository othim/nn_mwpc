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
#include "gsl_complex.h"
#include "gsl_complex_math.h"
#include "gsl_matrix_complex_double.h"
#include <math.h>
#include <complex>
#include "scattering.h"
#include "physics_helpers.h"

class LS_Solver
{
private:
    Potential_mwpc* pot_V_;

    double* p_grid_;
    double* w_grid_;
    std::size_t mom_grid_size_;

    bool finite_grid_;
    double finite_grid_max_;

    gsl_vector*  setup_D_vector(double q_on_shell, bool coupled, double mu);
    gsl_matrix* setup_F_matrix(bool coupled, gsl_vector* D_vector, 
        gsl_matrix* V_mtx);


public:

    LS_Solver(unsigned int mom_grid_size, double* p_grid, double* w_grid, 
            bool finite_grid);

    ~LS_Solver();
    void gauss_legendre_inf_mesh(unsigned int Numper_of_points, 
            double scale,double** p,double** w);

    // Returns an array of phase shifts in the convention: Stapp
    Phase_shifts_chn solve_in_chn_R(double T_lab, qs::quantum_channel chn, 
            gsl_matrix* pot_V_mtx);

    double* solve_in_chn_R_Relem(double T_lab, 
            qs::quantum_channel chn, gsl_matrix* pot_V_mtx);

    gsl_vector_complex* setup_D_vector_complex(double q_on_shell, bool coupled, 
            double mu);
    gsl_matrix_complex* setup_F_matrix_complex(bool coupled, 
            gsl_vector_complex* D_vector, gsl_matrix* V_mtx);
    Phase_shifts_chn solve_in_chn_T(double T_lab, qs::quantum_channel chn, 
            gsl_matrix* pot_V_mtx);
   
    std::complex<double>* solve_in_chn_T_Telem(double T_lab, 
            qs::quantum_channel chn, gsl_matrix* pot_V_mtx);

    std::complex<double>* T_matrix_from_R_matrix(double Rmm, 
            double Rmp, double Rpp, double mu, double q_on_shell);

    static void get_mu_q_on_shell(double T_lab, qs::quantum_channel chn, 
            double* mu, double* q_on_shell);
};

#endif
