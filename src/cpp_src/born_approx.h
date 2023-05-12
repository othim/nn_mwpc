/*
    born_approx.h
    This file is part of the nn_mwpc project.
    
    This header contains functions to compute amplitude of scattering 
    processes in Born- and distorted-wave Born approximation.

    Oliver Thim 2022-09 --
    Department of Physics, Chalmers
*/
#ifndef BORN_APPROX
#define BORN_APPROX

#include "quantum_states.h"
#include "physics_helpers.h"
#include "Constants.h"
#include "wigxjpf.h"
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "gsl_matrix.h"
#include "gsl_blas.h"
#include "gsl_linalg.h"
#include "gsl_complex.h"
#include "potential_mwpc.h"
#include "LS_Solver.h"
#include <string>
namespace dwba 
{
 
/*
 * This function compute the Born approximation to the scattering amplitude.
 * 
 * start_order: is the start position in the series, this is inculded
 * stop_order : is the highest order that is included, this is included
 * e.g. start_order=1,stop_order=2 returns (V*G0*V + V*G0*V*G0*V)
 * V    : is the potential matrix
 * G0   : is the free Greens function
 *
 * V and G0 has to be in the form that the momentum gid weights and momenta
 * are included in the definition of the potential and propagator to make 
 * the resoultion of identity integrals to be just ordinary matrix 
 * multiplications.
 *
 * start_order = 0,
 *
 * stop_order=0: T = V
 * stop_order=1: T = V + V*G0*V
 * stop_order=2: T = V + V*G0*V + V*G0*V*G0*V
 * ...
 */

gsl_matrix_complex* pw_T_BA(int start_order,int stop_order, gsl_matrix_complex* V, 
        gsl_matrix_complex* G0);



/* This function computes the distorted-wave Born approximation to a given 
 * order for the input parameters
 *
 * order: is the order that the approximation is computed at. 0,1,2.
 * T_I  : is the full T-matrix as given by the solution to
 *        T_I = V_I + V_I*G0*T_I.
 * V_I  : is the leading order potential
 * V_II : is the correction to the leading order potential
 * G0   : is the free Greens function
 *
 * T_I, V_I, V_II and G0 has to be in the same form as accepted into 
 * 'dwba::pw_T_BA()' i.e. including the momentum space factors and 
 * weights.
 *
 *
 * order=0: returns T_I
 * order=1: returns T_I + Omega^dagger_- V_II Omega_+
 * order=2: returns T_I + Omega^dagger_- (V_II + V_II*G1*V_II) Omega_+ 
 * ...
 *
 * where G1 is the full greens function, related to the free Greens function as
 * G1 = G0 + G0*T_I*G0
 */

gsl_matrix_complex* pw_T_DWBA(int order, gsl_matrix_complex* T_I, 
        gsl_matrix_complex* V_I, gsl_matrix_complex* V_II, 
        gsl_matrix_complex* G0);

/*
 * ****************************************************************************
 * Functions to compute the corrections to T-matrices in MWPC
 * ****************************************************************************
 */


/* This function computes the DWBA according to 
 * T_NLO = \Omega_m^\dagger V_NLO \Omega_p.
 *
 * T_I  : The full off shell leading order T-matrix
 * G0   : The free propagator
 * V_NLO: The correction potential
 *
 * Note that T_NLO computed by this function is the NLO correction, to get the 
 * full NLO amplitude the you must add T_LO + T_NLO.
 */
gsl_matrix_complex* pw_T_DWBA_PC_NLO(gsl_matrix_complex* T_I,
        gsl_matrix_complex* G0, gsl_matrix_complex* V_NLO);


/* This function computes the DWBA according to 
 * T_N2LO = \Omega_m^\dagger V_N2LO \Omega_p + \Omega_m^\dagger 
 * (V_NLO*G1*V_NLO) \Omega_p
 *
 * where G1 is the propagator for H_LO = H_0 + V_LO.
 *
 * T_I   : The full off shell leading order T-matrix
 * G0    : The free propagator
 * V_NLO : The correction potential at NLO
 * V_N2LO: The correction potential at N2LO
 *
 * Note that T_N2LO computed by this function is the N2LO correction, to get the 
 * full N2LO amplitude the you must add T_LO + T_NLO + T_N2LO.
 */
gsl_matrix_complex* pw_T_DWBA_PC_N2LO(gsl_matrix_complex* T_I,
        gsl_matrix_complex* G0, gsl_matrix_complex* V_NLO,
        gsl_matrix_complex* V_N2LO);
/*
 * This function computes the DWBA according to 
 * T_N3LO = ... see NN notes.
 *
 * T_I    : Full LO T-matrix
 * G0     : Free propagator
 * V_NLO  : Correction potential at NLO
 * V_N2LO : Correction potential at N2LO
 * V_N3LO : Correction potential at N3LO
 *
 * Note that this is just the N3LO correction to the T-matrix and that the full
 * sum off LO and correction need to be compute the get the full T-matrix.
 */
gsl_matrix_complex* pw_T_DWBA_PC_N3LO(gsl_matrix_complex* T_I,
        gsl_matrix_complex* G0, gsl_matrix_complex* V_NLO,
        gsl_matrix_complex* V_N2LO, gsl_matrix_complex* V_N3LO);



/*
 * ****************************************************************************
 * Helperfunctions to compute the Möller wave operators
 * ****************************************************************************
 */
gsl_matrix_complex* pw_moller_plus(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0);
gsl_matrix_complex* pw_moller_minus_dagger(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0);


gsl_matrix_complex* pw_T_mwpc_DWBA();

gsl_matrix_complex* full_BA_T_matrix();
gsl_matrix_complex* full_DWBA_T_matrix();

void make_tests(std::string chn_string);
void make_tests_DWBA(std::string chn_string);
void make_tests_DWBA_2(std::string chn_string);
void make_tests_DWBA_3(std::string chn_string);

void solve_DWB_from_potentials(Potential_mwpc<gsl_matrix>& pot1_real_noweights,
    Potential_mwpc<gsl_matrix_complex>& pot1_complex_weights,
    Potential_mwpc<gsl_matrix>& pot2_real_noweights,
    Potential_mwpc<gsl_matrix_complex>& pot2_complex_weights, 
    double number_of_p_points,
    double* p_grid, double* w_grid, bool FINITE_GRID, double Tl, 
    qs::quantum_channel chn, double REL_CORR);
}

/*
 * Method to compute sum_{i=0}^pow M^i with the on_shell multiplication
 */
void pow_matrix_on_shell_mult(gsl_matrix_complex* M,int pow,gsl_matrix_complex* res);

/*
 * This function performs returns \Omega^\dagger_- * M * \Omega^\dagger_+
 * using the on-shell multiplication
 */
void F(gsl_matrix_complex* omega_p,gsl_matrix_complex* omega_m_dagger,
        gsl_matrix_complex* M,gsl_matrix_complex* res);

#endif
