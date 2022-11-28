/*
    born_approx.h
    This file is part of the nn_mwpc project.
    
    This header contains functions to compute amplitude od scattering 
    processes in Born- and distorted-wave Born approximation.

    Oliver Thim 2022-09 --
    Department of Physics, Chalmers
*/
#ifndef BORN_APPROX
#define BORN_APPROX

#include "quantum_states.h"
#include "physics_helpers.h"
#include "gsl_blas.h"
#include "Constants.h"
#include "wigxjpf.h"
#include <vector>
#include <complex>
#include <iostream>
#include "gsl_matrix.h"
#include "gsl_blas.h"
#include "gsl_linalg.h"
#include "gsl_complex.h"

namespace dwba 
{
 
/*
 * This function compute the Born approximation to the scattering amplitude.
 * 
 * start_order: is the start position in the series, this is inculded
 * stop_order : is the highest order that is included, this is included
 * e.g. start_order=1,stop_order=2 returns (V*G0*V + V*G0*V*G0*V)
 * V    : is the potential matrix
 * G0   : is the free greens function
 *
 * V and G0 has to be in the form that the momentum gid weights and momenta
 * are included in the definition of the potential and propagator to make 
 * the resoultion of identity integrals to be just matrix multiplications
 * with the <star> matrix operation that I have defined.
 *
 * order=0: T = V
 * order=1: T = V + V*G0*V
 * order=2: T = V + V*G0*V + V*G0*V*G0*V
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
 * V_I, V_II and G0 has to be in the same form as accepted into 
 * 'dwba::pw_compute_BA()'.
 *
 * Define A = T_I*V^{-1}_I, B = V_II*V^{-1}_I*T_I.
 *
 * order=0: returns T_I
 * order=1: returns T_I + A*B
 * order=2: returns T_I + A*B + A*B*G0*B
 */

gsl_matrix_complex* pw_T_DWBA(int order, gsl_matrix_complex* T_I, 
        gsl_matrix_complex* V_I, gsl_matrix_complex* V_II, 
        gsl_matrix_complex* G0);


/*
 * Helperfunctions to compute the Möller wave operators
 */
gsl_matrix_complex* pw_moller_plus(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* V_I);
gsl_matrix_complex* pw_moller_minus_dagger(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* V_I);


gsl_matrix_complex* pw_T_mwpc_DWBA();

gsl_matrix_complex* full_BA_T_matrix();
gsl_matrix_complex* full_DWBA_T_matrix();

}

void pow_matrix_on_shell_mult(gsl_matrix_complex* M,int pow,gsl_matrix_complex* res);


#endif
