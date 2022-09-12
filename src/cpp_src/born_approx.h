/*
    born_approx.h
    This file is part of the nn_mwpc project.
    
    This header contains functions to compute amplitude od scattering 
    processes in Born- and distorted-wave Born approximation.

    Oliver Thim 2022-09 --
    Department of Physics, Chalmers
*/
#include "quantum_states.h"
#include "physics_helpers.h"
#include "gsl_blas.h"
#include "Constants.h"
#include "wigxjpf.h"
#include <vector>
#include <complex>
#include <iostream>
#include "gsl_matrix.h"

namespace dwba 
{
 
/*
 * This function compute the Born approximation to the scattering amplitude.
 * 
 * order: is the order of the approximation.
 * V    : is the potential matrix
 * G0   : is the free greens function
 *
 * order=0: returns T = V
 * order=1: returns T = V + V*G0*V
 * order=2: returns T = V + V*G0*V + V*G0*V*G0*V
 * ...
 */

gsl_matrix* pw_compute_BA(int order, gsl_matrix* V, gsl_matrix* G0);



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
 * Define A = T_I*V^{-1}_I, B = V_II*V^{-1}_I*T_I.
 *
 * order=0: returns T_I
 * order=1: returns T_I + A*B
 * order=2: returns T_I + A*B + A*B*G0*B
 */

gsl_matrix_complex* pw_compute_DWBA(int order, gsl_matrix_complex*
        T_I, gsl_matrix_complex* V_I, gsl_matrix_complex* V_II, 
        gsl_matrix_complex* G0);



gsl_matrix_complex* full_BA_T_matrix();
gsl_matrix_complex* full_DWBA_T_matrix();

}
