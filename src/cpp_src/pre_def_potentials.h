/* pre_def_potentials.h
 * This file is a part of the nn_mwpc project.
 *
 * This file contains potential definitions made from the Potential_mwpc class
 *
 * Oliver Thim 2023-04
 * Department of Physics, Chalmers
 */
#pragma once

#include "potential_mwpc.h"
#include "physics_helpers.h"

namespace pre_def_pot
{

    /*
     * Function that defines potentials and returns a 
     * Potential_mwpc<gls_matrix> object
     *
     */
    Potential_mwpc<gsl_matrix_complex>* create_pre_def_pot(
            std::string pre_def_name, int ang_int_points_, double* p_grid_, 
            double* w_grid_, int number_of_p_points_, int J_max_in_pot_,
            double cutoff_, int cut_pow_, bool sharp_cutoff_, 
            double sharp_cutoff_add_, bool cut_on_shell_,
            ph::constants_struct* program_const_);

}
