/* pre_def_potentials.h
 * This file is a part of the nn_mwpc project.
 *
 * This file contains potential definitions made from the potential_mwpc class
 *
 * Oliver Thim 2023-04
 * Department of Physics, Chalmers
 */
#pragma once

#include "potential_mwpc.h"



namespace pre_def_pot
{

    /*
     * Function that defines potentials and populate potentials_ and
     * potential_names_
     *
     */
    Pot_mwpc<gsl_matrix_complex>* load_pre_def_pot(
            std::string pre_def_name);




}
