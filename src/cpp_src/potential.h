/*
 * potential.h
 * This file is part of the nn_mwpc project.
 *
 * This class is a virtual base class that defines the necessary constituents
 * of the potential. Two important classes will be derived from this class,
 * pot_mwpc and pot_ext.
 *
 * Oliver Thim 2023
 * Department of Physics, Chalmers
 */
#pragma once

#include "quantum_states.h"
#include "Constants.h"
#include "gsl_matrix.h"
#include <iostream>
#include <cmath>
#include "gsl_pow_int.h"




template <class gsl_m>
class potential
{
private:


public:

    potential();
    ~potential();

    virtual gsl_m* get_matrix(double q_on_shell, qs::quantum_channel chn);
    
    virtual gsl_m* get_matrix_on_onshell(qs::quantum_channel chn);

    virtual gsl_m* get_saved_matrix(double q_on_shell, qs::quantum_channel chn);

    virtual void populate_save_mtx();

    virtual void print_LECs_and_params_info();
};
