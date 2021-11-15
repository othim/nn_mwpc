/*
 * pot_ext.h
 * This file is part of the nn_mwpc project.
 *
 * This class is a potential class that computes potential
 * matrices from a given function for computing matrix
 * elements that is given.
 *
 * Oliver Thim 2021-11 --
 * Department of Physics, Chalmers
 */

#ifndef POT_EXT
#define POT_EXT

#include "quantum_states.h"
#include "Constants.h"
#include "gsl_matrix.h"
#include <iostream>
#include <cmath>

class Potential_ext
{
private:
    /*
     * Function pointer to the function that computes the matrix elements of
     * the potential.
     */ 
    void (*my_element_V_arr)(double qi,double qo, bool coupled, int J, double* V_arr);
    
    double* p_grid_;
    int p_grid_length_;    
    double cutoff_Lambda_;

pubilc:

    Potential_ext(double* p_grid, int p_grid_length, void (*f)(double qi,double qo, bool coupled, int J, double* V_arr));
    ~Potential_ext();
    /*
     * This function returns the potential matrix
     */
    gsl_matrix* get_matrix(double q_on_shell, qs::quantum_channel chn);
    
    /*
     * This function returns the potneital matrix without an 
     * on-shell point.
     */
    gsl_matrix*  get_matrix_no_on_shell(qs::quantum_channel chn);
}
#endif
