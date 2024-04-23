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
#include "gsl_matrix.h"
#include <iostream>
#include <cmath>
#include "gsl_pow_int.h"
#include "potential.h"

template <class gsl_m>
class Potential_ext : public potential<gsl_m>
{
private:
    /*
     * Function pointer to the function that computes the matrix elements of
     * the potential.
     */ 
    void (*my_element_V_arr)(double qi,double qo, bool coupled, int S, int J, int T, int Tz, double* V_arr);
    
    double* p_grid_;
    int mom_grid_size_;    
    double cutoff_Lambda_;

public:

    Potential_ext(double* p_grid, int p_grid_length, double cutoff_Lambda, 
            void (*f)(double qi,double qo, bool coupled, int S, int J, int T, 
                int Tz,  double* V_arr));
    
    ~Potential_ext();
    
    /*
     * This function returns the potential matrix
     */
    gsl_matrix* get_matrix(double q_on_shell, qs::quantum_channel chn,
            bool rel_correction);
    
    /*
     * This function returns the potneital matrix without an 
     * on-shell point.
     */
    gsl_matrix*  get_matrix_no_onshell(qs::quantum_channel chn, 
            bool rel_correction);
    
    /*
     * These functions are trivial in this class.
     */
    void populate_saved_mtx(qs::quantum_channel chn, bool rel_correction);
    gsl_m* get_saved_matrix(double q_on_shell, qs::quantum_channel chn,
            bool rel_correction);
    
    void print_LECs_and_params_info();
};


void cdbonn_correct_arg(double qi, double qo, bool coupled, int S, int J, 
        int T, int Tz,  double* V_arr);
void nijm_correct_arg(double qi, double qo, bool coupled, int S, int J, int T, 
        int Tz,  double* V_arr);
void nijm_OPE_correct_arg(double qi, double qo, bool coupled, int S, int J, 
        int T, int Tz,  double* V_arr);

extern "C" {
    void cdbonn_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

extern "C" {
    void nijmegen_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}
#endif
