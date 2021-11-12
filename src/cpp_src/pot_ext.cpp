#include "pot_ext.h"



Potential_ext::Potential_ext(double* p_grid, int p_grid_length, void (*f)(double qi,double qo, bool coupled, int J, double* V_arr))
{
    // Set the function opinter to the correct function
    my_element_V_arr = f; 

    p_grid_ = p_grid;
    p_grid_length_ = p_grid_length;
}


gsl_matrix* Potential_ext::get_matrix(double q_on_shell, qs::quantum_channel chn)
{}
    
gsl_matrix*  Potential_ext::get_matrix_no_on_shell(qs::quantum_channel chn)
{}
