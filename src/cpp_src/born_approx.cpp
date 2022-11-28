/*
 * born_approx.cpp
 * This file is part of the nn_mwpc project.
 *
 * This file contatins the implementation of the function in the header
 * born_approx.h.
 *
 * Oliver Thim 2022-11 --
 * Department of Physics, Chalmers
 */
#include "born_approx.h"



gsl_matrix_complex* dwba::pw_T_BA(int start_order,int stop_order, gsl_matrix_complex* V, 
        gsl_matrix_complex* G0)
{
    gsl_matrix_complex* G0V = gsl_matrix_complex_alloc(V->size1,V->size2);

    ph::on_shell_mult(G0,V,G0V);
    
    gsl_matrix_complex* tmp1 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* tmp2 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* res = gsl_matrix_complex_alloc(V->size1,V->size2);

    gsl_matrix_complex_set_zero(res);

    for (int i=start_order; i<stop_order+1;i++)
    {
        // Compute (G0*V)^i
        pow_matrix_on_shell_mult(G0V,i,tmp1);
        ph::on_shell_mult(V,tmp1,tmp2);

        // Add to the result
        gsl_matrix_complex_add(res,tmp2);
            
    }
    // Deallocate 
    gsl_matrix_complex_free(G0V);
    gsl_matrix_complex_free(tmp1);
    gsl_matrix_complex_free(tmp2);

    return res;
}

gsl_matrix_complex* dwba::pw_T_DWBA(int order, 
        gsl_matrix_complex* T_I, gsl_matrix_complex* V_I,
        gsl_matrix_complex* V_II, gsl_matrix_complex* G0)
{
    // If the order is trivial
    if (order==0) 
    {
        return T_I;
    }
    // Get the Möller wave operators
    gsl_matrix_complex* omega_p = pw_moller_plus(T_I, V_I);
    gsl_matrix_complex* omega_m_dagger = pw_moller_minus_dagger(T_I, V_I);
    
    
    // Allocate the modified matrix
    gsl_matrix_complex* omega_p_G0 = gsl_matrix_complex_alloc(T_I->size1,
            T_I->size2);
        
    ph::on_shell_mult(omega_p,G0,omega_p_G0);

    // Perform the sum with the help of the Born approx code
    int start = 1;
    
    gsl_matrix_complex* VGV = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    pw_T_BA(start,order,V_II,omega_p_G0);
    
    // Multiply the VGVGV...V sum with the Möller operators from left and right
    // This is what I call F(...) in the notes
    
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    F(omega_p,omega_m_dagger,VGV,tmp);

    // Add the result to the leading order T_I matrix
    // T = T_I + tmp
    gsl_matrix_complex_add(tmp,T_I); // tmp <- tmp + T_I

    // Remove all temporary matrices
    gsl_matrix_complex_free(omega_p_G0);
    gsl_matrix_complex_free(omega_p);
    gsl_matrix_complex_free(omega_m_dagger);
    gsl_matrix_complex_free(VGV);

    return tmp;
}


/*
 *
 * Methods outside the namespace that is conisidered helper methods
 * that should not generally be acessed from outside of this file.
 *
 */

void pow_matrix_on_shell_mult(gsl_matrix_complex* M,int pow,gsl_matrix_complex* res)
{
    if (pow==0)
    {
        gsl_matrix_complex_set_identity(res);
        return;
    } else 
    {
        gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);
        gsl_matrix_complex_set_identity(tmp);
        
        gsl_matrix_complex_memcpy(res,M);
        for (int i = 0; i < pow; i++)
        {
            ph::on_shell_mult(res,M,tmp);
            gsl_matrix_complex_memcpy(res,tmp);

        }
        gsl_matrix_complex_free(tmp);
        return;
    }
}


void F(gsl_matrix_complex* omega_p,gsl_matrix_complex* omega_m_dagger,
        gsl_matrix_complex* M,gsl_matrix_complex* res)
{
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);

    // left
    ph::on_shell_mult(omega_m_dagger,M,tmp);
    // right
    ph::on_shell_mult(tmp,omega_p,res);

    gsl_matrix_complex_free(tmp);   
}
