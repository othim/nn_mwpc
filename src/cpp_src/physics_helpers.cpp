#include "physics_helpers.h"


gsl_matrix_complex* kronecker_product(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
{

    int new_size1 = m1->size1*m2->size1;
    int new_size2 = m1->size2*m2->size2;
    
    gsl_matrix_complex* prod = gsl_matrix_complex_alloc(new_size1, new_size2);
    
    // Create sub-matrices.

    for (int i = 0; i < new_size1; i++)
    {
        for (int j = 0; j < new_size2; j++)
        {
            /* First select the correct m1 matrix-element
             */
            int m1_r = (int)i / (int)m2->size1;
            int m1_c = (int)j / (int)m2->size2;
            gsl_complex m1_el = gsl_matrix_complex_get(m1,m1_r,m1_c);
            
            /* Get the element from the second matrix
             */
            int m2_r = (int)i % (int)m2->size1;
            int m2_c = (int)j % (int)m2->size2;
            gsl_complex m2_el = gsl_matrix_complex_get(m2,m2_r,m2_c);
            
            gsl_matrix_complex_set(prod,i,j,gsl_complex_mul(m1_el,m2_el));

        }
    }
    return prod;
}


gsl_matrix_complex* get_Pauli_matrix(char axis)
{
    gsl_matrix_complex* sigma = gsl_matrix_complex_alloc(2,2);

    if (axis=='x')
    {
        gsl_matrix_complex_set(sigma,0,0,gsl_complex_rect(0,0));
        gsl_matrix_complex_set(sigma,0,1,gsl_complex_rect(1,0));
        gsl_matrix_complex_set(sigma,1,0,gsl_complex_rect(1,0));
        gsl_matrix_complex_set(sigma,1,1,gsl_complex_rect(0,0));
    } else if (axis=='y')
    {
        gsl_matrix_complex_set(sigma,0,0,gsl_complex_rect(0,0));
        gsl_matrix_complex_set(sigma,0,1,gsl_complex_rect(0,-1));
        gsl_matrix_complex_set(sigma,1,0,gsl_complex_rect(0,1));
        gsl_matrix_complex_set(sigma,1,1,gsl_complex_rect(0,0));
  
    } else if (axis=='z')
    {
        gsl_matrix_complex_set(sigma,0,0,gsl_complex_rect(1,0));
        gsl_matrix_complex_set(sigma,0,1,gsl_complex_rect(0,0));
        gsl_matrix_complex_set(sigma,1,0,gsl_complex_rect(0,0));
        gsl_matrix_complex_set(sigma,1,1,gsl_complex_rect(-1,0));
  
    }
    return sigma;
}

gsl_complex trace(gsl_matrix_complex* m)
{
    if (m->size1 != m->size2)
    {
        std::cout << "invlid matrix" << std::endl;
        return gsl_complex_rect(-1,0);
    }
    gsl_complex trace = gsl_complex_rect(0,0);
    for (int i = 0; i < m->size1; i++)
    {
        // Add up all the diagonal elements
        gsl_complex_add(trace,gsl_matrix_complex_get(m,i,i));
    }
    return trace;
}

/* Note here that ALL quantum numbers are doubled so everythin
 * becomes integers.
 */

double CG_coeff(int J_2, int M_2, int j1_2, int j2_2, int m1_2, int m2_2)
{
    wig_table_init(2*100, 9);
    wig_temp_init(2*100);

    double wig = wig3jj(j1_2 , j2_2 , J_2 ,
                        m1_2 , m2_2 , -M_2 );

    wig_temp_free();
    wig_table_free();
    return std::pow(-1.0,j1_2/2.0+j2_2/2.0+M_2/2.0)*std::sqrt(J_2 + 1)*wig;
} 
