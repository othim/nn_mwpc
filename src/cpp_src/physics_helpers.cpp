#include "physics_helpers.h"

void ph::physics_helpers_init()
{
    // Init wiegner symbol tables
    wig_table_init(2*100, 9);
    wig_temp_init(2*100);
}

void ph::physics_helpers_free()
{
    // Deallocate tables
    wig_temp_free();
    wig_table_free();
}



gsl_matrix_complex* ph::kronecker_product(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
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


gsl_matrix_complex* ph::get_Pauli_matrix(char axis)
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

gsl_complex ph::trace(gsl_matrix_complex* m)
{
    if (m->size1 != m->size2)
    {
        std::cout << "invlid matrix" << std::endl;
        return gsl_complex_rect(-1,0);
    }
    gsl_complex trace = gsl_complex_rect(0,0);
    for (int i = 0; i < (int)m->size1; i++)
    {
        // Add up all the diagonal elements
        gsl_complex_add(trace,gsl_matrix_complex_get(m,i,i));
    }
    return trace;
}

/* Note here that ALL quantum numbers are doubled so everything
 * becomes integers. Remember that you have to run the init function
 * before you use this function and then run the free function when
 * you are done.
 */

double ph::CG_coeff(int J_2, int M_2, int j1_2, int j2_2, int m1_2, int m2_2)
{

    double wig = wig3jj(j1_2 , j2_2 , J_2 ,
                        m1_2 , m2_2 , -M_2 );
    std::cout << "wig" << wig << std::endl;
    std::cout << m1_2+ m2_2 + M_2 << std::endl;
    return std::pow(-1.0,j1_2/2.0+j2_2/2.0+M_2/2.0)*std::sqrt(J_2 + 1)*wig;
} 

void ph::gauss_legendre_finite_mesh(unsigned int Number_of_points, double min, 
        double max, double**p, double** w)
{
    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    // Make GL-grid from -1 to 1
    gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, 
            Number_of_points, min, max, 0, 0);
    
    double* p_grid = gsl_integration_fixed_nodes(int_ang_);
    double* w_grid = gsl_integration_fixed_weights(int_ang_);
    
    double* pp = (double*)malloc(Number_of_points*sizeof(double));
    double* ww = (double*)malloc(Number_of_points*sizeof(double));
  
    for (int i = 0; i < (int)Number_of_points; i++)
    {
        pp[i] = p_grid[i];
        ww[i] = w_grid[i];
        //std::cout << pp[i] << "   " << ww[i] << std::endl;
    }
    *p = pp;
    *w = ww;
}

void ph::gauss_legendre_inf_mesh(unsigned int Number_of_points, double scale,double** p,double** w)
{
    // Make grid from -1 to 1
    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, Number_of_points, -1.0, 1.0, 0, 0);
    
    //gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, 
    //        Number_of_points, 0, 4300, 0, 0);
    double* p_grid = gsl_integration_fixed_nodes(int_ang_);
    double* w_grid = gsl_integration_fixed_weights(int_ang_);
    // Make transformation
    double pi_4 = M_PI/4.0;

    double* pp = (double*)malloc(Number_of_points*sizeof(double));
    double* ww = (double*)malloc(Number_of_points*sizeof(double));
  
    for (int i = 0; i < (int)Number_of_points; i++)
    {
        double x = p_grid[i];
        pp[i] = scale*tan(pi_4*(x+1));
        ww[i] = (scale*pi_4/(cos(pi_4*(x+1))*cos(pi_4*(x+1))))*w_grid[i];
        //pp[i] = p_grid[i];
        //ww[i] = w_grid[i];
        //std::cout << pp[i] << "   " << ww[i] << std::endl;

    }
    *p = pp;
    *w = ww;
}


ph::eigen_t ph::solve_SE(double* p, double* w, unsigned int number_of_grid_points
        ,qs::quantum_channel chn, const gsl_matrix* V)
{
   // The potential is assumed to be in a partial wave basis with normalization 
   // <p'|p> = \delta(p'-p)/p^2, with \pi/2 factor from Landau removed.
   if (chn.coupled)
   {
      if (V->size1 != 2*number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   } else 
   {
     if (V->size1 != number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   }
   // Get reduced mass of system
   double mu = ph::get_mN(chn.Tz)/2.0;
   
   gsl_matrix* H = gsl_matrix_alloc(V->size1,V->size2);
   
   // Construct Hamiltonian
   for (int i = 0; i < (int)H->size1; i++)
   {
      for (int j=0; j < (int)H->size2; j++)
      {
         // This is to still use the same momenta
         int l = j;
         int k = i;
         if (!(j<(int)number_of_grid_points))
         {
            l = j-number_of_grid_points;
         }
         if (!(i<(int)number_of_grid_points))
         {
            k = i-number_of_grid_points;
         }
         double p2 = p[l]*p[l];
         
         // These two gives the same answer...
         double el = gsl_matrix_get(V,i,j)*p[l]*p[k]*sqrt(w[l]*w[k]);
         //double el = gsl_matrix_get(V,i,j)*p2*w[l];
         if (i==j) {
            el += p2/(2.0*mu);
         }
         gsl_matrix_set(H,i,j,el);
      }
   }
   //std::cout << "Hamiltonian matrix" << std::endl;
   //print_m(H);

   // Diagonalize the matrix
   gsl_vector_complex* eval = gsl_vector_complex_alloc(V->size1);
   gsl_matrix_complex* evec = gsl_matrix_complex_alloc(V->size1, V->size1);

   gsl_eigen_nonsymmv_workspace* ws = gsl_eigen_nonsymmv_alloc(V->size1);
   gsl_eigen_nonsymmv(H,eval,evec,ws);
   
   gsl_eigen_nonsymmv_free(ws);

   gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
   ph::eigen_t e;

   e.eigenvalues  = eval;
   e.eigenvectors = evec;
   
   // Use a test matrix
   /* finds eigenvalues to this...
   gsl_matrix* H_test = gsl_matrix_alloc(2,2);
   gsl_matrix_set(H_test,0,0,1);
   gsl_matrix_set(H_test,0,1,0);
   gsl_matrix_set(H_test,1,0,0);
   gsl_matrix_set(H_test,1,1,1);

   eval = gsl_vector_complex_alloc(2);
   gsl_eigen_nonsymm_workspace* ws_test = gsl_eigen_nonsymm_alloc(2);
   gsl_eigen_nonsymm(H_test,eval,ws_test);
   */
   
   return e;
}

ph::eigen_t_herm ph::solve_SE_complex_weights(double* p, double* w, 
        unsigned int number_of_grid_points
        ,qs::quantum_channel chn, const gsl_matrix_complex* V)
{
   // The potential is assumed to be in a partial wave basis with normalization 
   // <p'|p> = \delta(p'-p)/p^2, with \pi/2 factor from Landau removed.
   if (chn.coupled)
   {
      if (V->size1 != 2*number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   } else 
   {
     if (V->size1 != number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   }
   // Get reduced mass of system
   double mu = ph::get_mN(chn.Tz)/2.0;
   
   gsl_matrix_complex* H = gsl_matrix_complex_alloc(V->size1,V->size2);
   
   // Construct Hamiltonian
   for (int i = 0; i < (int)H->size1; i++)
   {
      for (int j=0; j < (int)H->size2; j++)
      {
         // This is to still use the same momenta
         int l = j;
         if (!(j<(int)number_of_grid_points))
         {
            l = j-number_of_grid_points;
         }
         double p2 = p[l]*p[l];
         
         // These two gives the same answer...
         gsl_complex el = gsl_matrix_complex_get(V,i,j);
         
         if (i==j) {
            el = gsl_complex_add(gsl_complex_rect(p2/(2.0*mu),0.0),el);
         }
         gsl_matrix_complex_set(H,i,j,el);
      }
   }
   //std::cout << "Hamiltonian matrix" << std::endl;
   //print_m(H);

   // Diagonalize the matrix
   gsl_vector* eval = gsl_vector_alloc(V->size1);
   gsl_matrix_complex* evec = gsl_matrix_complex_alloc(V->size1, V->size1);

   gsl_eigen_hermv_workspace* ws = gsl_eigen_hermv_alloc(V->size1);
   gsl_eigen_hermv(H,eval,evec,ws);
   
   gsl_eigen_hermv_free(ws);

   gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
   ph::eigen_t_herm e;

   e.eigenvalues  = eval;
   e.eigenvectors = evec;
   
   // Use a test matrix
   /* finds eigenvalues to this...
   gsl_matrix* H_test = gsl_matrix_alloc(2,2);
   gsl_matrix_set(H_test,0,0,1);
   gsl_matrix_set(H_test,0,1,0);
   gsl_matrix_set(H_test,1,0,0);
   gsl_matrix_set(H_test,1,1,1);

   eval = gsl_vector_complex_alloc(2);
   gsl_eigen_nonsymm_workspace* ws_test = gsl_eigen_nonsymm_alloc(2);
   gsl_eigen_nonsymm(H_test,eval,ws_test);
   */
   
   return e;
}


double ph::rad_to_deg(double in)
{
   return in*180.0/M_PI;
}


void ph::print_m(gsl_matrix* matrix)
{
   std::cout << "---------" << std::endl;
   for (int i = 0; i < (int)matrix->size1; i++)
   {
      for (int j = 0; j < (int)matrix->size1; j++)
      {
         std::cout << gsl_matrix_get(matrix,i,j) << " ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}

void print_v(gsl_vector* vec)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < vec->size; i++)
   {
        std::cout << gsl_vector_get(vec,i) << "\n";
   }
   std::cout << "---------" << std::endl;
}

void ph::print_m(gsl_matrix_complex* matrix)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < matrix->size1; i++)
   {
      for (std::size_t j = 0; j < matrix->size1; j++)
      {
          gsl_complex el = gsl_matrix_complex_get(matrix,i,j); 
         std::cout << "(" <<  GSL_REAL(el) << "," << GSL_IMAG(el) << ") ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}

void ph::print_v(gsl_vector_complex* vec)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < vec->size; i++)
   {
        std::cout << "(" << GSL_REAL(gsl_vector_complex_get(vec,i)) << "," << 
        GSL_IMAG(gsl_vector_complex_get(vec,i)) << ") \n";
   }
   std::cout << "---------" << std::endl;
}


void ph::print_m_complex_to_file(std::string path,gsl_matrix_complex* matrix)
{
    FILE* f = std::fopen(path.c_str(), "wb");
    gsl_matrix_complex_fprintf(f, matrix,"%e");
    std::fclose (f);
}


void ph::on_shell_mult(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res)
{
    // Multiply as usual
    gsl_complex alpha = gsl_complex_rect(1,0);
    gsl_complex beta = gsl_complex_rect(0,0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, m1,m2,beta,res);

    // Correct for the inclusion of the on-shell elements in the 
    // matrix-matrix product
    
    // Take the last column of m1
    gsl_vector_complex_view m1_vv = gsl_matrix_complex_column(m1,m1->size2-1);
    gsl_vector_complex* v1 = &m1_vv.vector;
    // Take the last row of m2
    gsl_vector_complex_view m2_vv = gsl_matrix_complex_row(m2,m2->size1-1);
    gsl_vector_complex* v2 = &m2_vv.vector;
    
    // Perform the outer product
    gsl_complex minus_one = gsl_complex_rect(-1,0);
    gsl_blas_zgeru(minus_one, v1,v2, res);
}

void ph::mult(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res)
{
    // Multiply as usual
    gsl_complex alpha = gsl_complex_rect(1,0);
    gsl_complex beta = gsl_complex_rect(0,0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, m1,m2,beta,res);
}

void ph::mult3(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* m3, gsl_matrix_complex* res)
{   
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(m1->size1,
            m1->size2);
    
    ph::mult(m2,m3,tmp);
    ph::mult(m1,tmp,res);

    gsl_matrix_complex_free(tmp);
}

void ph::mult5(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* m3, gsl_matrix_complex* m4, gsl_matrix_complex* m5,
        gsl_matrix_complex* res)
{
    gsl_matrix_complex* tmp1 = gsl_matrix_complex_alloc(m1->size1,
            m1->size2);
    gsl_matrix_complex* tmp2 = gsl_matrix_complex_alloc(m1->size1,
            m1->size2);
    
    ph::mult(m4,m5,tmp1);
    ph::mult(m3,tmp1,tmp2);
    ph::mult(m2,tmp2,tmp1);
    ph::mult(m1,tmp1,res);

    gsl_matrix_complex_free(tmp1);
    gsl_matrix_complex_free(tmp2);
}

void ph::mult7(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* m3, gsl_matrix_complex* m4, gsl_matrix_complex* m5,
        gsl_matrix_complex* m6, gsl_matrix_complex* m7,gsl_matrix_complex* res)
{
    gsl_matrix_complex* tmp1 = gsl_matrix_complex_alloc(m1->size1,
            m1->size2);
    gsl_matrix_complex* tmp2 = gsl_matrix_complex_alloc(m1->size1,
            m1->size2);
    
    ph::mult(m6,m7,tmp1);
    ph::mult(m5,tmp1,tmp2);
    ph::mult(m4,tmp2,tmp1);
    ph::mult(m3,tmp1,tmp2);
    ph::mult(m2,tmp2,tmp1);
    ph::mult(m1,tmp1,res);

    gsl_matrix_complex_free(tmp1);
    gsl_matrix_complex_free(tmp2);
}

void ph::on_shell_mult_bf(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res)
{
    for(int i=0; i<(int)m1->size1; i++)
    {
        for(int j=0; j<(int)m2->size2; j++)
        {
            gsl_complex el = gsl_complex_rect(0,0);
            for(int k=0; k<(int)m1->size1-1; k++)
            {
                el = gsl_complex_add(el, gsl_complex_mul(gsl_matrix_complex_get(m1,i,k),
                            gsl_matrix_complex_get(m2,k,j)));
            }
            gsl_matrix_complex_set(res,i,j,el);
        }
    }
}



void ph::make_matrix_complex(gsl_matrix_complex* Mz,gsl_matrix* M)
{
    for (int i = 0; i < (int)M->size1; i++)
    {
        for (int j = 0; j < (int)M->size2; j++)
        {
            double el = gsl_matrix_get(M,i,j);
            gsl_complex elz = gsl_complex_rect(el,0.0);
            gsl_matrix_complex_set(Mz,i,j,elz);
        }
    }
}

void ph::matrix_set(gsl_matrix* m, const size_t i, const size_t j, double el)
{
    gsl_matrix_set(m,i,j,el);
}

void ph::matrix_set(gsl_matrix_complex* m, const size_t i, const size_t j, 
        const gsl_complex el)
{
    gsl_matrix_complex_set(m,i,j,el);

}

void ph::matrix_set(gsl_matrix_complex* m, const size_t i, const size_t j, 
        const double el)
{
    gsl_matrix_complex_set(m,i,j,gsl_complex_rect(el,0.0));

}

gsl_matrix* ph::matrix_alloc(const size_t i,const size_t j,gsl_matrix* m)
{
    return gsl_matrix_alloc(i,j);
}

gsl_matrix_complex* ph::matrix_alloc(const size_t i,const size_t j,gsl_matrix_complex* m)
{
    return gsl_matrix_complex_alloc(i,j);
}

void ph::matrix_free(gsl_matrix* m)
{
    gsl_matrix_free(m);
}

void ph::matrix_free(gsl_matrix_complex* m)
{
    gsl_matrix_complex_free(m);
}

void ph::matrix_set_zero(gsl_matrix* m)
{
    gsl_matrix_set_zero(m);
}

void ph::matrix_set_zero(gsl_matrix_complex* m)
{
    gsl_matrix_complex_set_zero(m);
}

void ph::matrix_add(gsl_matrix* m1, gsl_matrix* m2)
{
    gsl_matrix_add(m1,m2);
}

void ph::matrix_add(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
{
    gsl_matrix_complex_add(m1,m2);
}

void ph::matrix_sub(gsl_matrix* m1, gsl_matrix* m2)
{
    gsl_matrix_sub(m1,m2);
}

void ph::matrix_sub(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
{
    gsl_matrix_complex_sub(m1,m2);
}

void ph::matrix_memcpy(gsl_matrix* m1, gsl_matrix* m2)
{
    gsl_matrix_memcpy(m1,m2);
}

void ph::matrix_memcpy(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
{
    gsl_matrix_complex_memcpy(m1,m2);
}

void ph::matrix_scale(gsl_matrix* m1, const double scale)
{
    gsl_matrix_scale(m1, scale);
}

void ph::matrix_scale(gsl_matrix_complex* m1, const double scale)
{
    gsl_complex alpha = gsl_complex_rect(0.0,0.0);
    gsl_complex beta  = gsl_complex_rect(scale,0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, m1, m1, beta, m1); 
}

void ph::matrix_from_vector(gsl_matrix_complex* M,gsl_vector_complex* vec)
{
    gsl_matrix_complex_set_zero(M);
    // Set the diagonal values to the vecotr values
    for (int i = 0; i < (int)M->size1; i++)
    {
        gsl_matrix_complex_set(M,i,i,gsl_vector_complex_get(vec,i));
    }
}

double ph::get_mN(int Tz)
{
    if (Tz == -1)
    {
        return constants::Mn; // nn
    } else if (Tz == 0)
    {
        return 2.0*constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np

    } else if (Tz == 1)
    {
        return constants::Mp; // pp
    } else 
    {
        std::cout << "Error in get_mN(): Unknown isospin" << std::endl;
        return -1;
    }
}
