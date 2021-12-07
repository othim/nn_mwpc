#include "LS_Solver.h"

//#define ENABLE_DEBUG

// Constructor
LS_Solver::LS_Solver(std::vector<qs::quantum_channel> channels, unsigned int mom_grid_size,
    double mom_grid_scale, bool cutoff_enabled, double cutoff_Lambda, bool relcorr_enabled)
{
    #ifdef ENABLE_DEBUG
        std::cout << "LS_Solver()" << std::endl;
    #endif
    // Init variables
    mom_grid_size_ = mom_grid_size;
    cutoff_enabled_ = cutoff_enabled;
    cutoff_Lambda_ = cutoff_Lambda;
    relcorr_enabled_ = relcorr_enabled;
    
    // Make GL-grid
    ph::gauss_legendre_inf_mesh(mom_grid_size_,mom_grid_scale,&p_grid_,&w_grid_);
}

// Destructor
LS_Solver::~LS_Solver()
{
    delete p_grid_;
    delete w_grid_;
}

gsl_vector* LS_Solver::setup_D_vector(double q_on_shell, bool coupled, double mu)
{
    // See eq (18.19) in Landau
    #ifdef ENABLE_DEBUG
        std::cerr << "setup_D_vector()" << std::endl;
    #endif

    // This is the factor in front of the LS equation.
    // In some cases it is 2/pi when that is a factor in
    // the partial wave momentum base normalization.
    // For now, the code is written without this factor, it will affect the 
    // rho value when raltiong the T/R matrix to the S-matrix/phase-shifts.
    double fac = 1.0;
    // D - energy denominator and weights
    gsl_vector* D_vector;
    double q2_on_shell = q_on_shell*q_on_shell;

    if (coupled)
    {
        // If the channel is coupled D is of length 
        // 2*mom_grid_size_ + 2
        D_vector = gsl_vector_alloc(2*mom_grid_size_ + 2);
    } else 
    {
        // Channel is uncoupled D is of length mom_grid_size_ + 1
        D_vector = gsl_vector_alloc(mom_grid_size_ + 1);
    }
   
    
    
    for (int i = 0; i < mom_grid_size_; i++)
    {
        double p2 = p_grid_[i]*p_grid_[i];
        double el = (2.0*mu)*(fac)*w_grid_[i]*p2/(p2-q2_on_shell); // NO CUTOFF

        if (coupled)
        {
            gsl_vector_set(D_vector,i,el);
            gsl_vector_set(D_vector,i+ mom_grid_size_+1,el);       
        } else
        {
            gsl_vector_set(D_vector,i,el);
        }
    }

    // Add the last element
    double sum = 0;
    for (int i=0; i < mom_grid_size_; i++)
    {
        sum += w_grid_[i]/(p_grid_[i]*p_grid_[i]-q2_on_shell);
    }

    double el = -(fac)*2.0*mu*q2_on_shell*sum; // NO CUTOFF

    if (coupled)
    {
        gsl_vector_set(D_vector,mom_grid_size_,el);
        gsl_vector_set(D_vector,2*mom_grid_size_+1,el);
    } else
    {
        gsl_vector_set(D_vector,mom_grid_size_,el);
    }
    return D_vector;
}

gsl_matrix * my_diag_alloc(const gsl_vector * X)
{
    gsl_matrix * mat = gsl_matrix_alloc(X->size, X->size);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_matrix_set_all(mat, 0.0); //or whatever number you like
    gsl_vector_memcpy(&diag.vector, X);
    return mat;
}

int my_scale_columns(gsl_matrix* a, const gsl_vector* b)
{
    gsl_matrix* m_diag = my_diag_alloc(b);
    
    gsl_matrix* r = gsl_matrix_alloc(a->size1, a->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, a, m_diag, 0.0, r);
    gsl_matrix_memcpy(a,r);
    gsl_matrix_free(m_diag);
    gsl_matrix_free(r);
    return 0;
}

gsl_matrix* LS_Solver::setup_F_matrix(bool coupled, gsl_vector* D_vector, gsl_matrix* V_mtx)
{
    #ifdef ENABLE_DEBUG
        std::cerr << "setup_F_matrix()" << std::endl;
    #endif
    // F_ij = \delta_ij + D_j V_ij (no sum over i)

    // Create matrix pointer and allocate the size for the 
    // coupled and uncoupled case
    gsl_matrix* F_mtx;
    gsl_matrix* V_copy;
    if (coupled)
    {
        F_mtx = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
        gsl_matrix_set_identity(F_mtx); // Set to unity
        V_copy = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
        gsl_matrix_memcpy(V_copy,V_mtx);
    } else 
    {
        F_mtx = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
        gsl_matrix_set_identity(F_mtx); // Set to unity
        V_copy = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
        gsl_matrix_memcpy(V_copy,V_mtx);
    }

    // Compute the matrix elements according to (18.22) in Landau
    
    /*
    // Test
    gsl_matrix* m = gsl_matrix_alloc(2,2);
    gsl_matrix_set(m,0,0,1);
    gsl_matrix_set(m,0,1,3);
    gsl_matrix_set(m,1,0,2);
    gsl_matrix_set(m,1,1,6);

    gsl_vector* v = gsl_vector_alloc(2);
    gsl_vector_set(v,0,3);
    gsl_vector_set(v,1,4);

    //gsl_matrix_scale_columns(m,v);
    my_scale_columns(m,v);
    ph::print_m(m);
    */


    // Scale V by D
    // OLD
    // ---
    //gsl_matrix_scale_columns(V_copy, D_vector);
    // ---
    my_scale_columns(V_copy, D_vector);

    // Add V D to F
    gsl_matrix_add(F_mtx,V_copy);
    gsl_matrix_free(V_copy);

    return F_mtx;
}

void print_matrix(gsl_matrix* matrix)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < matrix->size1; i++)
   {
      for (std::size_t j = 0; j < matrix->size1; j++)
      {
         std::cout << gsl_matrix_get(matrix,i,j) << " ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}

/*
    Converts the phase shifts in the BB convention to the Stapp
    convention
*/
Phase_shifts_chn BB_to_Stapp(Phase_shifts_chn ps)
{
    Phase_shifts_chn phases;
    phases.epsilon = 0.5*asin(sin(2*ps.epsilon)* sin(ps.delta_m - ps.delta_p));
    //std::cout << phases.epsilon << std::endl;
    phases.delta_m = 0.5*(ps.delta_p + ps.delta_m + asin(tan(2*phases.epsilon)/tan(2*ps.epsilon)));
    //std::cout << tan(2*phases.epsilon)/tan(2*ps.epsilon) << std::endl;
    phases.delta_p = 0.5*(ps.delta_p + ps.delta_m - asin(tan(2*phases.epsilon)/tan(2*ps.epsilon)));

    phases.delta_uncoupled = ps.delta_uncoupled;

    return phases;
}

void LS_Solver::get_mu_q_on_shell(double T_lab, qs::quantum_channel chn, double* mu, double* q_on_shell)
{
    if (chn.Tz == -1)
    {
        *mu = constants::Mn/2.0; // nn
        *q_on_shell = sqrt(*mu*T_lab);
    } else if (chn.Tz == 0)
    {
        *mu = constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np
        *q_on_shell = sqrt(constants::Mp*constants::Mp*T_lab*(T_lab + 2.0*constants::Mn)/
            ((constants::Mp + constants::Mn)*(constants::Mp + constants::Mn) + 2.0*T_lab*constants::Mp));
    } else if (chn.Tz == 1)
    {
        *mu = constants::Mp/2.0; // pp
        *q_on_shell = sqrt(*mu*T_lab);
    } else 
    {
        #ifdef ENABLE_DEBUG
            std::cerr << "Error in get_mu_q_on_shell(): Unknown isospin" << std::endl;
        #endif
    }
}

/*
    This function returns the phase shifts in at the desired lab energy T_lab for the channel chn.
    The phase shifts are returned in the Stapp convention
*/
Phase_shifts_chn LS_Solver::solve_in_chn_R(double T_lab, qs::quantum_channel chn, gsl_matrix* pot_V_mtx)
{
    #ifdef ENABLE_DEBUG
        std::cerr << "solve_in_chn()" << std::endl;
    #endif
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    //std::cout << "q_on_shell=" << q_on_shell << std::endl;
    // rho = 2*q*mu   
    // This is a convecntion dependent parameter that relates the 
    // R-matrix to the T/S matrices. rho will be different if a different 
    // normalization for the |klm> quantum states is choosen. To see the conventions
    // Used in this code see the README.md file.
    double rho = (M_PI/2.0)*2.0*q_on_shell*mu;

    //std::cout << "Potential" << std::endl;
    //print_matrix(pot_V_mtx);
   
    // Setup D-vector
    gsl_vector* D_vector = setup_D_vector(q_on_shell,chn.coupled,mu);
    //print_vector(D_vector);
    // Setup F-matrix
    gsl_matrix* F_matrix = setup_F_matrix(chn.coupled,D_vector,pot_V_mtx);
    
    // Solve the matrix equation F*R = V
        
    // LU decompose
    gsl_permutation* perm = gsl_permutation_alloc(F_matrix->size1);
    int signum;
    
    gsl_linalg_LU_decomp(F_matrix,perm,&signum);
   
    // Invert from LU decompusition
    gsl_matrix* inverse = gsl_matrix_alloc(F_matrix->size1,F_matrix->size2);
    gsl_linalg_LU_invert(F_matrix,perm,inverse);
   
    gsl_matrix* R_result = gsl_matrix_alloc(F_matrix->size1,F_matrix->size2);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse, pot_V_mtx, 0.0, R_result); 
   
    // Get the matrix elements corresponding to the on-shell R-matrix

    Phase_shifts_chn phase_shifts;
    if (chn.coupled) 
    {
        // on-shell R-matrix is 2x2
        double R_pp,R_mm,R_mp;
        R_mm = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        R_mp = gsl_matrix_get(R_result,2*mom_grid_size_+1,mom_grid_size_);
        R_pp = gsl_matrix_get(R_result,2*mom_grid_size_+1,2*mom_grid_size_+1);
        //std::cout << "R-matrix elements" << std::endl;
        //std::cout << R_mm << " " << R_mp << " " << R_pp << " " << std::endl;

        // Compute phase shifts in BB convention in radians
        phase_shifts.epsilon = atan(2.0*R_mp/(R_mm-R_pp))/2.0;
        phase_shifts.delta_p = atan((-rho/2.0)*(R_mm + R_pp - (R_mm - R_pp)/gsl_sf_cos(2*phase_shifts.epsilon)));
        phase_shifts.delta_m = atan((-rho/2.0)*(R_mm + R_pp + (R_mm - R_pp)/gsl_sf_cos(2*phase_shifts.epsilon)));
        //std::cout << "Phase shifts in BB" << std::endl;
        //std::cout << phase_shifts.delta_m << " " << phase_shifts.delta_p << " " << phase_shifts.epsilon << " " << std::endl;

        phase_shifts.delta_uncoupled = 0;
        phase_shifts = BB_to_Stapp(phase_shifts);

        //std::cout << "T:" << std::endl;
        //gsl_matrix_complex* T = T_matrix_from_R_matrix(R_result,rho);

    } else 
    {
        // on-shell R-matrix is 1x1
        double R = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        //std::cout << "R= " << R << std::endl;
        //std::cout << "-rho*R= " << -rho*R << std::endl;
        
        // Compute phase shift in radians BB and Stapp is the same for uncoupled channels
        phase_shifts.epsilon = 0;
        phase_shifts.delta_p = 0;
        phase_shifts.delta_m = 0;

        phase_shifts.delta_uncoupled = atan(-rho*R);
    }

    

    // Delete temporary pointers
    gsl_vector_free(D_vector);
    gsl_matrix_free(R_result);

    gsl_matrix_free(F_matrix);
    gsl_matrix_free(inverse);
    gsl_permutation_free(perm);

    return phase_shifts; // This conversion just affects the coupled channels
}

/*
    This function returns the T on-shell T-matrix from the phase shifts phase_shifts
    which needs to be entered in the Stapp convention. If the channel is uncoupled the 
    matrix is of size 1x1 otherwise it is of size 2x2.
*/
//gsl_matrix* get_T_matrix(Phase_shifts_chn phase_shifts, qs::quantum_channel chn)
//{
    // Observe! The phase chifts needs to be in the Stapp convention!
//}

// --------------------------------------------------------
// Here are methods for sloving LS for the complex T-matrix
// --------------------------------------------------------

// The method follows Chapter 18 in Landau QM

gsl_vector_complex* LS_Solver::setup_D_vector_complex(double q_on_shell, bool coupled, double mu)
{
    // See eq (18.19) in Landau
    #ifdef ENABLE_DEBUG
        std::cerr << "setup_D_vector_complex()" << std::endl;
    #endif

    // This is the factor in front of the LS equation.
    // In some cases it is 2/pi when that is a factor in
    // the partial wave momentum base normalization.
    // For now, the code is written without this factor, it will affect the 
    // rho value when raltiong the T/R matrix to the S-matrix/phase-shifts.
    double fac = 1.0; // Depends on convention
    //std::cout << "2mu: " << 2*mu << std::endl;
    //std::cout << "g_on_chell" << q_on_shell << std::endl;
    gsl_vector_complex* D_vector;
    double q2_on_shell = q_on_shell*q_on_shell;

    if (coupled)
    {
        // If the channel is coupled D is of length 
        // 2*mom_grid_size_ + 2
        D_vector = gsl_vector_complex_alloc(2*mom_grid_size_ + 2);
    } else 
    {
        D_vector = gsl_vector_complex_alloc(mom_grid_size_ + 1);
    }

    for (int i = 0; i < mom_grid_size_; i++)
    {
        double p2 = p_grid_[i]*p_grid_[i];
        double el = (2.0*mu)*(fac)*w_grid_[i]*p2/(p2-q2_on_shell); // NO CUTOFF
        
        gsl_complex comp_el = gsl_complex_rect(el,0.0);
        if (coupled)
        {
            gsl_vector_complex_set(D_vector,i,comp_el);
            gsl_vector_complex_set(D_vector,i+ mom_grid_size_+1,comp_el);       
        } else
        {
            gsl_vector_complex_set(D_vector,i,comp_el);
        }
    }

    // Add the last element. The first part is the same as in the non-complex case.
    double sum = 0;
    for (int i=0; i < mom_grid_size_; i++)
    {
        sum += w_grid_[i]/(p_grid_[i]*p_grid_[i]-q2_on_shell);
    }
    double re_el = -(fac)*2.0*mu*q2_on_shell*sum; // NO CUTOFF
    //double im_el = -2*mu*q_on_shell; // PROBALY A FACTOR pi/2 missing
    double im_el = M_PI*mu*q_on_shell; // PROBALY A FACTOR pi/2 missing

    gsl_complex comp_el = gsl_complex_rect(re_el,im_el);

    if (coupled)
    {
        gsl_vector_complex_set(D_vector,mom_grid_size_,comp_el);
        gsl_vector_complex_set(D_vector,2*mom_grid_size_+1,comp_el);
    } else
    {
        gsl_vector_complex_set(D_vector,mom_grid_size_,comp_el);
    }
    return D_vector;
}

gsl_matrix_complex* LS_Solver::setup_F_matrix_complex(bool coupled, gsl_vector_complex* D_vector, gsl_matrix* V_mtx)
{
    #ifdef ENABLE_DEBUG
        std::cerr << "setup_F_matrix_complex()" << std::endl;
    #endif
    
  
    gsl_matrix_complex* F_mtx;
    if (coupled)
    {
        F_mtx = gsl_matrix_complex_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
    } else 
    {
        F_mtx = gsl_matrix_complex_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
    }

    // Construct F manually with two loops. These are very costly,
    // would like to have a better solution
    // F_ij = \delta_ij + D_j V_ij
    for (int i = 0; i < F_mtx->size1; i++)
    {
        for (int j = 0; j < F_mtx->size2; j++)
        {
            gsl_complex matrix_el = gsl_complex_mul(gsl_complex_rect(gsl_matrix_get(V_mtx, i,j),0),
                gsl_vector_complex_get(D_vector,j));
            if (i==j) {
                matrix_el = gsl_complex_add(matrix_el,gsl_complex_rect(1.0,0.0));
            }
            gsl_matrix_complex_set(F_mtx,i,j,matrix_el);
        }
    }
    
    return F_mtx;
}

/*
    This function is the same as solve_in_chn_R() with the difference that is works with complex types.
    The potential matrix can be complex and is solvec the LS equation in complex form, thereby
    obtaining the full of-shell somplex T-matrix.
*/
Phase_shifts_chn LS_Solver::solve_in_chn_T(double T_lab, qs::quantum_channel chn, gsl_matrix* pot_V_mtx)
{
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    double rho = (2.0/M_PI)*M_PI*q_on_shell*mu; 
    // Depends on convention conncted to the factor 'fac' in D.
    // The transformation of T-matrix elements becomes correct when
    // a factor of 2.0/M_PI is added...

    gsl_vector_complex* D_vector = setup_D_vector_complex(q_on_shell,chn.coupled,mu);

    gsl_matrix_complex* F_matrix = setup_F_matrix_complex(chn.coupled,D_vector,pot_V_mtx);
    
    /*
    print_matrix(pot_V_mtx);
    std::cout << "D-vector" << std::endl;
    print_vector_complex(D_vector);
    std::cout << std::endl << "F-matrix" << std::endl;
    print_matrix_complex(F_matrix);
    */

    // Solve matrix equation F*R = V

    // LU decompose
    gsl_permutation* perm = gsl_permutation_alloc(F_matrix->size1);
    int signum;
    gsl_linalg_complex_LU_decomp(F_matrix,perm,&signum);

    // Invert from LU decompusition
 
    gsl_matrix_complex* inverse = gsl_matrix_complex_alloc(F_matrix->size1,F_matrix->size2);
    gsl_linalg_complex_LU_invert(F_matrix,perm,inverse); // Some error here

    gsl_matrix_complex* T_result = gsl_matrix_complex_alloc(F_matrix->size1,F_matrix->size2);

    gsl_complex alpha = gsl_complex_rect(1.0,0.0);
    gsl_complex beta  = gsl_complex_rect(0.0,0.0);

    gsl_matrix_complex* pot_complex = gsl_matrix_complex_alloc(F_matrix->size1,F_matrix->size2);
    for (int i=0; i < F_matrix->size1; i++)
    {
        for (int j=0; j < F_matrix->size2; j++)
        {
            gsl_matrix_complex_set(pot_complex,i,j,gsl_complex_rect(gsl_matrix_get(pot_V_mtx,i,j),0.0));
        }
    }
 
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, inverse, pot_complex, beta, T_result); 


    Phase_shifts_chn phase_shifts;
    if (chn.coupled) 
    {
        gsl_complex T_pp,T_mm,T_mp;
        T_mm = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_));
        T_mp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,mom_grid_size_));
        T_pp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,2*mom_grid_size_+1));
        
        // Correct values, just to test the calculation of phase shifts from T-elements
        /*
        T_mm = gsl_complex_rect(3.806699169483742e-06,-1.4527989826039648e-05);
        T_mp = gsl_complex_rect(-1.0676319790660848e-07,4.312313423841913e-07);
        T_pp = gsl_complex_rect(2.1377153170276125e-07,-1.565749889370409e-08);
        */
        /*
        std::cout << "T-matrix elements" << std::endl;
        std::cout << GSL_REAL(T_mm) << "," << GSL_IMAG(T_mm) << std::endl;
        std::cout << GSL_REAL(T_mp) << "," << GSL_IMAG(T_mp) << std::endl;
        std::cout << GSL_REAL(T_pp) << "," << GSL_IMAG(T_pp) << std::endl;
        */
        // Know that this is correct
        // ----
        // Compute phase shifts in BB convention in radians
        gsl_complex epsilon = gsl_complex_mul(gsl_complex_rect(0.5,0.0),gsl_complex_arctan(gsl_complex_div(
            gsl_complex_mul(gsl_complex_rect(2.0,0.0),T_mp),gsl_complex_sub(T_mm,T_pp))));
        
        // Some helper variables
        gsl_complex tmp1 = gsl_complex_sub(gsl_complex_rect(1.0,0.0),gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0.0,1.0),gsl_complex_rect(rho,0.0)),
            gsl_complex_add(T_mm,T_pp)));
        gsl_complex tmp2 = gsl_complex_mul(gsl_complex_rect(0.0,2.0),gsl_complex_mul(gsl_complex_rect(rho,0.0),
            gsl_complex_div(T_mp,gsl_complex_sin(gsl_complex_mul(gsl_complex_rect(2.0,0.0),epsilon)))));
        
        gsl_complex delta_p = gsl_complex_mul(gsl_complex_rect(0.0,-0.5),gsl_complex_log(gsl_complex_add(tmp1,tmp2)));
        gsl_complex delta_m = gsl_complex_mul(gsl_complex_rect(0.0,-0.5),gsl_complex_log(gsl_complex_sub(tmp1,tmp2)));
        // ----

        /*
        std::cout << "Printing phase shifts in BB from T:" << std::endl;
        std::cout << "(" << GSL_REAL(delta_m) << "," << GSL_IMAG(delta_m) << ")\n";
        std::cout << "(" << GSL_REAL(delta_p) << "," << GSL_IMAG(delta_p) << ")\n";
        std::cout << "(" << GSL_REAL(epsilon) << "," << GSL_IMAG(epsilon) << ")\n";
        */

        phase_shifts.epsilon = GSL_REAL(epsilon);
        phase_shifts.delta_p = GSL_REAL(delta_p);
        phase_shifts.delta_m = GSL_REAL(delta_m);
        //std::cout << "Phase shifts in BB" << std::endl;
        //std::cout << phase_shifts.delta_m << " " << phase_shifts.delta_p << " " << phase_shifts.epsilon << " " << std::endl;

        phase_shifts.delta_uncoupled = 0;
        phase_shifts = BB_to_Stapp(phase_shifts);

    } else 
    {
        gsl_complex T = gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_);
        //std::cout << GSL_REAL(T) << std::endl;
        phase_shifts.epsilon = 0;
        phase_shifts.delta_p = 0;
        phase_shifts.delta_m = 0;

        gsl_complex tmp_c = gsl_complex_sub(gsl_complex_rect(1.0,0.0),gsl_complex_mul(gsl_complex_rect(0.0,M_PI*rho),T));
        phase_shifts.delta_uncoupled = GSL_REAL(gsl_complex_mul(gsl_complex_rect(0.0,-0.5),gsl_complex_log(tmp_c)));
    }

    gsl_vector_complex_free(D_vector);
    gsl_matrix_complex_free(T_result);

    gsl_matrix_complex_free(F_matrix);
    gsl_matrix_complex_free(inverse);
    gsl_permutation_free(perm);

    return phase_shifts;
}


gsl_matrix_complex* LS_Solver::T_matrix_from_R_matrix(const gsl_matrix* R_matrix,double rho)
{
    gsl_matrix_complex* T = gsl_matrix_complex_alloc(R_matrix->size1,R_matrix->size2);

    // Apply the formula
    // R = T + i*pi*T*\delta(E-H_0)R
    // in the momentum basis R is expressed (no 2/pi factor in normalization).
    // This gives T = R(1+i\rhoR)^{-1}.

    // Make the R-matrix complex
    gsl_matrix_complex* R_complex = gsl_matrix_complex_alloc(R_matrix->size1,R_matrix->size2);

    for (int i = 0; i < R_matrix->size1; i++)
    {
        for (int j = 0; j < R_matrix->size2; j++)
        {
            gsl_complex R_el = gsl_complex_rect(gsl_matrix_get(R_matrix, i,j),0);
            gsl_matrix_complex_set(R_complex,i,j,R_el);
        }
    }

    gsl_matrix_complex* temp_matrix = gsl_matrix_complex_alloc(R_matrix->size1,R_matrix->size2);

    gsl_matrix_complex_set_identity(temp_matrix);
    
    gsl_complex alpha = gsl_complex_rect(0.0,rho);
    gsl_complex beta  = gsl_complex_rect(1.0,0.0);
    
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,alpha,
        temp_matrix, R_complex,beta, temp_matrix);

    // Now tmp_matrix is (1+i\rho R)

    // Invert tmp_matrix
    // LU decompose
    gsl_permutation* perm = gsl_permutation_alloc(temp_matrix->size1);
    int signum;
    gsl_linalg_complex_LU_decomp(temp_matrix,perm,&signum);

    // Invert from LU decompusition
 
    gsl_matrix_complex* inverse = gsl_matrix_complex_alloc(temp_matrix->size1,temp_matrix->size2);
    gsl_linalg_complex_LU_invert(temp_matrix,perm,inverse); // Some error here

    gsl_complex alpha2 = gsl_complex_rect(1.0,0.0);
    gsl_complex beta2  = gsl_complex_rect(0.0,0.0);

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha2, R_complex, inverse, beta2, T); 

    gsl_matrix_complex_free(inverse);

    gsl_matrix_complex_free(R_complex);
    gsl_matrix_complex_free(temp_matrix);
    gsl_permutation_free(perm);

    gsl_complex T_pp,T_mm,T_mp;
    T_mm = gsl_matrix_complex_get(T,mom_grid_size_,mom_grid_size_);
    T_mp = gsl_matrix_complex_get(T,2*mom_grid_size_+1,mom_grid_size_);
    T_pp = gsl_matrix_complex_get(T,2*mom_grid_size_+1,2*mom_grid_size_+1);
    std::cout << GSL_REAL(T_mm) << "," << GSL_IMAG(T_mm) << std::endl;
    std::cout << GSL_REAL(T_mp) << "," << GSL_IMAG(T_mp) << std::endl;
    std::cout << GSL_REAL(T_pp) << "," << GSL_IMAG(T_pp) << std::endl;

    return T;
}
