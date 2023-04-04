#include "LS_Solver.h"
//#define TIME
//#define ENABLE_DEBUG

// Constructor
LS_Solver::LS_Solver(unsigned int mom_grid_size,
    double* p_grid, double* w_grid, bool finite_grid)
{
    #ifdef ENABLE_DEBUG
        std::cout << "LS_Solver()" << std::endl;
    #endif
    // Init variables
    mom_grid_size_ = mom_grid_size;
    p_grid_        = p_grid;
    w_grid_        = w_grid;
    
    // These are important, since if the grid is finite the counterterm for the
    // integral needs to be added. It is important that the potential is
    // set to zero when the momenta are larger than finite_grid_max
    finite_grid_     = finite_grid;
    finite_grid_max_ = 0.0; // Default value
    if (finite_grid)
    {
        // The max-value of the grid is the last element
        finite_grid_max_ = p_grid[mom_grid_size_-1];
    }
}

// Destructor
LS_Solver::~LS_Solver()
{
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
   
    
    for (int i = 0; i < (int)mom_grid_size_; i++)
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
    for (int i=0; i < (int)mom_grid_size_; i++)
    {
        sum += w_grid_[i]/(p_grid_[i]*p_grid_[i]-q2_on_shell);
    }
    // This is the counterterm that comes fromt that the integral is not computed
    // to infinity for a finite grid. This part compenstates for this by
    // adding the last part of the integral analytically.
    if (finite_grid_)
    {
        double counterterm = std::atanh(q_on_shell/finite_grid_max_)/q_on_shell;
        sum += counterterm;
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
/*TODO: fix these functions
 
gsl_matrix_complex * my_diag_alloc_complex(const gsl_vector * X)
{
    gsl_matrix * mat = gsl_matrix_alloc(X->size, X->size);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_matrix_set_all(mat, 0.0); //or whatever number you like
    gsl_vector_memcpy(&diag.vector, X);
    return mat;
}

int my_scale_columns_complex(gsl_matrix_complex* a, const gsl_vector* b)
{
    gsl_matrix_complex* m_diag = my_diag_alloc_complex(b);
    
    gsl_matrix* r = gsl_matrix_alloc(a->size1, a->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, a, m_diag, 0.0, r);
    gsl_matrix_memcpy(a,r);
    gsl_matrix_free(m_diag);
    gsl_matrix_free(r);
    return 0;
}
*/
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
    my_scale_columns(V_copy, D_vector);

    // Add V D to F
    gsl_matrix_add(F_mtx,V_copy);
    gsl_matrix_free(V_copy);

    return F_mtx;
}

void print_matrix(gsl_matrix* matrix)
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

/*
    Converts the phase shifts in the BB convention to the Stapp
    convention. This version does not have a (known) numerical instability.
*/
Phase_shifts_chn BB_to_Stapp(Phase_shifts_chn ps)
{
    Phase_shifts_chn phases;
    /*double diff = ps.delta_m-ps.delta_p;
    
    phases.epsilon = 0.5*asin((double)sin((double)2.0*ps.epsilon)* sin((double)diff));
    phases.delta_m = 0.5*(ps.delta_p + ps.delta_m + asin(tan(2.0*phases.epsilon)/tan(2.0*ps.epsilon)));
    phases.delta_p = 0.5*(ps.delta_p + ps.delta_m - asin(tan(2.0*phases.epsilon)/tan(2.0*ps.epsilon)));
    */
    phases.delta_uncoupled = ps.delta_uncoupled;


    // New from Boris code. There is a numerical instability in the above
    // equations
    // ----
    double cos2eps = std::cos(ps.epsilon)*std::cos(ps.epsilon);
    double cos_2delta_plus  = std::cos(2.*ps.delta_p);
    double sin_2delta_plus  = std::sin(2.*ps.delta_p);
    double cos_2delta_minus = std::cos(2.*ps.delta_m);
    double sin_2delta_minus = std::sin(2.*ps.delta_m);
    
    double aR, aI, tmp;
    aR = cos2eps*cos_2delta_minus + (1.-cos2eps)*cos_2delta_plus;
    aI = cos2eps*sin_2delta_minus + (1.-cos2eps)*sin_2delta_plus;
    
    phases.delta_m = 0.5*std::atan2(aI, aR);
    aR = cos2eps*cos_2delta_plus + (1.-cos2eps)*cos_2delta_minus;
    aI = cos2eps*sin_2delta_plus + (1.-cos2eps)*sin_2delta_minus;
    phases.delta_p = 0.5*std::atan2(aI, aR);
    
    tmp = 0.5*std::sin(2.0*ps.epsilon);
    aR  = tmp*(cos_2delta_minus - cos_2delta_plus);
    aI  = tmp*(sin_2delta_minus - sin_2delta_plus);
    tmp = (phases.delta_p + phases.delta_m);
    phases.epsilon = 0.5*std::asin(aI*cos(tmp) - aR*sin(tmp));
    // ----
    return phases;
}

Phase_shifts_chn BB_to_Stapp_2(Phase_shifts_chn ps,double tan_p, double tan_m, double tan_e)
{
    Phase_shifts_chn phases;
    double diff = ps.delta_m-ps.delta_p;
    std::cout << "--------" << std::endl;
    //std::cout << "Diff: " << diff << std::endl; 
    phases.epsilon = 0.5*asin(sin(2*ps.epsilon)* sin(diff));
    
    phases.delta_m = 0.5*(ps.delta_p + ps.delta_m + asin(tan(2*phases.epsilon)/tan(2*ps.epsilon)));
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
    #ifdef TIME
        std::clock_t start, end;
        start = std::clock(); 
    #endif    
    
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    // std::cout << "q_on_shell=" << q_on_shell << std::endl;
    // rho = 2*q*mu   
    // This is a convecntion dependent parameter that relates the 
    // R-matrix to the T/S matrices. rho will be different if a different 
    // normalization for the |klm> quantum states is choosen. To see the conventions
    // Used in this code see the README.md file.
    double rho = M_PI*q_on_shell*mu;

    //std::cout << "Potential" << std::endl;
    //print_matrix(pot_V_mtx);
   
    // Setup D-vector
    gsl_vector* D_vector = setup_D_vector(q_on_shell,chn.coupled,mu);
    //print_vector(D_vector);
    // Setup F-matrix
    gsl_matrix* F_matrix = setup_F_matrix(chn.coupled,D_vector,pot_V_mtx);
    
    #ifdef TIME
        end = std::clock();
        std::cout << "Time D,F: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << std::endl;
        start = std::clock();
    #endif  
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
    #ifdef TIME
        end = std::clock();
        std::cout << "Time invert: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << std::endl;
        start = std::clock();
    #endif

    Phase_shifts_chn phase_shifts;
    if (chn.coupled) 
    {
        // on-shell R-matrix is 2x2
        double R_pp,R_mm,R_mp;
        R_mm = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        R_mp = gsl_matrix_get(R_result,2*mom_grid_size_+1,mom_grid_size_);
        R_pp = gsl_matrix_get(R_result,2*mom_grid_size_+1,2*mom_grid_size_+1);

        // Compute phase shifts in BB convention in radians
        /*
        double tm = (R_mm-R_pp);
        double tp = (R_mm+R_pp);
        double x = 2.0*R_mp/tm;
        //phase_shifts.epsilon = atan2(2.0*R_mp,tm)/2.0;
        phase_shifts.epsilon = atan(2.0*R_mp/tm)/2.0;

        double rr = tm*(sqrt(1+x*x));

        phase_shifts.delta_m = atan((-rho/2.0)*(tp + rr));
        phase_shifts.delta_p = atan((-rho/2.0)*(tp - rr));

        phase_shifts.delta_uncoupled = 0;
        */

        // Compute phase shifts
        phase_shifts = BB_phases_from_R_coup(R_mm,R_pp,R_mp,rho);
        
        // Convert to Stapp convention
        phase_shifts = BB_to_Stapp(phase_shifts);
        
    } else 
    {
        // on-shell R-matrix is 1x1
        double R = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        
        /*
        // Compute phase shift in radians BB and Stapp is the same for uncoupled channels
        phase_shifts.epsilon = 0;
        phase_shifts.delta_p = 0;
        phase_shifts.delta_m = 0;

        phase_shifts.delta_uncoupled = atan(-rho*R);
        */

        phase_shifts = BB_phases_from_R_uncoup(R,rho);
    }

    #ifdef TIME
        end = std::clock();
        std::cout << "Time compute phase shifts: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << std::endl;
    #endif

    // Delete temporary pointers
    gsl_vector_free(D_vector);
    gsl_matrix_free(R_result);

    gsl_matrix_free(F_matrix);
    gsl_matrix_free(inverse);
    gsl_permutation_free(perm);

    return phase_shifts; // This conversion just affects the coupled channels
}

double* LS_Solver::solve_in_chn_R_Relem(double T_lab, qs::quantum_channel chn, gsl_matrix* pot_V_mtx)
{
    #ifdef ENABLE_DEBUG
        std::cerr << "solve_in_chn()" << std::endl;
    #endif
    // Compute reduced mass mu, which depends on the isospin-prijection.
    #ifdef TIME
        std::clock_t start, end;
        start = std::clock(); 
    #endif    
    
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    //std::cout << "q_on_shell=" << q_on_shell << std::endl;
    // rho = 2*q*mu   
    // This is a convecntion dependent parameter that relates the 
    // R-matrix to the T/S matrices. rho will be different if a different 
    // normalization for the |klm> quantum states is choosen. To see the conventions
    // Used in this code see the README.md file.
    
    //double rho = M_PI*q_on_shell*mu;

    //std::cout << "Potential" << std::endl;
    //print_matrix(pot_V_mtx);
   
    // Setup D-vector
    gsl_vector* D_vector = setup_D_vector(q_on_shell,chn.coupled,mu);
    //print_vector(D_vector);
    // Setup F-matrix
    gsl_matrix* F_matrix = setup_F_matrix(chn.coupled,D_vector,pot_V_mtx);
    
    #ifdef TIME
        end = std::clock();
        std::cout << "Time D,F: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << std::endl;
        start = std::clock();
    #endif  
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
    #ifdef TIME
        end = std::clock();
        std::cout << "Time invert: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << std::endl;
        start = std::clock();
    #endif
    
    double* Relem = new double[4];
    if (chn.coupled) 
    {
        // on-shell R-matrix is 2x2
        double R_pp,R_mm,R_mp;
        R_mm = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        R_mp = gsl_matrix_get(R_result,2*mom_grid_size_+1,mom_grid_size_);
        R_pp = gsl_matrix_get(R_result,2*mom_grid_size_+1,2*mom_grid_size_+1);
        //std::cout << "R-matrix elements" << std::endl;
        //std::cout << R_mm << " " << R_mp << " " << R_pp << " " << std::endl;

        Relem[0] = R_mm;
        Relem[1] = R_mp;
        Relem[2] = R_pp;
        Relem[3] = 0.0;
        
    } else 
    {
        // on-shell R-matrix is 1x1
        double R = gsl_matrix_get(R_result,mom_grid_size_,mom_grid_size_);
        Relem[0] = 0.0;
        Relem[1] = 0.0;
        Relem[2] = 0.0;
        Relem[3] = R;
        
    }


    // Delete temporary pointers
    gsl_vector_free(D_vector);
    gsl_matrix_free(R_result);

    gsl_matrix_free(F_matrix);
    gsl_matrix_free(inverse);
    gsl_permutation_free(perm);

    return Relem; // This conversion just affects the coupled channels
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

    for (int i = 0; i < (int)mom_grid_size_; i++)
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
    for (int i=0; i < (int)mom_grid_size_; i++)
    {
        sum += w_grid_[i]/(p_grid_[i]*p_grid_[i]-q2_on_shell);
    }
    // This is the counterterm that comes fromt that the integral is not computed
    // to infinity for a finite grid. This part compenstates for this by
    // adding the last part of the integral analytically.
    if (finite_grid_)
    {
        double counterterm = std::atanh(q_on_shell/finite_grid_max_)/q_on_shell;
        sum += counterterm;
    }
    
    
    double re_el = -(fac)*2.0*mu*q2_on_shell*sum;

    double im_el = M_PI*mu*q_on_shell;

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

gsl_vector_complex* LS_Solver::setup_G0_vector_complex(double q_on_shell, bool coupled, double mu)
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

    for (int i = 0; i < (int)mom_grid_size_; i++)
    {
        double p2 = p_grid_[i]*p_grid_[i];
        double el = (2.0*mu)/(q2_on_shell-p2); // NO CUTOFF
        
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
    for (int i=0; i < (int)mom_grid_size_; i++)
    {
        double p2 = p_grid_[i]*p_grid_[i];
        sum += w_grid_[i]/(q2_on_shell-p2);
    }
    // This is the counterterm that comes fromt that the integral is not computed
    // to infinity for a finite grid. This part compenstates for this by
    // adding the last part of the integral analytically.
    if (finite_grid_)
    {
        double counterterm = -std::atanh(q_on_shell/finite_grid_max_)/q_on_shell;
        sum += counterterm;
    }
    
    
    double re_el = -2.0*mu*q2_on_shell*sum;

    double im_el = -M_PI*mu*q_on_shell;

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

gsl_matrix_complex* LS_Solver::setup_F_matrix_complex(bool coupled, 
        gsl_vector_complex* D_vector, gsl_matrix* V_mtx)
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
    for (int i = 0; i < (int)F_mtx->size1; i++)
    {
        for (int j = 0; j < (int)F_mtx->size2; j++)
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

gsl_matrix_complex* LS_Solver::setup_F_matrix_complex(bool coupled, 
        gsl_vector_complex* D_vector, gsl_matrix_complex* V_mtx)
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
    for (int i = 0; i < (int)F_mtx->size1; i++)
    {
        for (int j = 0; j < (int)F_mtx->size2; j++)
        {
            gsl_complex matrix_el = gsl_complex_mul(gsl_matrix_complex_get(V_mtx, i,j),
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
    gsl_matrix_complex* T_result = LS_Solver::solve_in_chn_T_fullT(T_lab, 
            chn, pot_V_mtx);
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    double rho = M_PI*q_on_shell*mu; 
    // Depends on convention conncted to the factor 'fac' in D.
    // The transformation of T-matrix elements becomes correct when
    // a factor of 2.0/M_PI is added...

    Phase_shifts_chn phase_shifts;
    if (chn.coupled) 
    {
        gsl_complex T_pp,T_mm,T_mp;
        //T_mm = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_));
        //T_mp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,mom_grid_size_));
        //T_pp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,2*mom_grid_size_+1));
        
        T_mm = gsl_complex_mul(gsl_complex_rect(1.0,0),gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_));
        T_mp = gsl_complex_mul(gsl_complex_rect(1.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,mom_grid_size_));
        T_pp = gsl_complex_mul(gsl_complex_rect(1.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,2*mom_grid_size_+1));
        
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

        gsl_complex tmp_c = gsl_complex_sub(gsl_complex_rect(1.0,0.0),gsl_complex_mul(gsl_complex_rect(0.0,2.0*rho),T));
        phase_shifts.delta_uncoupled = GSL_REAL(gsl_complex_mul(gsl_complex_rect(0.0,-0.5),gsl_complex_log(tmp_c)));
    }

    gsl_matrix_complex_free(T_result);

    return phase_shifts;
}


gsl_matrix_complex* LS_Solver::solve_in_chn_T_fullT(double T_lab, 
        qs::quantum_channel chn, gsl_matrix* pot_V_mtx)
{
    gsl_matrix_complex* pot_complex = 
        gsl_matrix_complex_alloc(pot_V_mtx->size1,pot_V_mtx->size2);
    for (int i=0; i < (int)pot_complex->size1; i++)
    {
        for (int j=0; j < (int)pot_complex->size2; j++)
        {
            gsl_matrix_complex_set(pot_complex,i,j,gsl_complex_rect(gsl_matrix_get(pot_V_mtx,i,j),0.0));
        }
    }
    return solve_in_chn_T_fullT(T_lab,chn,pot_complex);
}

gsl_matrix_complex* LS_Solver::solve_in_chn_T_fullT(double T_lab, 
        qs::quantum_channel chn, gsl_matrix_complex* pot_complex)
{
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    //double rho = (2.0/M_PI)*M_PI*q_on_shell*mu; 
    // Depends on convention conncted to the factor 'fac' in D.
    // The transformation of T-matrix elements becomes correct when
    // a factor of 2.0/M_PI is added...

    gsl_vector_complex* D_vector = 
        setup_D_vector_complex(q_on_shell,chn.coupled,mu);
    gsl_matrix_complex* F_matrix = 
        setup_F_matrix_complex(chn.coupled,D_vector,pot_complex);
    
    // LU decompose
    gsl_permutation* perm = gsl_permutation_alloc(F_matrix->size1);
    int signum;
    gsl_linalg_complex_LU_decomp(F_matrix,perm,&signum);

    // Invert from LU decompusition
 
    gsl_matrix_complex* inverse = gsl_matrix_complex_alloc(F_matrix->size1,F_matrix->size2);
    gsl_linalg_complex_LU_invert(F_matrix,perm,inverse);

    gsl_matrix_complex* T_result = gsl_matrix_complex_alloc(F_matrix->size1,F_matrix->size2);

    gsl_complex alpha = gsl_complex_rect(1.0,0.0);
    gsl_complex beta  = gsl_complex_rect(0.0,0.0);

 
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, inverse, pot_complex, beta, T_result); 

    gsl_matrix_complex_free(pot_complex);
    gsl_vector_complex_free(D_vector);
    gsl_matrix_complex_free(F_matrix);
    gsl_matrix_complex_free(inverse);
    gsl_permutation_free(perm);
    return T_result;
}

gsl_matrix_complex* LS_Solver::solve_in_chn_T_fullT_weights(double T_lab, 
        qs::quantum_channel chn, gsl_matrix_complex* pot_complex,
        gsl_matrix_complex* G0)
{
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    // Compute K = (1 - V*G0)
    gsl_matrix_complex* id = 
        gsl_matrix_complex_alloc(pot_complex->size1,pot_complex->size2);
    gsl_matrix_complex_set_identity(id);
    gsl_complex alpha = gsl_complex_rect(-1.0,0.0);
    gsl_complex beta  = gsl_complex_rect(1.0,0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, pot_complex, 
            G0, beta, id); 
    
    // Invert id=K now
    gsl_permutation* perm = gsl_permutation_alloc(id->size1);
    int signum;
    gsl_linalg_complex_LU_decomp(id,perm,&signum);

    // Invert from LU decompusition
    gsl_matrix_complex* K_inv = 
        gsl_matrix_complex_alloc(pot_complex->size1,pot_complex->size2);
    gsl_linalg_complex_LU_invert(id,perm,K_inv);

    // Multiply K_inv*V
    gsl_matrix_complex* T_result = 
        gsl_matrix_complex_alloc(pot_complex->size1,pot_complex->size2);
    alpha = gsl_complex_rect(1.0,0.0);
    beta  = gsl_complex_rect(0.0,0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, K_inv, 
            pot_complex, beta, T_result); 

    gsl_permutation_free(perm);
    gsl_matrix_complex_free(id);
    gsl_matrix_complex_free(K_inv);
    
    return T_result;
}

std::complex<double>* LS_Solver::solve_in_chn_T_Telem(double T_lab, 
        qs::quantum_channel chn, gsl_matrix* pot_V_mtx)
{

    gsl_matrix_complex* T_result = LS_Solver::solve_in_chn_T_fullT(T_lab, 
            chn, pot_V_mtx);

    std::complex<double>* T = new std::complex<double>[4];
    for (int  i=0; i < 4; i++)
    {
        T[i] = (std::complex<double>)0.0;
    }

    if (chn.coupled) 
    {
        gsl_complex T_pp,T_mm,T_mp;
        T_mm = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_));
        T_mp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,mom_grid_size_));
        T_pp = gsl_complex_mul(gsl_complex_rect(M_PI/2.0,0),gsl_matrix_complex_get(T_result,2*mom_grid_size_+1,2*mom_grid_size_+1));
          
        T[0].real(GSL_REAL(T_mm));
        T[0].imag(GSL_IMAG(T_mm));

        T[1].real(GSL_REAL(T_mp));
        T[1].imag(GSL_IMAG(T_mp));

        T[2].real(GSL_REAL(T_pp));
        T[2].imag(GSL_IMAG(T_pp));

    } else 
    {
        gsl_complex Tt = gsl_matrix_complex_get(T_result,mom_grid_size_,mom_grid_size_);
        //std::cout << GSL_REAL(T) << std::endl;
        T[3].real(GSL_REAL(Tt));
        T[3].imag(GSL_IMAG(Tt));

    }

    gsl_matrix_complex_free(T_result);
    
    return T;
}

std::complex<double>* LS_Solver::T_matrix_from_R_matrix(double Rmm, double Rmp, 
        double Rpp, double mu, double q_on_shell)
{
    // In the basis convention without 2/pi factor the relation between the 
    // R and T matrix in the partial wave basis reads
    // R = T + i*pi*\mu*q(E) R T
    // where there is a matrix product over l'l implicit.
    
    std::cout << Rmm << "  " << Rmp << "  " << Rpp << std::endl;
    // Make the R-matrix complex
    const int size = 2;
    gsl_matrix_complex* R_complex = gsl_matrix_complex_alloc(size,size);
    
    // Set the matrix elements
    double f = 1.0;
    gsl_matrix_complex_set(R_complex,0,0,gsl_complex_rect(f*Rmm,0));
    gsl_matrix_complex_set(R_complex,0,1,gsl_complex_rect(1.0*f*Rmp,0));
    gsl_matrix_complex_set(R_complex,1,0,gsl_complex_rect(1.0*f*Rmp,0));
    gsl_matrix_complex_set(R_complex,1,1,gsl_complex_rect(f*Rpp,0));

    std::cout << "R_complex:" << std::endl;
    ph::print_m(R_complex);
    
    // Make the other matrix
    gsl_matrix_complex* temp_matrix = gsl_matrix_complex_alloc(R_complex->size1
            ,R_complex->size2);
    gsl_matrix_complex_set_identity(temp_matrix);
    
    double fac = M_PI*mu*q_on_shell;
    gsl_complex alpha = gsl_complex_rect(0.0,fac);
    gsl_complex beta  = gsl_complex_rect(1.0,0.0);
    
    // temp_matrix = 1 + fac*R
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,alpha,
        temp_matrix, R_complex,beta, temp_matrix);

    // Invert tmp_matrix
    // LU decompose
    gsl_permutation* perm = gsl_permutation_alloc(temp_matrix->size1);
    int signum;
    gsl_linalg_complex_LU_decomp(temp_matrix,perm,&signum);

    // Invert from LU decompusition
    gsl_matrix_complex* inverse = 
        gsl_matrix_complex_alloc(temp_matrix->size1,temp_matrix->size2);
    gsl_linalg_complex_LU_invert(temp_matrix,perm,inverse);

    std::cout << "inverse:" << std::endl;
    ph::print_m(inverse);

    // Performs inverse <- inverse*R_complex
    gsl_complex alpha2 = gsl_complex_rect(1.0,0.0);
    gsl_complex beta2  = gsl_complex_rect(0.0,0.0);
    gsl_matrix_complex* tmp = 
        gsl_matrix_complex_alloc(temp_matrix->size1,temp_matrix->size2);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha2,inverse,R_complex, beta2, 
            tmp); 
    
    std::cout << "tmp:" << std::endl;
    ph::print_m(tmp);
    // Now inverse is the T-matrix
    
    std::complex<double>* T = new std::complex<double>[4];
    
    gsl_complex T_pp,T_mm,T_mp;
    T_mm = gsl_matrix_complex_get(tmp,0,0);
    T_mp = gsl_matrix_complex_get(tmp,1,0);
    T_pp = gsl_matrix_complex_get(tmp,1,1);
    
    T[0].real(GSL_REAL(T_mm));
    T[0].imag(GSL_IMAG(T_mm));

    T[1].real(GSL_REAL(T_mp));
    T[1].imag(GSL_IMAG(T_mp));

    T[2].real(GSL_REAL(T_pp));
    T[2].imag(GSL_IMAG(T_pp));

    std::cout << GSL_REAL(T_mm) << "," << GSL_IMAG(T_mm) << std::endl;
    std::cout << GSL_REAL(T_mp) << "," << GSL_IMAG(T_mp) << std::endl;
    std::cout << GSL_REAL(T_pp) << "," << GSL_IMAG(T_pp) << std::endl;
    
    gsl_matrix_complex_free(inverse);
    gsl_matrix_complex_free(R_complex);
    gsl_matrix_complex_free(temp_matrix);
    gsl_matrix_complex_free(tmp);
    gsl_permutation_free(perm);

    return T;
}

Phase_shifts_chn LS_Solver::BB_phases_from_R_coup(double R_mm, double R_pp, 
        double R_mp, double rho)
{
    Phase_shifts_chn phase_shifts;

    double tm = (R_mm-R_pp);
    double tp = (R_mm+R_pp);
    double x = 2.0*R_mp/tm;
    //phase_shifts.epsilon = atan2(2.0*R_mp,tm)/2.0;
    phase_shifts.epsilon = atan(2.0*R_mp/tm)/2.0;

    double rr = tm*(sqrt(1+x*x));

    phase_shifts.delta_m = atan((-rho/2.0)*(tp + rr));
    phase_shifts.delta_p = atan((-rho/2.0)*(tp - rr));

    phase_shifts.delta_uncoupled = 0;

    return phase_shifts;
}
    
Phase_shifts_chn LS_Solver::BB_phases_from_R_uncoup(double R, double rho)
{
    Phase_shifts_chn phase_shifts;

    // Compute phase shift in radians BB and Stapp is the same for uncoupled channels
    phase_shifts.epsilon = 0;
    phase_shifts.delta_p = 0;
    phase_shifts.delta_m = 0;

    phase_shifts.delta_uncoupled = atan(-rho*R);

    return phase_shifts;
}
