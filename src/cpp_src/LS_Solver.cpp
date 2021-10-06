#include "LS_Solver.h"

//#define ENABLE_DEBUG

// Constructor
LS_Solver::LS_Solver(std::vector<qs::quantum_channel> channels, Potential_mwpc* pot_V, unsigned int mom_grid_size,
    double mom_grid_scale, bool cutoff_enabled, double cutoff_Lambda, bool relcorr_enabled)
{
    #ifdef ENABLE_DEBUG
        std::cout << "LS_Solver()" << std::endl;
    #endif
    // Init variables
    pot_V_ = pot_V;
    channels_ = channels;
    mom_grid_size_ = mom_grid_size;
    cutoff_enabled_ = cutoff_enabled;
    cutoff_Lambda_ = cutoff_Lambda;
    relcorr_enabled_ = relcorr_enabled;
    
    // Make GL-grid
    gauss_legendre_inf_mesh(mom_grid_size_,mom_grid_scale,&p_grid_,&w_grid_);
}

// Destructor
LS_Solver::~LS_Solver()
{
    delete p_grid_;
    delete w_grid_;
}

/*
    This function creates a GL-grid on the interval [0,\infty) with the scale <scale>
*/
void LS_Solver::gauss_legendre_inf_mesh(unsigned int Numper_of_points, double scale,double** p,double** w)
{
    #ifdef ENABLE_DEBUG
        std::cout << "gauss_legendre_inf_mesh()" << std::endl;
    #endif
    // Make grid from -1 to 1

    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, Numper_of_points, -1.0, 1.0, 0, 0);

    double* p_grid = gsl_integration_fixed_nodes(int_ang_);
    double* w_grid = gsl_integration_fixed_weights(int_ang_);

    // Make transformation
    double pi_4 = M_PI/4.0;

    for (int i = 0; i < Numper_of_points; i++)
    {
        double x = p_grid[i];
        p_grid[i] = scale*tan(pi_4*(x+1));
        w_grid[i] = (scale*pi_4/(cos(pi_4*(x+1))*cos(pi_4*(x+1))))*w_grid[i];
    }

    // ???
    *p = p_grid;
    *w = w_grid;
   
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
        //double cutoff_regulator = exp(-gsl_pow_uint(p_grid_[i]/cutoff_Lambda_,6));
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
    
    // Scale V by D
    gsl_matrix_scale_columns(V_copy, D_vector);
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

void get_mu_q_on_shell(double T_lab, qs::quantum_channel chn, double* mu, double* q_on_shell)
{
    if (chn.tz == -1)
    {
        *mu = constants::Mn/2.0; // nn
        *q_on_shell = sqrt(*mu*T_lab);
    } else if (chn.tz == 0)
    {
        *mu = constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np
        *q_on_shell = sqrt(constants::Mp*constants::Mp*T_lab*(T_lab + 2.0*constants::Mn)/
            ((constants::Mp + constants::Mn)*(constants::Mp + constants::Mn) + 2.0*T_lab*constants::Mp));
    } else if (chn.tz == 1)
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
Phase_shifts_chn LS_Solver::solve_in_chn_R(double T_lab, qs::quantum_channel chn, bool rel_correction,bool get_saved_potential)
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

    
    // Get potential matrix with the correct on-shell momentum
    // The size of the matrix depends on if the channels is coupled or not
    gsl_matrix* pot_V_mtx;
    if (get_saved_potential) {
        pot_V_mtx = pot_V_->get_saved_matrix(q_on_shell,chn,rel_correction);
    } else {
        pot_V_mtx = pot_V_->get_matrix(q_on_shell,chn,rel_correction);
    }
    
    
    //std::cout << "Potential" << std::endl;
    //print_matrix(pot_V_mtx);
   
    // Setup D-vector
    gsl_vector* D_vector;
    D_vector = setup_D_vector(q_on_shell,chn.coupled,mu);
    
    // Setup F-matrix
    gsl_matrix* F_matrix;
    F_matrix = setup_F_matrix(chn.coupled,D_vector,pot_V_mtx);
    
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
        //std::cout << R_mm << " " << R_mp << " " << R_pp << " " << std::endl;

        // Compute phase shifts in BB convention in radians
        phase_shifts.epsilon = atan(2.0*R_mp/(R_mm-R_pp))/2.0;
        phase_shifts.delta_p = atan((-rho/2.0)*(R_mm + R_pp - (R_mm - R_pp)/gsl_sf_cos(2*phase_shifts.epsilon)));
        phase_shifts.delta_m = atan((-rho/2.0)*(R_mm + R_pp + (R_mm - R_pp)/gsl_sf_cos(2*phase_shifts.epsilon)));
        //std::cout << "Phase shifts in BB" << std::endl;
        //std::cout <<  std::setprecision(16) << phase_shifts.delta_m << " " << phase_shifts.delta_p << " " << phase_shifts.epsilon << " " << std::endl;

        phase_shifts.delta_uncoupled = 0;
        phase_shifts = BB_to_Stapp(phase_shifts);
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

    gsl_matrix_free(pot_V_mtx);
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
        //double cutoff_regulator = exp(-gsl_pow_uint(p_grid_[i]/cutoff_Lambda_,6));
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
    double im_el = -2*mu*q_on_shell; // PROBALY A FACTOR pi/2 missing

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

gsl_matrix_complex* LS_Solver::setup_F_matrix_complex(bool coupled, gsl_vector_complex* D_vector, gsl_matrix_complex* V_mtx)
{
    // F_ij = \delta_ij + D_j V_ij

    // Must be able to scale columns
    // Do scaling manually since the potential matrix
    // Anyhow needs to be reconstructed (this will be slow)

}


/*
    This function is the same as solve_in_chn_R() with the difference that is works with complex types.
    The potential matrix can be complex and is solvec the LS equation in complex form, thereby
    obtaining the full of-shell somplex T-matrix.
*/
Phase_shifts_chn LS_Solver::solve_in_chn_T(double T_lab, qs::quantum_channel chn, bool rel_correction,bool get_saved_potential)
{
    // Compute reduced mass mu, which depends on the isospin-prijection.
    double mu;
    double q_on_shell;
    get_mu_q_on_shell(T_lab,chn,&mu,&q_on_shell);

    gsl_vector_complex* D_vector = setup_D_vector_complex(q_on_shell,chn.coupled,mu);

    double rho = M_PI*q_on_shell*mu; // Depends on convention conncted to the factor 'fac' in D.

    gsl_matrix_complex* pot_V_mtx;
    if (get_saved_potential) {
        //pot_V_mtx = pot_V_->get_saved_matrix(q_on_shell,chn,rel_correction);
    } else {
        //pot_V_mtx = pot_V_->get_matrix(q_on_shell,chn,rel_correction);
    }

        

}