#include "LS_Solver.h"


// Constructor
LS_Solver::LS_Solver()
{

}

// Destructor
LS_Solver::~LS_Solver()
{

}


gsl_vector* setup_D_vector(double q_on_shell, bool coupled, double mu)
{
    // See eq (18.19) in Landau

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
        D_vector = gsl_vector_alloc(mom_grid_size_ + 1);
    }
    // Channel is uncoupled D is of length mom_grid_size_ + 1
    
    
    for (int i = 0; i < mom_grid_size_; i++)
    {
        double p2 = p_grid_[i]*p_grid_[i];
        double el = (2.0*mu)*(2.0/M_PI)*w_grid_[i]*p2/(p2-q2_on_shell);

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
    double el = -(2.0/M_PI)*2.0*mu*q2_on_shell*sum;

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

gsl_matrix* setup_F_matrix(bool coupled, gsl_vector* D_vector, gsl_matrix* V_mtx)
{
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

Phase_shifts_chn LS_Solver::solve_in_chn(double q_on_shell, qs::quantum_channel chn, bool rel_correction, bool cutoff_on)
{
    // Compute reduced mass mu
    double mu =     

    gsl_matrix* pot_V_mtx;
    // Get potential matrix with the correct on-shell momentum
    // The size of the matrix depends on if the channels is coupled or not
    pot_V_->get_get_saved_matrix(q_on_shell,chn,rel_correction,cutoff_on, pot_V_mtx);
    // Setup D-vector
    gsl_vector* D_vector;
    D_vector = setup_D_vector(q_on_shell,coupled,mu);
    // Setup F-matrix
    gsl_matrix* F_matrix;
    F_matrix = setup_F_matrix(q_on_shell,pot_V_);
    // Solve the matrix equation F*R = V
        
    // LU decompose
    gsl_permutation* perm;
    int* signum;
    gsl_linalg_LU_decomp(F_matrix,perm,signum);

    // Invert from LU decompusition
    gsl_vector* inverse;
    gsl_linalg_LU_invert(F_matrix,perm,V_last_col,inverse);
    
    gsl_matrix* R_result;

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse, pot_V_mtx, 0.0, R_result); 

    // Get the matrix elements corresponding to the on-shell R-matrix

    Phase_shifts_chn phase_shifts;

    if (coupled) 
    {
        // on-shell R-matrix is 2x2
        double R[2][2];
        R[0][0] = gsl_matrix_get_element(R_result,mom_grid_size_,mom_grid_size_);
        R[1][0] = gsl_matrix_get_element(R_result,2*mom_grid_size_+1,mom_grid_size_);
        R[0][1] = R[1][0];
        R[1][1] = gsl_matrix_get_element(R_result,2*mom_grid_size_+1,2*mom_grid_size_+1);

        // Compute phase shifts

        
    } else 
    {
        // on-shell R-matrix is 1x1
        double R = gsl_matrix_get_element(R_result,mom_grid_size_,mom_grid_size_);

        // Compute phase shift
    }


    // Delete temporary pointers
    gsl_vector_free(D_vector);
    gsl_vector_free(R_result);

    gsl_matrix_free(pot_V_mtx;
    gsl_matrix_free(F_matrix);

    gsl_permutation_free(perm);
}