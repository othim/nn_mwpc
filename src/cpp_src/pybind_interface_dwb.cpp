#include "pybind_interface_dwb.h"

#ifdef PYBIND

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"
    
#endif


nn_mwpc_dwb_interface::nn_mwpc_dwb_interface(const std::string& model_name, 
        int J_max_chn, double cutoff, int cut_pow, 
        bool sharp_cutoff, bool pre_comp_pot, bool rel_corr,
        int number_of_p_points,bool finite_grid,bool inc_weights_in_pot,
        bool cut_on_shell)
{
    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    scale_ = 100.0; // Scale of momenutm grid MeV (100)
    number_of_p_points_ = number_of_p_points; // Number of momentum-grid points (60)
    ang_int_points_ = 76; // Number of points in angular integration
    J_max_in_pot_ = 50; // Maximum J that is stored for L-polynomials
    cutoff_ = cutoff; // Cutoff in LS-equation
    cut_pow_ = cut_pow;
    sharp_cutoff_ = sharp_cutoff;
    pre_comp_pot_ = pre_comp_pot; // If pre-computations should be made
    rel_corr_ = rel_corr;
    finite_grid_ = finite_grid;
    finite_grid_max_ = 0.0;
    inc_weights_in_pot_ = inc_weights_in_pot;
    cut_on_shell_ = cut_on_shell;
    
    // For the quantum states
    int J_max = J_max_chn;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = false;
    // ---------------------------------
    // ---------------------------------
    
    // Initialize physics helpers
    ph::physics_helpers_init();
    
    // Make GL grid
    if (!finite_grid)
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points_,scale_,&p_grid_,&w_grid_);
    } else 
    {
        // Make GL-grid that is finite. Since we use a sharp cutoff this grid is fine
        double sharp_cut_add = 300.0;
        finite_grid_max_ = cutoff_+sharp_cut_add;
        ph::gauss_legendre_finite_mesh(number_of_p_points_,0,finite_grid_max_,&p_grid_,&w_grid_);
    }
    
    // Construct the quantum states
    if (print) {
        std::cout << "Constructing quantum states..." << std::endl;
    }
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    if (print) {
        std::cout << "Contruction scattering channels..." << std::endl;
    }
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    


    // Construct a LS_Solver
    
    
    // Constructing the potentials 
    

}


nn_mwpc_dwb_interface::~nn_mwpc_dwb_interface()
{
    delete LS_Solver;

    // TODO: also delete potentals and other stuff

    delete[] p_grid;
    delete[] w_grid;
}
    
gsl_matrix_complex* nn_mwpc_dwb_interface::solve_exact_pot_sum_full_T(
        double T_lab, qs::quantum_channel chn)
{
    // Get G0 and potentials for the correct on-shell point
    gsl_matrix_complex* G0, VI, VII;
    get_G0_and_potentials(T_lab, chn, G0, VI, VII);

    // Solve LS equation for SUM of potentials, NOTE: VII is now VI+VII!!!
    gsl_matrix_complex_add(VII,VI);
    
    gsl_matrix_complex* T_tmp = 
        solver.solve_in_chn_T_fullT_weights(Tl,chn,VII,G0);
    
    gsl_matrix_complex_free(VI);
    gsl_matrix_complex_free(VII);
    gsl_matrix_complex_free(G0);

    return T_tmp;
}

gsl_matrix_complex* nn_mwpc_dwb_interface::solve_full_TI(double T_lab, 
        qs::quantum_channel chn)
{

    return T_tmp;
}

gsl_matrix_complex* nn_mwpc_dwb_interface::solve_DWBA_full_TI(double T_lab, 
        qs::quantum_channel chn, int order)
{
    // Get G0 and potentials for the correct on-shell point
    gsl_matrix_complex* G0, VI, VII;
    get_G0_and_potentials(T_lab, chn, G0, VI, VII);

    // Solve for TI
    gsl_matrix_complex* TI = 
        solver.solve_in_chn_T_fullT_weights(Tl,chn,VI,G0);
    
    // Solve DWB series
    gsl_matrix_complex* T_DWBA = dwba::pw_T_DWBA(order,TI,VI,VII,G0);
    
    // Delete and return
    gsl_matrix_complex_free(VI);
    gsl_matrix_complex_free(VII);
    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(TI);
    return T_tmp;
}

void get_G0_and_potentials(double T_lab, qs::quantum_channel chn,
        gsl_matrix_complex** G0, gsl_matrix_complex** VI,
        gsl_matrix_complex** VII)
{
    // Get on-shell momentum and reduced mass
    double q_on_shell,mu;
    LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
    
    // Get G0 vector matrix
    gsl_vector_complex* prop_vec = 
        solver.setup_G0_vector_complex(q_on_shell,chn.coupled,mu);
    // Make it to a diagonal matrix
    gsl_matrix_complex* G0_loc = 
        gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    matrix_from_vector(G0,prop_vec);
    gsl_vector_complex_free(prop_vec);

    // Get the sum of VI and VII with weights
    gsl_matrix_complex* VI_loc = 
        VI_complex_weights.get_saved_matrix(q_on_shell, chn, REL_CORR);

    gsl_matrix_complex* VII_loc = 
        VII_complex_weights.get_saved_matrix(q_on_shell, chn, REL_CORR);
    
    // Set the value of the external pointers
    *G0  = G0_loc;
    *VI  = VI_loc;
    *VII = VII_loc;
}


void nn_mwpc_dwb_interface::solve_DWBA_PC_full_T(int order)
{

}

void nn_mwpc_dwb_interface::print_LECs_in_use()
{
}

void nn_mwpc_dwb_interface::print_LEC_values()
{
}


double nn_mwpc_dwb_interface::get_on_shell_momentum(double T_lab)
{
}

double nn_mwpc_dwb_interface::get_scale()
{
}

int nn_mwpc_dwb_interface::get_ang_int_points()
{
}

int nn_mwpc_dwb_interface::get_momentum_grid_points()
{
}

int nn_mwpc_dwb_interface::get_chn_len()
{
} 

std::string nn_mwpc_dwb_interface::get_chn_LS_term(int chn_number)
{
}
double nn_mwpc_dwb_interface::get_gA()
{
}

double nn_mwpc_dwb_interface::get_fpi()
{
}

double nn_mwpc_dwb_interface::get_mpi()
{
}

double nn_mwpc_dwb_interface::get_Mp()
{
}

double nn_mwpc_dwb_interface::get_Mn()
{
}

