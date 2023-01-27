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
}
    
void nn_mwpc_dwb_interface::solve_exact_full_T()
{

}

void nn_mwpc_dwb_interface::solve_LO_full_T()
{

}

void nn_mwpc_dwb_interface::solve_DWBA_full_T(int order)
{

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

