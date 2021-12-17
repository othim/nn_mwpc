#include "pybind_interface.h"

nn_mwpc_interface::nn_mwpc_interface(const std::string& model_name, 
        int J_max_chn, double cutoff, bool pre_comp_pot, bool rel_corr)
{

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale_ = 100.0; // Scale of momenutm grid MeV
    int number_of_p_points_ = 100; // Number of momentum-grid points
    int ang_int_points_ = 96; // Number of points in angular integration
    int J_max_in_pot_ = 50; // Maximum J that is stored for L-polynomials
    double cutoff_ = cutoff; // Cutoff in LS-equation
    bool pre_comp_pot_ = pre_comp_pot; // If pre-computations should be made
    bool rel_corr_ = rel_corr;

    // For the quantum states
    int J_max = J_max_chn;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;
    // ---------------------------------
    // ---------------------------------

    // Initialize physics helpers
    ph::physics_helper_init();

    // Construct the quantum states and quantum channels
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
            Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    chns_ = get_channels(states, print);   

    // Make GL grid
    double* p_grid;
    double* w_grid;
    ph::gauss_legendre_inf_mesh(number_of_p_points_,scale_,&p_grid,&w_grid);
    
    // These are the pre-determined models
    if ("WPC_LO"==model_name)
    {
        // Choose terms in potential
        std::vector<std::string terms> terms;
        terms.push_back("OPEP");
        terms.push_back("C1S0");
        terms.push_back("C3S1");

        // Construct potential
        Pot_ = new Potential_mwpc(terms,ang_int_points_,p_grid,w_grid,
                number_of_p_points_,J_max_in_pot_,cutoff_);
        Pot_ext_ = nullptr;

        if (pre_comp_pot)
        {
            // Save potential
            for (auto chn : chns_)
            {
                Pot->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(chns_,number_of_p_points_, scale_
                cutoff_,rel_corr_);

    } else
    {
        std::cout << "Error, not a valid model_name" << std::endl;
    }
}

nn_mwpc_interface::initialize()
{
    // Delete all pointers
}

nn_mwpc_interface::~nn_mwpc_interface()
{

}
    
std::vector<double> nn_mwpc_interface::compute_observable(const std::string& name, 
        std::vector<double> angles, double T_lab, std::vector<double> LECs)
{
    // Compute phase shifts in all channels for the given energy
    
    // For the specified angles (not 90...) 
    std::vector<double> obs_ang;
    /*for (int i = 0; i < 1000; i++)
    {
        obs_ang.push_back((double)i);
    }*/

    return obs_ang;
}
