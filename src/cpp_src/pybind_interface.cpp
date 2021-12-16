#include "pybind_interface.h"

nn_mwpc_interface::nn_mwpc_interface(const std::string& model_name, double scale_gl_grid, 
        int gl_grid_size, int ang_int_size, int J_max_in_pot, int J_max_chn)
{
    // Construct potential
    
    // Construct LS Solver
    
    // Make pre-computations
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
    obs_ang.push_back(10.0);
    obs_ang.push_back(20.0);

    return obs_ang;
}
