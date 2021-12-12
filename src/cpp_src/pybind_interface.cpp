#include "pybind_interface.h"

nn_mwpc_interface::nn_mwpc_interface(int J)
{
    // Define constants
    
    // Make pre-computations

    J_ = J;
}

double* nn_mwpc_interface::compute_phase_shifts()
{

}

double nn_mwpc_interface::compute_obs()
{
    return (double)J_;
}
