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
    std::cout << "Done!" << std::endl;
}
nn_mwpc_dwb_interface::~nn_mwpc_dwb_interface()
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

