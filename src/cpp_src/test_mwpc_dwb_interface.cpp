
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include "pybind_interface.h"
#include "pybind_interface_dwb.h"
#include "potential_ext.h"
#include "born_approx.h"
#include "potential_mwpc.h"
#include "physics_helpers.h"

void compute_spectrum(std::string name, std::string pot_name, 
        std::vector<double> LECs, std::vector<double> params, 
        nn_mwpc_dwb_interface& obj);


int main(int argc, char** argv)
{
    std::cout << "\033[1;31mbold FAILED\033[0m\n" << std::endl;
    std::cout << "\033[1;32mbold OK\033[0m\n" << std::endl;
    // ************************************************************************
    // ****** CONSTANTS TO CHANGE *********************************************
    // ************************************************************************
    double scale_              = 100; // Scale of momenutm grid MeV (100)
    int number_of_p_points_    = 100; // Number of momentum-grid points (60)
    int ang_int_points_        = 96; // Number of points in angular integration
    int J_max_in_pot_          = 50; // Maximum J that is stored for L-polynomials
    int J_max_chn_             = 2;
    double cutoff_             = 500; // Cutoff in LS-equation
    int cut_pow_               = 6;
    bool sharp_cutoff_         = false;
    bool rel_corr_             = true;
    bool finite_grid_          = false;
    bool finite_grid_max_      = 0.0; // Just default value
    bool cut_on_shell_         = true;
    bool print_                = true;
    
    
    // Construct an object
    
    nn_mwpc_dwb_interface obj = 
        nn_mwpc_dwb_interface(scale_,J_max_chn_,cutoff_, 
            cut_pow_,sharp_cutoff_,rel_corr_,number_of_p_points_, finite_grid_,
            cut_on_shell_,print_);
    
    std::cout << "Constructed object" << std::endl;
    
    /*
     * ***************
     * LO   
     * ***************
     */
    double C1S0 = -0.10768e-2;
    double C3S1 = -0.07172e-2;

    double gA   = 1.276;

    // Set parameters before saving!
    std::vector<double> params = {gA};
    // LECs need to be set after saving!
    std::vector<double> LECs = {C1S0,C3S1};
    
    std::string name = "LO";
    std::string pre_def_name = "WPC_LO";
    //compute_spectrum(name,pre_def_name,LECs,params, obj);
    
    /*
     * ***************
     * NLO   
     * ***************
     */
    C1S0 = -0.150533e-2;
    C3S1 = -0.1742e-2;
    
    double D1P1 = 0.849e-8;
    double D1S0 = 1.6926e-8;
    double D3P0 = 1.3085e-8;
    double D3P1 = -0.3409e-8;
    double D3P2 = -0.2011e-8;
    double D3S1 = -0.408e-8;
    double D_SD = -0.238e-8;

    gA   = 1.29;
    //gA   = 1.276;

    // Set parameters before saving!
    params = {gA};
    // LECs need to be set after saving!
    LECs = {C1S0,C3S1,D1P1,D1S0,D3P0,D3P1,D3P2,D3S1,D_SD};
    
    name = "NLO";
    pre_def_name = "WPC_NLO_SFR";
    //compute_spectrum(name,pre_def_name,LECs,params, obj);
    
    //C1S0 = -0.150533e-2;
    C1S0 = -0.0015062256361672;
    C3S1 = -0.0018431235715802;
    
    D1P1 = 0.849e-8;
    D1S0 = 1.6926e-8;
    D3P0 = 1.3085e-8;
    D3P1 = -0.3409e-8;
    D3P2 = -0.2011e-8;
    D3S1 = -0.0000000021789842;
    D_SD = -0.0000000026342609;

    gA   = 1.29;
    //gA   = 1.276;

    // Set parameters before saving!
    params = {gA};
    // LECs need to be set after saving!
    LECs = {C1S0,C3S1,D1P1,D1S0,D3P0,D3P1,D3P2,D3S1,D_SD};
    
    name = "NLO2";
    pre_def_name = "WPC_NLO_SFR";
    compute_spectrum(name,pre_def_name,LECs,params, obj);


    /*
     * ***************
     * N2LO   
     * ***************
     */
    C1S0 = -0.15387e-2;
    C3S1 = -0.1671e-2;
    
    D1P1 = 0.521e-8;
    D1S0 = 2.7442e-8;
    D3P0 = 1.2782e-8;
    D3P1 = -0.9378e-8;
    D3P2 = -0.68645e-8;
    D3S1 = +0.8738e-8;
    D_SD = -0.6899e-8;

    gA   = 1.29;

    double c1 = -0.69*1e-3;
    double c3 = -4.12*1e-3;
    double c4 = 5.35*1e-3;
    // Set parameters before saving!
    params = {gA};
    // LECs need to be set after saving!
    LECs = {C1S0,C3S1,D1P1,D1S0,D3P0,D3P1,D3P2,D3S1,D_SD,c1,c3,c4};
    name = "N2LO";
    pre_def_name = "WPC_N2LO_SFR";
    //compute_spectrum(name,pre_def_name,LECs,params, obj);
    
    //std::vector<std::complex<double>> T = obj.solve_DWBA_T_PC(25.0,0,2,"LO","NLO","NLO2","");

    //for (auto t : T)
    //{
    //    std::cout << std::real(t) << "," << std::imag(t) << std::endl;
    //}
    /*for (int q=1; q<300;q++)
    {
        std::cout << Term::L_DR((double)q,constants::mpi) << ", ";
        std::cout << Term::L_SFR((double)q,constants::mpi,1000000.0) << std::endl;
    }*/

    return 0;
}


void compute_spectrum(std::string name, std::string pot_name, std::vector<double> LECs, 
        std::vector<double> params, nn_mwpc_dwb_interface& obj)
{
    double lam_SFR = 700.0;
    obj.create_new_potential(name,pot_name,lam_SFR);
    
    obj.print_potential_info(name);

    obj.set_params_in_potential(name,params);
    obj.save_potential_decomposition(name);
    obj.set_LECs_in_potential(name,LECs);
    
    obj.print_potential_info(name);
 
    int chn_number = 3;
    std::vector<double> res = obj.compute_binding_energy(chn_number, true, name);
    std::cout << name << std::endl;
    for (auto a : res)
    {   
        if (a<0) 
        {
            std::cout << a << std::endl;
        }
    }
}
