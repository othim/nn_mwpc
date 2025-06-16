/*
 * 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <vector>
#include "pybind_interface.h"
#include "pybind_interface_dwb.h"

int main(int argc, char** argv)
{
    bool print = false;
    /*
     * ************************************************************************
     *                             Constants
     * ************************************************************************
     */

    double       scale              = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points     = 76;    // Number of points in angular integration
    unsigned int number_of_p_points = 80;   // Number of momentum-grid points
    unsigned int J_max_in_pot       = 50;    // Maximum J that is stored for L-polynomials
    int          J_max              = 20;
    int          J_min              = 0;
    int          Tz_min             = 0;
    int          Tz_max             = 0;
    double       Lambda             = 450.0;
    bool         FINITE_GRID        = true;
    double       finite_grid_max    = 2500.0;


    // Define the program constants (All in MeV)
    // ****************************
    ph::constants_struct* program_const = new ph::constants_struct;

    program_const->fpi            = 92.4;
    program_const->mpi            = 139.5675;        // Average of +,-,0 (PDG 2022-03) 
    program_const->Mp             = 938.27231; 
    program_const->Mn             = 939.56563; 
    program_const->inv_fm_to_MeV  = 197.326971941683;    
    program_const->MeVm2_to_mbarn = program_const->inv_fm_to_MeV
        *program_const->inv_fm_to_MeV*10.0;

    // ****************************
    // Do pre-computations for the physics helper functions
    ph::physics_helpers_init();
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
            Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    
    std::cout << "Creating grid and potential..." << std::endl;
    double* p_grid;
    double* w_grid;

    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,
                finite_grid_max,&p_grid,&w_grid);
    } else
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }
    
    /*
     *
     * Define some potentials
     *
     */

    Potential_ext<gsl_matrix> idaho_n3lo = Potential_ext<gsl_matrix>(
            p_grid, number_of_p_points, 50000.0, &idaho_n3lo_correct_arg);
    
    // Make a new potential of this type
    std::vector<std::string> terms;
    terms.push_back("W_T_1pi_nu_0");
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    
    bool inc_weights_in_pot = false; 
    std::string loop_reg    = "DR";
    double lam_SFR = 700.0;
    
    std::cout << "Creating WPC_LO potential" << std::endl;
    // Make the potential complex
    Potential_mwpc<gsl_matrix_complex>* pot= 
            new Potential_mwpc<gsl_matrix_complex>(terms,96,p_grid,
            w_grid, number_of_p_points,J_max_in_pot,
            Lambda, 6, false, 1500.0, inc_weights_in_pot, 
            true,loop_reg,lam_SFR,program_const);

    double C1S0	= -0.1/100.0; 
    double C3S1	= -0.13/100.0;
    double gA   = 1.29;
    pot->params_["gA"]      = gA;
    // Save potential
    for (auto chn : chns)
    {
        pot->populate_saved_mtx(chn,false);
    }
    pot->LECs_["C1S0"] = C1S0;
    pot->LECs_["C3S1"] = C3S1;
    
    // Save HO me
    //pot->save_ho_me("tmp/LO_WPC.txt",20,24);

    //Potential_ext<gsl_matrix> cdbonn = Potential_ext<gsl_matrix>(
    //        p_grid, number_of_p_points, Lambda, &cdbonn_correct_arg);
    
    std::cout << "Constructing HO NN quantum states..." << std::endl;
    std::vector<qs::quantum_NN_HO_state> states_HO = get_states_NN_HO(20,0,0,true);
    for (int i=0; i < number_of_p_points; i++)
    {   
        std::cout << p_grid[i] << "\t" << w_grid[i] << std::endl;
    }
    
    double mN = 2*program_const->Mp*program_const->Mn/(program_const->Mp+program_const->Mn);
    double Omega = 24;
    std::cout << "mN=" << mN << ", Omega=" << Omega << std::endl; 
    
    /*
    for (int i=0; i < number_of_p_points; i++)
    {
        for (int j=0; j < number_of_p_points; j++)
        {
            double el_1 = idaho_n3lo.get_pot_element_LSJ(p_grid[j],p_grid[i],0,0,0,0,1,1);
            double el_0 = idaho_n3lo.get_pot_element_LSJ(p_grid[j],p_grid[i],0,0,0,0,1,0);
            std::cout << (1.0/3.0)*el_0 + (2.0/3.0)*el_1 << std::endl;
        }
    }
    //return 0;
    double* R;
    R = ph::get_mom_HO_R(p_grid, number_of_p_points, 1, 1, mN, Omega);
    for (int i=0; i < number_of_p_points; i++)
    {
        std::cout << p_grid[i] << "\t" << R[i] << std::endl;
    }
    free(R);
    */

    // Print elements to file
    //std::string save_dir = "tmp";
    //pot->save_ho_me_decomp(save_dir, 4, 24);  
    
    
    std::cout << "n" << "\t" <<  "l" << "\t" << "np" << "\t" <<
        "lp" << "\t" << "S" << "\t" << "J" << "\t" << "Tz" << 
        "\t" << "ME" << std::endl;

    for (auto bra : states_HO)
    {   
        for (auto ket : states_HO)
        {
            if (bra.S == ket.S && bra.J == ket.J && bra.T == ket.T && 
                    bra.Tz==ket.Tz)
            {
                double el = idaho_n3lo.compute_HO_matrix_el(bra.n, bra.L, ket.n, 
                    ket.L, bra.S, bra.J, bra.T, bra.Tz, p_grid, w_grid, number_of_p_points, mN, Omega);
                
                //double el = pot->compute_HO_matrix_el(bra.n, bra.L, ket.n, 
                //    ket.L, bra.S, bra.J, bra.T, bra.Tz, p_grid, w_grid, number_of_p_points, mN, Omega);
                
                // n l n’ l’ S J Tz ME
                std::cout << std::setprecision(16);
                std::cout << bra.n << "\t" <<  bra.L << "\t" << ket.n << "\t" <<
                    ket.L << "\t" << bra.S << "\t" << bra.J << "\t" << bra.Tz << 
                    "\t" << el << std::endl;
            
            }
        }
    
    }
    
    // Free allocated memory
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}
