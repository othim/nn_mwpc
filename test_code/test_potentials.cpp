/*
 *
 *
 *
 * This file contains code that computes phase shifts and the PB observable 
 * in WPC.
 *
 *
 * Oliver Thim 2024
 *
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <vector>
#include "physics_helpers.h"
#include "pre_def_potentials.h"
#include "potential_mwpc.h"
#include <ctime>
/*
 * Function declarations
 */
void check_potentials(int ang_int_points_, double* p_grid_, 
        double* w_grid_, int number_of_p_points_, int J_max_in_pot_,
        double cutoff_, int cut_pow_, bool sharp_cutoff_, 
        bool sharp_cutoff_add_, bool inc_weights_in_pot_, bool cut_on_shell_,
        ph::constants_struct* program_const_, bool rel_corr);
void new_log();
void print_to_log(std::string data);

void print_color(std::string data, std::string color);

struct arg_struct
{
    std::vector<std::string> args;
    std::vector<std::string> arg_value;
};
arg_struct parse_arguments(int argc, char** argv);


/*
 * Constants
 */

std::string DATA_DIR = "../data/"; 

// Where test logfiles are written
std::string OUT_DIR  = "./out/";
std::string LOG_DIR  = "./test_logfiles/";

int main(int argc, char** argv)
{
    arg_struct ar = parse_arguments(argc,argv);
    
    // Set log dir
    LOG_DIR += ar.arg_value[1];

    // Create new log-file
    new_log();

    /*
     * ************************************************************************
     *                             Constants
     * ************************************************************************
     */

    unsigned int ang_int_points     = 76;    // Number of points in angular integration
    unsigned int number_of_p_points = 100;   // Number of momentum-grid points
    unsigned int J_max_in_pot       = 50;    // Maximum J that is stored for L-polynomials
 
    double       Lambda             = 500.0;
    int          cut_pow            = 6;
    bool         sharp_cutoff       = false;
    bool         rel_corr           = true;
    bool         inc_weights_in_pot = false;
    bool         cut_on_shell       = true;
    double sharp_cutoff_add = 0;


    // Define the program constants (All in MeV)
    // ****************************
    ph::constants_struct* program_const = new ph::constants_struct;


    program_const->fpi            = 92.4;
    program_const->mpi            = 138.039;        // Average of +,-,0 (PDG 2022-03) 
    program_const->Mp             = 938.2720; 
    program_const->Mn             = 939.5653; 
    program_const->inv_fm_to_MeV  = 197.326971941683;    
    program_const->MeVm2_to_mbarn = program_const->inv_fm_to_MeV
        *program_const->inv_fm_to_MeV*10.0;

    // ****************************

    // Do pre-computations for the physics helper functions
    ph::physics_helpers_init();
    

    // Do stuff
    check_potentials(ang_int_points, nullptr, 
        nullptr, number_of_p_points, J_max_in_pot,
        Lambda, cut_pow, sharp_cutoff, 
        sharp_cutoff_add, inc_weights_in_pot, cut_on_shell,
        program_const, rel_corr);

    // Free allocated memory
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}

void check_potentials(int ang_int_points_, double* p_grid_, 
    double* w_grid_, int number_of_p_points_, int J_max_in_pot_,
    double cutoff_, int cut_pow_, bool sharp_cutoff_, 
    bool sharp_cutoff_add_, bool inc_weights_in_pot_, bool cut_on_shell_,
    ph::constants_struct* program_const_, bool rel_corr)
{

    /*
     *
     * ************************************************************************
     * Minimal example
     * ************************************************************************
     *
     */

    // Create a potential object
    Potential_mwpc<gsl_matrix_complex>* pot = pre_def_pot::create_pre_def_pot(
            "WPC_LO", ang_int_points_, p_grid_, 
            w_grid_, number_of_p_points_, J_max_in_pot_,
            cutoff_, cut_pow_, sharp_cutoff_, 
            sharp_cutoff_add_, cut_on_shell_,
            program_const_);
    
    pot->print_LECs_and_params_info();
    
    double V_arr[6]; // To hold the potential matrix elements
    
    // Define a quantum channel (3S-D1)
    qs::quantum_channel chn_3S1 = {chn_3S1.J = 1, chn_3S1.S=1,chn_3S1.T=0,chn_3S1.Tz=0,
            chn_3S1.coupled=true};
    qs::quantum_channel chn_1S0 = {chn_1S0.J = 0, chn_1S0.S=0,chn_1S0.T=1,
        chn_1S0.Tz=0, chn_1S0.coupled=false};
    
    // This variable sets if the rel corrections and cutoff is included
    // in the calc_element_V_arr_full. If false, only the raw potential
    // is returned.
    bool inc_reg_cut_and_rel = true;
    double qi = 10.0;
    double qo = 10.0;

    
    // Set params and LECs, this info can be collected from the 
    // print statement above.
    pot->params_["gA"] = 1.29;
    pot->LECs_["C1S0"] = -0.1/100.0;
    pot->LECs_["C3S1"] = -0.13/100.0;


    // Compute the potential elements
    pot->calc_element_V_arr_full(qi,qo,chn_3S1,rel_corr,inc_reg_cut_and_rel,&V_arr[0]);
    std::cout << "V_arr = [";
    for (int i=0; i<6; i++)
    {
        std::cout << V_arr[i] << " ";
    }
    std::cout << "]" << std::endl;

    
    /*
     *
     * ************************************************************************
     * Check agains the N2LO sim potential to check the potential expressions.
     * ************************************************************************
     *
     */
    
    // Create a potential object
    Potential_mwpc<gsl_matrix_complex>* pot2 = pre_def_pot::create_pre_def_pot(
            "WPC_N2LO_SFR_NO_REL", ang_int_points_, p_grid_, 
            w_grid_, number_of_p_points_, J_max_in_pot_,
            cutoff_, cut_pow_, sharp_cutoff_, 
            sharp_cutoff_add_, cut_on_shell_,
            program_const_);

    pot2->print_LECs_and_params_info();
    double C1S0 = -0.0014735727027159;
    double C3S1 = -0.0016868234304907;
    
    double D1P1 = 0.0000000051968115;
    double D1S0 = 0.0000000254759102;
    double D3P0 = 0.0000000116063756;
    double D3P1 = -0.0000000095529083;
    double D3P2 = -0.0000000065750021;
    double D3S1 = 0.0000000070477598;
    double D_SD = -0.0000000059727742;
    
    double c1 = 0.0002233095898695;
    double c3 = -0.0035553709664838;
    double c4 = 0.0039329551987581;

    pot2->params_["gA"] = 1.29;
    pot2->LECs_["C1S0"] = C1S0;
    pot2->LECs_["C3S1"] = C3S1;
    
    pot2->LECs_["D1P1"] = D1P1;
    pot2->LECs_["D1S0"] = D1S0;
    pot2->LECs_["D3P0"] = D3P0;
    pot2->LECs_["D3P1"] = D3P1;
    pot2->LECs_["D3P2"] = D3P2;
    pot2->LECs_["D3S1"] = D3S1;
    pot2->LECs_["D_SD"] = D_SD;
    
    pot2->LECs_["c1"] = c1;
    pot2->LECs_["c3"] = c3;
    pot2->LECs_["c4"] = c4;
    
    // Compute the potential elements
    pot2->print_LECs_and_params_info();
    
    inc_reg_cut_and_rel = true;
    qi = 100.0;
    qo = 100.0;
    
    pot2->calc_element_V_arr_full(qi,qo,chn_1S0,rel_corr,inc_reg_cut_and_rel,&V_arr[0]);
    std::cout << "V_arr (1S0) = [";
    std::cout << std::setprecision(10);
    for (int i=0; i<6; i++)
    {
        std::cout << V_arr[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Diff: " << std::endl;

    double diff_1S0 = std::abs(V_arr[0]+4.4866673483505186e-6);
    std::cout << diff_1S0 << std::endl;
    
    pot2->calc_element_V_arr_full(qi,qo,chn_3S1,rel_corr,inc_reg_cut_and_rel,&V_arr[0]);
    std::cout << "V_arr (3S1) = [";
    std::cout << std::setprecision(10);
    for (int i=0; i<6; i++)
    {
        std::cout << V_arr[i] << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Diff: " << std::endl;
    double diff_1 = std::abs(V_arr[2]-1.9469495194399775e-7);
    std::cout << diff_1 << std::endl;
    double diff_2 = std::abs(V_arr[3]+4.9702002132140208e-6); 
    std::cout << diff_2 << std::endl;
    double diff_3 = std::abs(V_arr[4]+1.1177779071221942e-6);
    std::cout << diff_3 << std::endl;
    double diff_4 = std::abs(V_arr[5]+1.1177779071221942e-6); 
    std::cout << diff_4 << std::endl;
    
    double tol = 3e-13;
    bool SUCCESS = false;
    if (diff_1S0 < tol && diff_1 < tol && diff_2 < tol && diff_3 < tol &&
            diff_4 < tol) 
    {
        SUCCESS = true;
    }

    if (SUCCESS)
    {
        std::cout << "TEST: SUCCESS" << std::endl;
        print_to_log("TEST AGAINST N2LO SIM POTENTIAL 1S0 and 3S1: SUCCESS");
    } else
    {
        std::cout << "TEST: FAILED" << std::endl;
        print_to_log("TEST AGAINST N2LO SIM POTENTIAL 1S0 and 3S1: FAILED");
        
    }

    

    //std::cout << "\n\nComparison to chp_preset_nnlo_sim_500_290 without CIB in OPE:" << std::endl;

    //std::cout << "qo,qi=10.0" << std::endl;
    //std::cout << "  -6.5304563871181584E-006   0.0        0.0        0.0        0.0        0.0\n " << 
    //"0.0        0.0        8.7829995132189200E-011  -6.6500620808995674E-006  -3.3882139099931637E-008  -3.3882139099931637E-008" << std::endl;
    
    
    //std::cout << "qo,qi=100.0" << std::endl;
    //std::cout <<   " -4.4866673483505186E-006   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000\n" << 
   //"0.0000000000000000        0.0000000000000000        1.9469495194399775E-007  -4.9702002132140208E-006  -1.1177779071221942E-006  -1.1177779071221942E-006";
}

arg_struct parse_arguments(int argc, char** argv)
{

    arg_struct ar;
    ar.args = {"-t","-f","-p","-h"};
    ar.arg_value = {"phase","test_potentials.out","false","false"};

    bool help = false;
    for (int i=0; i<ar.args.size(); i++)
    {
        for (int j=0; j<argc; j++)
        {
            std::string s(argv[j]);
            if (s == ar.args[i])
            {
                if (s == "-h")
                {
                    help = true;
                } else {
                    ar.arg_value[i] = argv[j+1];
                }
            }
        }
    }
    
    if (help)
    {
        std::cout << "Arguments:" << std::endl;
        for (int i=0; i<ar.args.size(); i++)
        {
            std::cout << ar.args[i] << " : " << ar.arg_value[i] << std::endl;
        }
        std::exit(0);
    }
    std::cout << std::endl;
    return ar;
}

void new_log()
{
    std::ofstream myfile;
    myfile.open(LOG_DIR);
    myfile << "New log: " << LOG_DIR << std::endl;
    std::time_t result = std::time(nullptr);
    myfile << std::asctime(std::localtime(&result)) << std::endl;
    myfile.close();
}

void print_to_log(std::string data)
{
    std::ofstream myfile;
    myfile.open(LOG_DIR, std::ios_base::app);
    myfile << data;
    myfile.close();
}

void print_color(std::string data, std::string color)
{
    if (color == "green")
    {
        std::cout << "\033[0;32m"+data+"\033[0m\n";
    } else if (color == "red")
    {
        std::cout << "\033[0;31m"+data+"\033[0m\n";
    }
}
