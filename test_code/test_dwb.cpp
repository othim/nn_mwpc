

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
#include "pybind_interface.h"
#include "pybind_interface_dwb.h"
#include <chrono>

/*
 * Function declarations
 */

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

    bool print = false;
    if (ar.arg_value[2] == "true")
    {
        print = true;
    }
    /*
     * ************************************************************************
     *                             Constants
     * ************************************************************************
     */

    double scale = 100.0;
    unsigned int number_of_p_points = 10;   // Number of momentum-grid points
    int          J_max              = 15;
 
    double       Lambda             = 500.0;
    int          cut_pow            = 6;
    bool         sharp_cutoff       = false;
    double       sharp_cutoff_add   = 300.0;
    bool         rel_corr           = false;
    bool         cut_on_shell       = true;
    bool         finite_grid        = true;
    double       finite_grid_add    = 1500.0;


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
    
    nn_mwpc_dwb_interface nn_dwb = nn_mwpc_dwb_interface(scale,J_max,Lambda,
            cut_pow,sharp_cutoff,sharp_cutoff_add, rel_corr, number_of_p_points, finite_grid,
            finite_grid_add, cut_on_shell, true, program_const->fpi,
            program_const->mpi,program_const->Mp,program_const->Mn,
            program_const->inv_fm_to_MeV);


    nn_dwb.create_new_potential("MLO","MWPC_LO_SP","DR",700.0);
    nn_dwb.print_potential_info("MLO");
    std::vector<double> params = {1.29};
    nn_dwb.set_params_in_potential("MLO",params);
    
    nn_dwb.create_new_potential("MNLO","MWPC_NLO_SP","DR",700.0);
    nn_dwb.print_potential_info("MNLO");
    params = {1.29};
    nn_dwb.set_params_in_potential("MNLO",params);
    
    nn_dwb.create_new_potential("MN2LO","MWPC_N2LO_SP","DR",700.0);
    nn_dwb.print_potential_info("MN2LO");
    params = {1.29};
    nn_dwb.set_params_in_potential("MN2LO",params);
    
    nn_dwb.create_new_potential("MN3LO","MWPC_N3LO_SP","DR",700.0);
    nn_dwb.print_potential_info("MN3LO");
    params = {1.29};
    nn_dwb.set_params_in_potential("MN3LO",params);


    //std::chrono::time_point<std::chrono::system_clock> start_time, end_time, time;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    nn_dwb.save_potential_decomposition("MLO");
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << "Time to save LO potential: " << 
        time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
    
    start_time = std::chrono::high_resolution_clock::now();
    nn_dwb.save_potential_decomposition("MNLO");
    end_time = std::chrono::high_resolution_clock::now();
    time = end_time - start_time;
    std::cout << "Time to save NLO potential: " << 
        time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
    
    start_time = std::chrono::high_resolution_clock::now();
    nn_dwb.save_potential_decomposition("MN2LO");
    end_time = std::chrono::high_resolution_clock::now();
    time = end_time - start_time;
    std::cout << "Time to save N2LO potential: " << 
        time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
    
    start_time = std::chrono::high_resolution_clock::now();
    nn_dwb.save_potential_decomposition("MN3LO");
    end_time = std::chrono::high_resolution_clock::now();
    time = end_time - start_time;
    std::cout << "Time to save N3LO potential: " << 
        time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
    
    
    //nn_dwb.save_potential_decomposition("MLO");
    std::vector<double> LECs = {-0.1/100.0,-0.13/100,1e-8,2e-8};
    nn_dwb.set_LECs_in_potential("MLO",LECs);
    nn_dwb.print_potential_info("MLO");
    
    double T_lab    = 100.0;
    double theta_cm = 80.0;
    int order = 0;

    std::vector<int> chn_index_LO = {0,1,2,3,4,7};
    std::vector<int> orders       = {0,1,2,3};
    
    for (int i=0; i<4; i++)
    {
        start_time = std::chrono::high_resolution_clock::now();
        nn_dwb.solve_save_T_chn_PC(T_lab, chn_index_LO, orders, "MLO", "MNLO","MN2LO","MN3LO");
        end_time = std::chrono::high_resolution_clock::now();
        time = end_time - start_time;
        std::cout << "Time to compute all T's: " << 
            time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
    }


    for (auto order : orders)
    {
        start_time = std::chrono::high_resolution_clock::now();
        nn_dwb.set_saved_T_vec_from_saved_orders(order);
        std::complex<double> obs_comp = nn_dwb.observable_from_saved_T_vec("I 0000",theta_cm);
        end_time = std::chrono::high_resolution_clock::now();
        time = end_time - start_time;
        std::cout <<"order: " << order <<  ". Time to compute obs    : " << 
            time/(std::chrono::microseconds(1)*1000.0) << " ms" << std::endl;
        std::cout << "obs: " << obs_comp << std::endl;
    }


    std::vector<double> T_lab_arr = {1.0,10.0,100.0};
    int order_decomp = 3;
    nn_dwb.save_decomp_T_chn_PC(T_lab_arr, chn_index_LO, order_decomp, "MLO", "MNLO","MN2LO","MN3LO");
    nn_dwb.save_decomp_T_chn_PC(T_lab_arr, chn_index_LO, order_decomp, "MLO", "MNLO","MN2LO","MN3LO");
    
    

    // Free allocated memory
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}


arg_struct parse_arguments(int argc, char** argv)
{

    arg_struct ar;
    ar.args = {"-t","-f","-p","-h"};
    ar.arg_value = {"phase","test_LO_WPC.out","false","false"};

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
    myfile << "New log:" << LOG_DIR << std::endl;
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
