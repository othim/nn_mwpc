/*
 * This file contains the main code for calling the nn_mwpc library. This file
 * contains tests that can be run and when compiled it needs to be linked both
 * with the libnn_mwpc.co file in src/cpp_src/ and other libraries that
 * nn_mwpc depend on.
 *
 * help: run ./test -h
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


/*
 * Function declarations
 */

void check_phase_shifts(std::vector<qs::quantum_channel> chns, 
        unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
        unsigned int J_max_in_pot, ph::constants_struct* program_const,
        bool FINITE_GRID, double Lambda, double finite_grid_max,
        bool print_all);


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

double TOL = 0.03;
bool print = true;
// Location of test data
std::string DATA_DIR = "../data/"; 

// Where test logfiles are written
std::string LOG_DIR  = "test_logfiles/";

int main(int argc, char** argv)
{
    arg_struct ar = parse_arguments(argc,argv);
    
    // Set log dir
    LOG_DIR += ar.arg_value[1];

    // Create new log-file
    new_log();

    bool print_all = false;
    if (ar.arg_value[2] == "true")
    {
        print_all = true;
    }
    /*
     * ************************************************************************
     *                             Constants
     * ************************************************************************
     */

    double       scale              = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points     = 76;    // Number of points in angular integration
    unsigned int number_of_p_points = 120;   // Number of momentum-grid points
    unsigned int J_max_in_pot       = 50;    // Maximum J that is stored for L-polynomials
    int          J_max              = 5;
    int          J_min              = 0;
    int          Tz_min             = 0;
    int          Tz_max             = 0;
    double       Lambda             = 5000.0;
    bool         FINITE_GRID        = false;
    double       finite_grid_max    = 10000.0;


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

    print_color("Hello in green","green");
    print_color("Hello in red","red");
    // Do pre-computations for the physics helper functions
    ph::physics_helpers_init();
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
            Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    
    // Make phase shift checks
    check_phase_shifts(chns, number_of_p_points,scale, ang_int_points, 
            J_max_in_pot,program_const,FINITE_GRID,Lambda,finite_grid_max,
            print_all);
    
    // Free allocated memory
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}

void check_phase_shifts(std::vector<qs::quantum_channel> chns, 
        unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
        unsigned int J_max_in_pot, ph::constants_struct* program_const,
        bool FINITE_GRID, double Lambda, double finite_grid_max, bool print_all)
{
    std::cout << "Testing phase shifts with the nijmegen1 potential" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl << std::endl;

    print_to_log("Testing phase shifts with nijmegen1 potential\n");
    print_to_log("*********************************************\n");

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

    Potential_ext<gsl_matrix> nijmegen = Potential_ext<gsl_matrix>(
            p_grid, number_of_p_points, Lambda, &nijm_correct_arg);

    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid, FINITE_GRID,
            program_const);

    double mu;
    double q_on_shell;
    
    for (int c_i = 0; c_i < (int)chns.size(); c_i++)
    {
        qs::quantum_channel chn = chns[c_i]; // 1S0 channel
        // Open a file
        std::ofstream myfile;
        std::string filename = DATA_DIR + "out_" + quantum_channel_to_string(chn) + ".txt"; 
        myfile.open(filename);
        //std::cout << "File_name: " << filename << " | ";
        

        // Read in the correct file of data
        std::string data = DATA_DIR + "np_" + quantum_channel_to_string(chn) + "_nijm1.txt";   
        std::cout << "Load path: " << data << std::endl;
        std::cout << "Channel: " << quantum_channel_to_string(chn) << " | ";
        print_to_log("Channel: " + quantum_channel_to_string(chn) + " | ");
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "Load: OK | ";
        } else
        {
            std::cout << "Load: Failed | ";
        }
        double D_delta_uncoupled[300];
        double D_delta_m[300];
        double D_delta_p[300];
        double D_eps[300];

        if (!chn.coupled)
        {
            double E, d;
            int i = 0;
            while(infile >> E >> d)
            {
                D_delta_uncoupled[i] = d;
                i++;
            }
        } else
        {
            double E, dp, dm, e;
            int i = 0;
            while(infile >> E >> dm >> dp >> e)
            {
                D_delta_m[i] = dm;
                D_delta_p[i] = dp;
                D_eps[i] = e;
                i++;
            }
        }

        
        double error = 0;
        double error_m = 0;
        double error_p = 0;
        double error_eps = 0;
        
        int E_MAX = 40;
        for (int E = 1; E < E_MAX; E++)
        {
            double T_lab = (double)E;

            solver.get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
         
            gsl_matrix* pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn, false);
         
            Phase_shifts_chn phases   = solver.solve_in_chn_R(T_lab,chn,pot_V_mtx);
            
            if (print_all) {
                Phase_shifts_chn phases_T = solver.solve_in_chn_T(T_lab,chn,pot_V_mtx);
                std::cout << std::setprecision(10);
                std::cout <<"R: " <<  phases.delta_uncoupled << ",   T:";
                std::cout << phases_T.delta_uncoupled << std::endl;
            }

            gsl_matrix_free(pot_V_mtx);
         
            if (!chn.coupled)
            {
                double err =std::abs(phases.delta_uncoupled*180.0/M_PI - D_delta_uncoupled[E-1]);
                
                if (print_all)
                {
                    std::cout << T_lab << "   " << err << std::endl;
                }
                error = std::max(error,err);
            } else 
            {
                if (phases.delta_m < 0.0 && quantum_channel_to_string(chn) =="3S-D1") {
                    phases.delta_m = phases.delta_m + M_PI;
                    phases.epsilon = -phases.epsilon;
                }
                double em = std::abs(phases.delta_m*180.0/M_PI - D_delta_m[E-1]);
                double ep = std::abs(phases.delta_p*180.0/M_PI   - D_delta_p[E-1]);         
                double eps = std::abs(phases.epsilon*180.0/M_PI  - D_eps[E-1]);         
                
                if (print_all)
                {        
                    std::cout << phases.delta_m*180/M_PI << "  " << D_delta_m[E-1] << std::endl;
                    std::cout << phases.delta_p*180/M_PI << "  " << D_delta_p[E-1] << std::endl;
                    std::cout << phases.epsilon*180/M_PI << "  " << D_eps[E-1] << std::endl;
                    std::cout << std::cos(2*phases.epsilon)*std::sin(phases.delta_m + phases.delta_p) << std::endl;
                    std::cout << T_lab << "   -   " << em << "   " << ep << "   " << eps << std::endl;
                }
                error_m = std::max(em,error_m);
                error_p = std::max(ep,error_p);
                error_eps = std::max(eps,error_eps);
            } 
        }
        bool passed = false;
        if (error < TOL && error_m < TOL && error_p < TOL && error_eps < TOL) {
            passed = true;
        }
        if (passed) {
            print_color(" SUCCESS","green");
            print_to_log(" SUCCESS\n");
        } else {
            print_color(" FAILED","red");
            print_to_log(" FAILED\n");
        }

        std::cout << "                      d_uncoup \t dm \t dp \t eps " << std::endl; 
        std::cout << "Max. abs err. (deg): " << error << "   " << error_m 
            << "   " << error_p << "   " << error_eps << std::endl;
        print_to_log("Max. abs err. (deg): " + std::to_string(error) + "   " + 
                std::to_string(error_m) + "   " + std::to_string(error_p) 
                + "   " + std::to_string(error_eps) +"\n");
        

        myfile.close();
    }

    delete[] p_grid;
    delete[] w_grid;
}


arg_struct parse_arguments(int argc, char** argv)
{

    arg_struct ar;
    ar.args = {"-t","-f","-p","-h"};
    ar.arg_value = {"phase","test.out","false","false"};

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
