/*
 * File to compute observalbes and check the speed of the code.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include "pybind_interface.h"
#include "pybind_interface_dwb.h"
/*
 * This function can be called if this file is linked with 
 * the .o files from the fortran libray compiled.
 */



/*
 * Function declarations
 */

/*
void check_observable(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale, std::string obs_string, bool ope_J_geq9,
   ph::constants_struct* program_const);
*/
void create_ext_pot();

void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot,bool TEST, ph::constants_struct* program_const);

/*
void check_speed(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot);

void check_interface();

void check_MWPC(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string);

bool check_chn(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, std::string chn_string, bool print_all,double** computed);

bool check_chn_all(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, bool print_all);

void check_binding_energies(std::vector<qs::quantum_channel> chns, 
        unsigned int number_of_p_points, double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot);

void check_T_matrix(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string);

bool check_observable_LO_WPC(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, bool print_all);
void check_born(std::string chn_string);
void check_DWBA(std::string chn_string);
void check_NPOT(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string);

void test_f()
{

    double V_arr[6];
    double qo = 1;
    double qi = 1;
    int S = 1;
    int J = 0;
    int T = 0;
    int Tz = 0; // np

    
    int coup = 0;
    nijmegen_fort_interface(&qi, &qo, &coup, &S, &J, &T, &Tz, &V_arr[0]); 
     
    std::cout << "my_f" << std::endl;
    
    for (int i = 0; i < 6; i++)
    {
        std::cout << V_arr[i] << " "; 
    } 
}
*/
bool TEST = false; // This will regulate the printout in the testing
double TOL = 0.03;



int main(int argc, char** argv)
{
    std::cout << "\033[1;31mbold red text\033[0m\n";
    std::cout << "\033[0;31mred text\033[0m\n";
    std::cout << "\033[1;32mbold green text\033[0m\n";
    using namespace sc;
    // test_f();
    // int a = 0;
    // std::cin >> a;
    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 76; // Number of points in angular integration
    unsigned int number_of_p_points = 120; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Define the program constants
    // ****************************
    ph::constants_struct* program_const = new ph::constants_struct;

    program_const->fpi            = 92.1;
    program_const->mpi            = 138.039; // Average of +,-,0 (PDG 2022-03) 
    program_const->Mp             = 938.2720880259; // NIST
    program_const->Mn             = 939.5654203856; // NIST
    program_const->inv_fm_to_MeV  = 197.3269804; // NIST
    program_const->MeVm2_to_mbarn = program_const->inv_fm_to_MeV
        *program_const->inv_fm_to_MeV*10.0;
    // ****************************


    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    // When checking agains Andreas code Jmax=8
    int J_max = 8;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;
    bool OPE_inclue = true;
    
 
    if (std::string(argv[2]) == "nijm")
    {
        TEST = false;
        if (std::string(argv[1]) == "phase") {
            J_max = 5;
        }
    } else {
        TEST = false;
    }
   
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    
    // Computing observables
    // compute_observables(chns,number_of_p_points,ang_int_points,J_max_in_pot,scale,Lambda,C1S0,C3S1);
    if (std::string(argv[1]) == "phase") {
        //check_phase_shifts(chns, number_of_p_points,scale, ang_int_points, 
        //        J_max_in_pot,TEST,program_const);
    }
    /*
    // Check observables
    if (std::string(argv[1]) == "DSG") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"I 0000", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "PB") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"P n000", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "CKK") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"A 00mm", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "AYY") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"C nn00", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "DT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"K 0nn0", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "PT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"P 0n00", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "AZZ") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"A 00kk", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "SPEED") {
        check_speed(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot);
    } else if (std::string(argv[1]) == "INT") {
        check_interface();
    } else if (std::string(argv[1]) == "SGT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGT", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "SGTL") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGTL", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "SGTT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGTT", OPE_inclue, program_const);
    } else if (std::string(argv[1]) == "WPC_p") {
        double* all_data[10]; 
        check_chn(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, std::string(argv[2]),true,&all_data[0]);
    } else if (std::string(argv[1]) == "WPC_p_all") {
        check_chn_all(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, false);
    } else if (std::string(argv[1]) == "WPC_PB") {
        check_observable_LO_WPC(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, false);
    } else if (std::string(argv[1]) == "DIAG") {
        check_binding_energies(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot);
    } else if (std::string(argv[1]) == "MWPC") {
        check_MWPC(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, std::string(argv[2]));
    } else if (std::string(argv[1]) == "T") {
        check_T_matrix(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, std::string(argv[2]));
    } else if (std::string(argv[1]) == "BA") {
        check_born(std::string(argv[2]));
    } else if (std::string(argv[1]) == "DWBA") {
        check_DWBA(std::string(argv[2]));
    } else if (std::string(argv[1]) == "NPOT") {
        check_NPOT(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot, std::string(argv[2]));
    }
    */
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}

/*
void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot,bool TEST, ph::constants_struct* program_const)
{
    std::cout << "Testing phase shifts with the nijmegen1 potential" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl << std::endl;
    //std::cout << "This test is now with a finite mesh" << std::endl;
    double* p_grid;
    double* w_grid;
    double Lambda = 5000.0;
    bool FINITE_GRID = false;
    double finite_grid_max = 10000.0;
    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,
                finite_grid_max,&p_grid,&w_grid);
    } else
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }

    std::vector<std::string> terms;
    terms.push_back("W_T_1pi_nu_0");
    terms.push_back("C1S0");
    terms.push_back("C3S1");
  
    Potential_mwpc<gsl_matrix> Pot = Potential_mwpc<gsl_matrix>(
            terms,ang_int_points,p_grid,w_grid,number_of_p_points,
            J_max_in_pot,450.0, 6);
    
    
    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    Pot.params_["gA"]  = 1.275; // Set correct LEC
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;

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
        std::string filename = "../../data/out_" + quantum_channel_to_string(chn) + ".txt"; 
        myfile.open(filename);
        //std::cout << "File_name: " << filename << " | ";
        

        // Read in the correct file of data
        std::string data = "../../data/np_" + quantum_channel_to_string(chn) + "_nijm1.txt";   
        std::cout << "Channel: " << quantum_channel_to_string(chn) << " | ";
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "Load: OK | ";
        } else
        {
            std::cout << "Load: Failed | ";
        }
        if (!TEST) {
            std::cout << std::endl;
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

        if (TEST)
        {
            std::cout << "E_lab d-uncoup \t dm \t dp \t eps " << std::endl; 
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
         
            //gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
            //ph::print_m(pot_V_mtx);
            gsl_matrix* pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn, false);
         
            Phase_shifts_chn phases = solver.solve_in_chn_R(T_lab,chn,pot_V_mtx);
            //Phase_shifts_chn phases = solver.solve_in_chn_T(T_lab,chn,pot_V_mtx);
            
            gsl_matrix_free(pot_V_mtx);
         
            //std::cout << T_lab << "   " << phases.delta_uncoupled*180.0/M_PI << 
            //    "   " << phases.delta_m*180.0/M_PI << "   " << phases.delta_p*180.0/M_PI<< "   " << 
            //    phases.epsilon*180.0/M_PI << std::endl;
            
            //myfile << T_lab << "   " << phases.delta_uncoupled*180.0/M_PI << 
            //    "   " << phases.delta_m*180.0/M_PI << "   " << phases.delta_p*180.0/M_PI << "   " << 
            //    phases.epsilon*180.0/M_PI << "\n"; 
            if (!chn.coupled)
            {
                //double err =std::abs((phases.delta_uncoupled*180.0/M_PI - D_delta_uncoupled[E-1])/D_delta_uncoupled[E-1]);
                double err =std::abs(phases.delta_uncoupled*180.0/M_PI - D_delta_uncoupled[E-1]);
                if (TEST) {
                    std::cout << T_lab << "   " << err << std::endl;
                }
                error = std::max(error,err);
            } else 
            {
                //if (phases.delta_m < 0) {
                //    phases.delta_m = phases.delta_m + M_PI;
                //}

                //double em = std::abs((phases.delta_m*180.0/M_PI - D_delta_m[E-1])/D_delta_m[E-1]);         
                //double ep = std::abs((phases.delta_p*180.0/M_PI   - D_delta_p[E-1])/D_delta_p[E-1]);         
                //double eps = std::abs((phases.epsilon*180.0/M_PI  - D_eps[E-1])/D_eps[E-1]);         
                if (phases.delta_m < 0.0) {
                    phases.delta_m = phases.delta_m + M_PI;
                    phases.epsilon = -phases.epsilon;
                }
                double em = std::abs(phases.delta_m*180.0/M_PI - D_delta_m[E-1]);
                double ep = std::abs(phases.delta_p*180.0/M_PI   - D_delta_p[E-1]);         
                double eps = std::abs(phases.epsilon*180.0/M_PI  - D_eps[E-1]);         
                
                if (TEST)
                { 
                    std::cout << phases.delta_m*180/M_PI << "  " << D_delta_m[E-1] << std::endl;
                    std::cout << phases.delta_p*180/M_PI << "  " << D_delta_p[E-1] << std::endl;
                    std::cout << phases.epsilon*180/M_PI << "  " << D_eps[E-1] << std::endl;
                    std::cout << std::cos(2*phases.epsilon)*std::sin(phases.delta_m + phases.delta_p) << std::endl;
                    std::cout << T_lab << "   -   " << em << "   " << ep << "   " << eps << std::endl;
                    //std::cout << T_lab << "   " << phases.delta_m << "   " << phases.delta_p << "    " << phases.epsilon << std::endl; 
                }
                error_m = std::max(em,error_m);
                error_p = std::max(ep,error_p);
                error_eps = std::max(eps,error_eps);
            } 
        
        }

        std::cout << "Max. abs err. (deg): " << error << "   " << error_m << "   " << error_p << "   " << error_eps;
        
        bool passed = false;
        if (error < TOL && error_m < TOL && error_p < TOL && error_eps < TOL) {
            passed = true;
        }
        if (passed) {
            std::cout << " | Passed: YES" << std::endl;
        } else {
            std::cout << " | Passed: NO" << std::endl;
        }

        myfile.close();
        //double a;
        //std::cout << "One channel done" << std::endl;
        //std::cin >> a;
    }
    delete[] p_grid;
    delete[] w_grid;
}
*/
