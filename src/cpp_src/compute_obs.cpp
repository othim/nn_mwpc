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
#include "pot_ext.h"
#include "born_approx.h"
#include "potential_mwpc.h"
#include "physics_helpers.h"
/*
 * This function can be called if this file is linked with 
 * the .o files from the fortran libray compiled.
 */



/*
 * Function declarations
 */

void check_observable(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale, std::string obs_string, bool ope_J_geq9);

void create_ext_pot();

void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot,bool TEST);


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

bool TEST = false; // This will regulate the printout in the testing
double TOL = 0.03;

int main(int argc, char** argv)
{

    using namespace sc;
    // test_f();
    //int a = 0;
    //std::cin >> a;
    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 76; // Number of points in angular integration
    unsigned int number_of_p_points = 120; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    
    // When checking agains Andreas code Jmax=8
    int J_max = 14;
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
        check_phase_shifts(chns, number_of_p_points,scale, ang_int_points, J_max_in_pot,TEST);
    }
    
    // Check observables
    if (std::string(argv[1]) == "DSG") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"I 0000", OPE_inclue);
    } else if (std::string(argv[1]) == "PB") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"P n000", OPE_inclue);
    } else if (std::string(argv[1]) == "CKK") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"A 00mm", OPE_inclue);
    } else if (std::string(argv[1]) == "AYY") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"C nn00", OPE_inclue);
    } else if (std::string(argv[1]) == "DT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"K 0nn0", OPE_inclue);
    } else if (std::string(argv[1]) == "PT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"P 0n00", OPE_inclue);
    } else if (std::string(argv[1]) == "AZZ") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"A 00kk", OPE_inclue);
    } else if (std::string(argv[1]) == "SPEED") {
        check_speed(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot);
    } else if (std::string(argv[1]) == "INT") {
        check_interface();
    } else if (std::string(argv[1]) == "SGT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGT", OPE_inclue);
    } else if (std::string(argv[1]) == "SGTL") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGTL", OPE_inclue);
    } else if (std::string(argv[1]) == "SGTT") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, 
                J_max_in_pot,"SGTT", OPE_inclue);
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
    ph::physics_helpers_free();
    return 0;
}

void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot,bool TEST)
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
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
  
    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0, 6);
/* 
   for (auto chn : chns)
   {
      Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
   }
*/
    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;

    Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, Lambda, &nijm_correct_arg);

    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid, FINITE_GRID);

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

            LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
         
            //gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
            //ph::print_m(pot_V_mtx);
            gsl_matrix* pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn);
         
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

void check_observable(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale,std::string obs_string, bool ope_J_geq_9)
{
    std::cout << "Testing observables with the nijmegen1 potential" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl << std::endl;
    std::cout << "This test is now with an inf mesh" << std::endl;
    //std::clock_t start, end;   
    double Lambda = 5000.0;
    // Make grid
    double* p_grid;
    double* w_grid;
    bool FINITE_GRID = false;
    
    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,8000,&p_grid,&w_grid);
    } else
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }
    

    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    // Choose terms in LO WPC potential
    
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");

    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0,6);
   
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }
  
    Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;

    // Construct potential and LS-Solver 

    // Choose terms in OPE potential
    
    int l_max = 50;
    std::vector<std::string> terms2;
    terms2.push_back("OPEP"); // To just test elements use just OPEP

    Potential_mwpc OPE = Potential_mwpc(terms2,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,Lambda,6);
  
    
    for (auto chn : chns)
    {
        if (chn.J > 9)
        {
            OPE.populate_saved_mtx(chn,true); // Realtivistic factor on
        }
    }
  
    //OPE.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
    OPE.LECs_["gA2"]  = 1.29*1.29; // Set correct LEC

    // Compute the observables
    std::string obs_string2;

    if (obs_string == "I 0000") {
        obs_string2 = "DSG";
    } else if (obs_string == "P n000") {
        obs_string2 = "PB";
    } else if (obs_string == "K 0nn0") {
        obs_string2 = "DT";
    } else if (obs_string == "P 0n00") {
        obs_string2 = "PT";
    } else if (obs_string == "A 00mm") {
        obs_string2 = "CKK";
    } else if (obs_string == "C nn00") {
        obs_string2 = "AYY";
    } else if (obs_string == "A 00kk") {
        obs_string2 = "AZZ";
    } else if (obs_string == "SGT") {
        obs_string2 = "SGT"; // Special
    } else {
        obs_string2 = obs_string;
    }

    Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, Lambda, 
            &nijm_correct_arg);
    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    double rho_T;
   
    const int len = 3;

    double energies1[len] = {10.0,50.0, 200.0};
    double energies2[len] = {125.0, 200.0,350.0};
    double* energies; 
    if (obs_string2 == "DSG")
    {
        energies = energies2;
    } else 
    {
        energies = energies1;
    }

    for (int i = 0; i < len; i++)
    {
        double Tl = energies[i];
        // Compute all the phase shifts in the channels
        std::vector<Phase_shifts_chn> phases_vec;
        //start = std::clock();
        std::string data;
        if (obs_string2 != "SGT" && obs_string2 != "SGTL" && obs_string2 != "SGTT")
        {
            std::cout << "Testing " + obs_string2 + " with T_lab=" << Tl << " MeV" << std::endl;
    
            // Read in the correct file of data
            data = "../../data/np_" + obs_string2 + "_" + std::to_string((int)Tl) + "_nijm1.txt";   
        } else
        {
            data = "../../data/np_" + obs_string2 +  "_nijm1.txt";
        }
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "File" + data + " loaded: OK" << std::endl;
        } else
        {
            std::cout << "File" + data + ": Failed" << std::endl;
        }

        // Read and save the data to arrays
        double D_obs[180];
        double theta, obs;
        int k = 0;
        while(infile >> theta >> obs)
        {
            D_obs[k] = obs;
            k++;
            //std::cout << obs << std::endl;
        }
        
        for (auto chn : chns)
        {
            LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
                      
            //gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
            
            gsl_matrix* pot_V_mtx;
            if (ope_J_geq_9)
            {
                if (chn.J < 10)
                {
                    pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn);
                } else 
                {
                    pot_V_mtx = OPE.get_saved_matrix(q_on_shell, chn, true);
                }
            } else
            {
                pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn);
            }
            //Phase_shifts_chn phases = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
            Phase_shifts_chn phases = solver.solve_in_chn_T(Tl,chn,pot_V_mtx);
            
            gsl_matrix_free(pot_V_mtx);
            phases_vec.push_back(phases);
        }
        
        std::ofstream myfile;
        std::string filename;
        
        if (obs_string2 == "SGT" || obs_string2 == "SGTT" || obs_string2 == "SGTL")
        {
            filename = "../../data/out_" + obs_string2 + ".txt"; 
        } else 
        {
            filename = "../../data/out_" + obs_string2 + "_" + std::to_string((int)Tl) + ".txt"; 
        }
        myfile.open(filename);


        double errors[180];
        double mean_error = 0;
        if (obs_string == "SGT" || obs_string == "SGTT" || obs_string == "SGTL")
        {
            // Get Saclay amplitudes
            std::vector<std::complex<double> > saclay_amplitudes;
            std::cout << "Energy \t obs \t correct \t abs. rel. error" << std::endl;    
            int e = i+1;
            // Loop over energies
            LS_Solver::get_mu_q_on_shell(energies[e-1],chns[0], &mu,&q_on_shell);
            rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);
            saclay_amplitudes = sc::compute_Saclay_amplitudes(chns, phases_vec, 0, q_on_shell, rho_T, l_max);

            // Compute the observable from the amplitudes
            double obs = sc::compute_observable(saclay_amplitudes, q_on_shell, obs_string);
            
            if (D_obs[e-1] != 0) {
                errors[e-1] = std::abs((D_obs[e-1] - obs)/D_obs[e-1]);
            } else {
                errors[e-1] = 0;
            }
            std::cout << energies[e-1] << "\t" << obs << "\t" << D_obs[e-1] <<"\t" <<  errors[e-1] << std::endl << std::endl;
            myfile << energies[e-1] << "\t" << obs << "\t" << D_obs[e-1] << "\t" << errors[e-1] << std::endl;
        } 
        else
        {

            std::cout << "Angle \t obs \t correct \t abs. rel. error" << std::endl;    
            for (int ang = 1; ang < 181; ang++)
            {
                double angle = (double)ang;
                if (ang == 90) {
                    angle = 90.001;
                }
                // Get Saclay amplitudes
                std::vector<std::complex<double> > saclay_amplitudes;

                LS_Solver::get_mu_q_on_shell(Tl,chns[0], &mu,&q_on_shell);

                rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);
                saclay_amplitudes = sc::compute_Saclay_amplitudes(chns, phases_vec, angle*M_PI/180.0, q_on_shell, rho_T, l_max);

                // Compute the observable from the amplitudes
                double obs = 0;
                if (obs_string2 == "AZZ")
                {
                    obs = sc::compute_observable_lab(saclay_amplitudes, q_on_shell, obs_string, angle*M_PI/180.0);
                } else
                {
                    obs = sc::compute_observable(saclay_amplitudes, q_on_shell, obs_string);
                }

                //std::cout << D_obs[ang-1] << " " << obs << std::endl; 
                if (D_obs[ang-1] != 0) {
                    errors[ang-1] = std::abs((D_obs[ang-1] - obs)/D_obs[ang-1]);
                } else {
                    errors[ang-1] = 0;
                }
                //std::cout << angle << " a: " << saclay_amplitudes[0] << " b: " << saclay_amplitudes[1] <<
                //    " c: " << saclay_amplitudes[2] << " d: " << saclay_amplitudes[3] << " e: " << saclay_amplitudes[4] << std::endl; 
            
                //std::cout << angle << "\t" << obs << "\t" << D_obs[ang-1] << "\t" << errors[ang-1]  << std::endl;
                myfile << angle << "\t" << obs << "\t" << D_obs[ang-1] << "\t" << errors[ang-1] << std::endl;
                mean_error += errors[ang-1];

                /*
                // With the M-matrix
                double obs_M;
                if (obs_string == "I 0000")
                {
                g   sl_matrix_complex* M_matrix = get_M_matrix(chns, phases_vec, q_on_shell, angle*M_PI/180.0, rho_T, l_max);
                gsl_matrix_complex* eig = gsl_matrix_complex_alloc(2,2);
                gsl_matrix_complex_set_identity(eig);
                ph::print_m(eig);
                ph::print_m(M_matrix);
                obs_M = get_observables(eig,eig,eig,eig,M_matrix); 
            }            
            std::cout << angle << "   " << obs_M << " mb" << std::endl; */
            }
            std::cout << "Mean absolute relative error: " << mean_error/180.0  << std::endl;
            std::cout << "Maximum error: " << *(std::max_element(errors, errors + 180)) << std::endl;
        }
        myfile.close();
        // Now all the pahse shifts in the relevent channels are known.
      // Now we can compute the total cross section for some on_shell
      // lab energy

      //LS_Solver::get_mu_q_on_shell(Tl,chns[0], &mu,&q_on_shell);
   
      //rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);

      //start = std::clock();
      //double cross_section = compute_total_cross_section(chns,phases_vec,q_on_shell,rho_T,30);
      
      // Compute defferential cross section 
      // ----------------------------------
      
      // Get the Saclay amplitudes
      
      // Get the M-matrix 
      
      // Compute observable from the M-matrix

      // ---------------------------------- 
     // end = std::clock();
      //std::cout << "Time to compute cross section: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;

   
      //std::cout << Tl << "\t \t" << cross_section << std::endl;
    } 
    //delete[] p_grid;
    //delete[] w_grid;
}

void check_speed(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot)
{
    std::cout << "Testing speed of code with LO WPC potential and the observable DSG" << std::endl << std::endl;
    int cut_pow = 6;
    std::clock_t start, end;   
    
    // Make grid
    double* p_grid;
    double* w_grid;
    bool FINITE_GRID = false;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");

    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0,cut_pow,false);
    
    std::cout << "Saving potential matrices" << std::endl;
    start = std::clock();   
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }
    end = std::clock();
    std::cout << "Time to save matrices: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 

    // Set correct LECs
    Pot.LECs_["gA2"]  = constants::gA*constants::gA;
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;

    int l_max = 50;

    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    double rho_T;
    
    double Tl = 50.0; // MeV

    // Compute all the phase shifts in the channels
    std::vector<Phase_shifts_chn> phases_vec;
    int count = 0;
    for (auto chn : chns)
    {
        std::cout << "Channel " << count << std::endl;
        count++;
        
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        
        start = std::clock();        
        Pot.LECs_["C1S0"] = C1S0;
        Pot.LECs_["C3S1"] = C3S1;
        gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);
        end = std::clock();
        std::cout << "Get V: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 

        start = std::clock();   
        Phase_shifts_chn phases = solver.solve_in_chn_T(Tl,chn,pot_V_mtx);
        end = std::clock();
        std::cout << "Solve LS (T): " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl;
        
        start = std::clock();   
        phases = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
        end = std::clock();
        std::cout << "Solve LS (R): " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl << std::endl;
        
        gsl_matrix_free(pot_V_mtx);
        phases_vec.push_back(phases);
    }
    phases_vec.clear();

    start = std::clock();
    for (auto chn : chns)
    {
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        Pot.LECs_["C1S0"] = C1S0;
        Pot.LECs_["C3S1"] = C3S1;
        gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);

        Phase_shifts_chn phases = solver.solve_in_chn_T(Tl,chn,pot_V_mtx);
        
        gsl_matrix_free(pot_V_mtx);
        phases_vec.push_back(phases);
    }
    end = std::clock();
    std::cout << "Total time to solve LS (T): " << 1e3*(double)(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
    phases_vec.clear();

    start = std::clock();
    for (auto chn : chns)
    {
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        Pot.LECs_["C1S0"] = C1S0;
        Pot.LECs_["C3S1"] = C3S1;
        gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);

        Phase_shifts_chn phases = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
        
        gsl_matrix_free(pot_V_mtx);
        phases_vec.push_back(phases);
    }
    end = std::clock();
    std::cout << "Total time to solve LS (R): " << 1e3*(double)(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl << std::endl;
    // Compute observable for angle 180 angles
    //start = std::clock();
    for (int ang = 1; ang < 181; ang++)
    {
        double angle = (double)ang;
        if (ang == 90) {
            angle = 90.001;
        }
        // Get Saclay amplitudes
        
        start = std::clock();
        std::vector<std::complex<double> > saclay_amplitudes;
        LS_Solver::get_mu_q_on_shell(Tl,chns[0], &mu,&q_on_shell);
        rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);
        saclay_amplitudes = sc::compute_Saclay_amplitudes(chns, phases_vec, angle*M_PI/180.0, q_on_shell, rho_T, l_max);
        end = std::clock();
        if (ang == 1) {
            std::cout << "Get Sac. amp.: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl;
        } 
        // Compute the observable from the amplitudes
        start = std::clock();
        double obs = sc::compute_observable(saclay_amplitudes, q_on_shell, "I 0000");
        end = std::clock();
        obs = obs+1.0; // Just to not get unused warning 
        if (ang == 1) {
            std::cout << "Compute OBS from Sac. amp.: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl << std::endl;
        }
    } 
    
    //end = std::clock();
    //std::cout << "Time calculate OBS for 180 angles: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl;

    delete[] p_grid;
    delete[] w_grid;
}

void check_interface()
{
    std::vector<double> angles = {1.0,5.0};
    std::vector<double> energies = {10.0,50.0};
    std::string observable = "I 0000";
    std::vector<double> LECs = {-0.112927/100.0,-0.087340/100.0,1.289*1.289};

    nn_mwpc_interface* obj = new nn_mwpc_interface("cdbonn",25,450.0,6,false, true, true,60,false);
    
    std::cout << "Object created" << std::endl;

    obj->solve_LS_ext_pot(10.0);
    double o = obj->compute_observable("SGT", 50.0);
    std::cout << o << std::endl;
    
    //std::vector<double> res = obj->compute_observable_l(observable,angles,energies,LECs);
   
    //obj->solve_LS(10.0,LECs);
    //double a = obj->compute_observable("I 0000", 80.0);
    //a = a + 1.0;
    delete obj;
    std::cout << "Done!" << std::endl; 
}

bool check_chn(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, std::string chn_string, bool print_all, double** computed)
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Testing  phase shifts with the WPC_LO potential." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Channel: " << chn_string << std::endl;
    std::vector<std::string> files;
    
    #if defined(ANDREAS_CONST)
    files.push_back("../../data/phase_shifts_Andreas_original_part.txt");
    #else
    files.push_back("../../data/data_gen_corr_Andreas.txt");
    #endif
    // Make constants the same
    // -----------------------
    // Set the constants to the values Andreas use
    
    double Lambda = 500.0;
    double C1S0	= -0.1/100.0; 
    double C3S1	= -0.13/100.0;
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // -----------------------
    
    int cut_pow = 6; 
    for (auto& data : files)
    {
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "File" + data + " loaded: OK" << std::endl;
        } else
        {
            std::cout << "File" + data + ": Failed" << std::endl;
        }

        // Read and save the data to arrays
        double D_E[350];
        double D_1S0[350];
        double D_3S1[350];
        double D_3D1[350];
        double D_E1[350];
        double D_1P1[350];
        double D_3P1[350];
        double D_3P0[350];
        double D_3P2[350];
        double D_3F2[350];
        double D_E2[350];
        
        const int nrow = 350;
        
        double* computed_tmp = new double[nrow];

        double E, a1S0, a3S1, aE1, a1P1, a3P1, a3P0, a3P2, aE2, a3F2, a3D1;
        int k = 0;
        //# Tlab,np 1S0,np 3S1, np E1, np 1P1, np 3P1, np 3P0, np 3P2, np E2,    
        //  np 3D1, np 3F2 

        // Read the first line, the heading
        std::string heading;
        std::getline(infile, heading);
        //std::cout << heading << std::endl;
        while(infile >> E >> a1S0 >> a3S1 >> aE1 >> a1P1 >> a3P1 >>
                a3P0 >> a3P2 >> aE2 >> a3D1 >> a3F2) 
        {
            D_E[k]   = E;
            D_1S0[k] = a1S0;
            
            D_3S1[k] = a3S1;
            D_3D1[k] = a3D1;
            D_E1[k]  = aE1;
            
            D_1P1[k] = a1P1;
            D_3P1[k] = a3P1;
            D_3P0[k] = a3P0;
            
            D_3P2[k] = a3P2;
            D_3F2[k] = a3F2;
            D_E2[k]  = aE2;

            /*std::cout << E << "   " << a1S0 << "   " << a3S1 << "   " << 
                aE1 << "   " << a1P1 << "   " << a3P1 << "   " << 
                a3P0 << "   " << a3P2 << "   " << "   " << aE2 << "   " 
                << a3D1 << "   " << a3F2 << std::endl;
            */
            k++;
            //std::cout << E << "   " << phase << std::endl;
        }
        
        // Make grid
        double* p_grid;
        double* w_grid;
        bool FINITE_GRID = false;
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

        // Choose terms in LO WPC potential
        std::vector<std::string> terms;
        terms.push_back("OPEP"); // To just test elements use just OPEP
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("C3P0");
        terms.push_back("C3P2");


        Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,
                w_grid,number_of_p_points,J_max_in_pot,Lambda,cut_pow,false);
        
        for (auto chn : chns)
        {
            Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
        }

        // Set correct LECs
        Pot.LECs_["gA2"]  = constants::gA*constants::gA;
        Pot.LECs_["C1S0"] = C1S0;
        Pot.LECs_["C3S1"] = C3S1;
        Pot.LECs_["C3P0"] = C3P0;
        Pot.LECs_["C3P2"] = C3P2;

        LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid, FINITE_GRID);
       
        double q_on_shell;
        double mu;
        
        std::string chn_string_full = chn_string;
        if (chn_string=="3S1" || chn_string=="3D1" ||chn_string=="E1")
        {
            chn_string_full = "3S-D1";
        } else if (chn_string=="3P2" || chn_string=="3F2" || chn_string=="E2")
        {
            chn_string_full = "3P-F2";
        }
        // Take just the relevant channel
        qs::quantum_channel chn;
        for (int i = 0; i < (int)chns.size(); i++)
        {
            chn = chns[i];
            if (quantum_channel_to_string(chn) == chn_string_full) {
                break;
            }
        }
        
        double max_err  = 0;
        double mean_err = 0;
        if (print_all)
        {
            std::cout << std::setw(5) << "T_lab" << "   " << std::setw(10) << 
                "data" << "   " << std::setw(10) << "C_phase" 
                << "   " << std::setw(8) << "err" << std::endl;
        }
        for (int i = 0; i < 350; i++)
        {
            LS_Solver::get_mu_q_on_shell(D_E[i], chn, &mu, &q_on_shell);
            
            // Use other convention to match the results of andreas code
            gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);
            Phase_shifts_chn phases = solver.solve_in_chn_R(D_E[i],chn,pot_V_mtx);
            
            
            double* D_phase = nullptr;
            double C_phase  = 0.0;
            if (chn_string == "1S0")
            {
                D_phase = D_1S0;
                C_phase = phases.delta_uncoupled*180.0/M_PI;
                
            } else if (chn_string == "1P1")
            {
                D_phase = D_1P1;
                C_phase = phases.delta_uncoupled*180.0/M_PI;
            } else if (chn_string == "3P0")
            {
                D_phase = D_3P0;
                C_phase = phases.delta_uncoupled*180.0/M_PI;
            } else if (chn_string == "3P1")
            {
                D_phase = D_3P1;
                C_phase = phases.delta_uncoupled*180.0/M_PI;
            } else if (chn_string == "3S1")
            {
                D_phase = D_3S1;
                C_phase = phases.delta_m*180.0/M_PI;
            } else if (chn_string == "3D1")
            {
                D_phase = D_3D1;
                C_phase = phases.delta_p*180.0/M_PI;
            } else if (chn_string == "E1")
            {
                D_phase = D_E1;
                C_phase = phases.epsilon*180.0/M_PI;
            } else if (chn_string == "3P2")
            {
                D_phase = D_3P2;
                C_phase = phases.delta_m*180.0/M_PI;
            } else if (chn_string == "3F2")
            {
                D_phase = D_3F2;
                C_phase = phases.delta_p*180.0/M_PI;
            } else if (chn_string == "E2")
            {
                D_phase = D_E2;
                C_phase = phases.epsilon*180.0/M_PI;
            }
            computed_tmp[i] = C_phase;
            double err = std::abs(D_phase[i] - C_phase);
            if (print_all)
            {
            std::cout << std::setw(5) << D_E[i] << "   " << std::setprecision(8) 
                << std::setw(10) << 
                D_phase[i] << "   " << std::setw(10) << C_phase 
                << "   " << std::setw(8) << err << std::endl;
            
            }
            max_err = std::max(max_err, err);
            mean_err += err; 

            delete[] pot_V_mtx;
        }

        *computed = computed_tmp;
        std::cout << "Mean error (deg): " << mean_err/350.0 << std::endl <<
            "Max error (deg): " << max_err << std::endl;
        double tol = 2e-4;
        if (mean_err/350.0<tol && max_err<tol)
        {
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "******** Test: OK (abs.tol=" << tol << ") ********" << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            return true;
        } else 
        {   
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "******** Test: FAILED (abs.tol=" << tol << ") ********" << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            return false;
        }

    }
    return false;
}

bool check_observable_LO_WPC(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, bool print_all)
{
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    int J_max = 20;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;
    
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    chns = get_channels(states, print);   
    
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "     Testing PB with the WPC_LO potential.      " << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::vector<std::string> files;
    
    #if defined(ANDREAS_CONST)
    files.push_back("../../data/PB_30_MeV_Andreas_original.txt");
    #else
    files.push_back("../../data/PB_30_MeV_Andreas_corr.txt");
    #endif
    // Make constants the same
    // -----------------------
    // Set the constants to the values Andreas use
    
    double Lambda = 500.0;
    double C1S0	= -0.1/100.0; 
    double C3S1	= -0.13/100.0;
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // -----------------------
    
    int cut_pow = 6; 
    for (auto& data : files)
    {
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "File" + data + " loaded: OK" << std::endl;
        } else
        {
            std::cout << "File" + data + ": Failed" << std::endl;
        }

        // Read the first line, the heading
        std::string heading;
        std::getline(infile, heading);
        //std::cout << heading << std::endl;
        double D_PB[180];
        int k=0;
        double t_cm, t_lab,PB, th;
        while(infile >> t_cm >> t_lab >> th >> PB >> th >> th >> th >> th) 
        {
            D_PB[k]   = PB;
            //std::cout << t_cm << "   " << PB << std::endl;
            k++;
        }
        
        // Make grid
        double* p_grid;
        double* w_grid;
        bool FINITE_GRID = false;
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

        // Choose terms in LO WPC potential
        std::vector<std::string> terms;
        terms.push_back("OPEP"); // To just test elements use just OPEP
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("C3P0");
        terms.push_back("C3P2");


        Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,
                w_grid,number_of_p_points,J_max_in_pot,Lambda,cut_pow,false);
        
        for (auto chn : chns)
        {
            Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
        }
        
        // Set correct LECs
        Pot.LECs_["gA2"]  = constants::gA*constants::gA;
        Pot.LECs_["C1S0"] = C1S0;
        Pot.LECs_["C3S1"] = C3S1;
        Pot.LECs_["C3P0"] = C3P0;
        Pot.LECs_["C3P2"] = C3P2;

        LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,
                FINITE_GRID);
       
        
        
        double Tl = 30.0; // MeV
        std::string obs_string = "P n000";
        double q_on_shell;
        double mu;
        std::vector <Phase_shifts_chn> phases_vec;
        for (auto chn : chns)
        {
            LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
            gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
            Phase_shifts_chn phases = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
            phases_vec.push_back(phases);
            gsl_matrix_free(pot_V_mtx);
        }

        double rho_T = M_PI*q_on_shell*mu; // In the convention used

        double max_err  = 0;
        double mean_err = 0;
        for (int i = 0; i < 179; i++)
        {
            double angle = (double)(i+1);
            std::vector<std::complex<double> > saclay_amplitudes;
            // Convert the angle to radians. This uses the pre-computed phase
            // shifts phase_shifts_ that are stored in the class.
            if (std::abs(angle-90.0)<1e-4)
            {
                angle = 90.001;
            }
            saclay_amplitudes = sc::compute_Saclay_amplitudes(chns, phases_vec, 
                angle*M_PI/180.0, q_on_shell, rho_T, J_max_in_pot);
            
            double obs = sc::compute_observable(saclay_amplitudes, q_on_shell, obs_string);
            
            double err = std::abs(D_PB[i] - obs)/std::abs(obs);
            if (print_all)
            {
                std::cout << std::setw(5) << angle << "   "<< std::setw(10) << D_PB[i] << "   " << std::setprecision(8) 
                    << std::setw(10) <<  obs << "   " << std::setw(10) << err << std::endl;
            }
            max_err = std::max(max_err, err);
            mean_err += err; 

        }
        std::cout << "Mean rel. error (%): " << 100.0*mean_err/179.0 << std::endl <<
            "Max rel. error (%): " << 100.0*max_err << std::endl;
        double tol = 2e-4;
        if (mean_err/350.0<tol && max_err<tol)
        {
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "******** Test: OK (abs.tol=" << tol << ") ********" << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            return true;
        } else 
        {   
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "******** Test: FAILED (abs.tol=" << tol << ") ********" << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            return false;
        }
    }
    return false;
}

bool check_chn_all(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot, bool print_all)
{
    std::cout << std::endl << "------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "******** STARTING TEST OF PHASE SHIFTS *********" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    //# Tlab,np 1S0,np 3S1, np E1, np 1P1, np 3P1, np 3P0, np 3P2, np E2,    
    //  np 3D1, np 3F2 
    std::vector<std::string> chn_strings;
    chn_strings.push_back("1S0");
    chn_strings.push_back("3S1");
    chn_strings.push_back("E1");
    chn_strings.push_back("1P1");
    chn_strings.push_back("3P1");
    chn_strings.push_back("3P0");
    chn_strings.push_back("3P2");
    chn_strings.push_back("E2");
    chn_strings.push_back("3D1");
    chn_strings.push_back("3F2");

    bool success = true;

    int nrows = 350;
    int ncol  = 10;
    
    double* all_data[ncol]; 
    for(int i=0; i<(int)chn_strings.size();i++)
    {
        std::string chn_string = chn_strings[i];
        bool suc = check_chn(chns, number_of_p_points, scale, ang_int_points, 
            J_max_in_pot, chn_string, print_all, &all_data[i]);
        if (suc==false) {success=false;}
        //std::cout << std::setprecision(8) << std::setw(10) << 
        //    all_data[0][349] << "   ";   
    }

    std::cout << std::endl << "------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "****************** TEST DONE *******************" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    if (success)
    {
        std::cout << "Tests OK in chns: ";
    } else 
    {
        std::cout << "Tests FAILED in chns: ";
    }
    for(int j=0; j<(int)chn_strings.size(); j++)
    {
        std::cout << chn_strings[j] << " ";
    }
    std::cout << std::endl;
    #ifdef ANDREAS_CONST
    std::cout << "This test was with the compile flag FLAGS=-DANDREAS_CONST \n"
        "which sets the C++ variable ANDREAS_CONST as defined as that \n"
        "selects the correct constants" << std::endl;
    #else
    std::cout << "These tests was wich my constants and the kinematics \n"
    "relation that is reversed w.r.t. Andreas code. For this i compare \n"
    "to files where I have generated the data with my code when it was \n"
    "verified." << std::endl;
    #endif
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    if (print_all)
    {
        for(int i=0; i<nrows; i++)
        {
            std::cout << std::setprecision(4) << std::setw(5) << i+1 << "   ";   
            for(int j=0; j<ncol; j++)
            {
                std::cout << std::setprecision(6) << std::setw(10) << 
                    (all_data[j])[i] << "   ";   
            }
            std::cout << std::endl;
        }
    }
    for(int i=0; i<ncol; i++)
    {
        delete[] all_data[i];
    }

    return success;
}

void check_binding_energies(std::vector<qs::quantum_channel> chns, 
        unsigned int number_of_p_points, double scale,unsigned int ang_int_points, 
        unsigned int J_max_in_pot)
{
    std::cout << "Testing binding energies with Nijmegen1 potential." << std::endl;
    std::cout << "-------------------------------------------------" << std::endl << std::endl;
    int cut_pow = 6;

    // Make grid
    double* p_grid;
    double* w_grid;
    double Lambda = 450.0;
    bool FINITE_GRID = false;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

    double C1S0	= -0.1/100.0; 
    //double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    double C3S1 = -0.081/100.0;
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    terms.push_back("C3P0");
    terms.push_back("C3P2");


    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,Lambda, cut_pow,false);
    
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }

    // Set correct LECs
    Pot.LECs_["gA2"]  = 1.29*1.29;
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;
    Pot.LECs_["C3P0"] = C3P0;
    Pot.LECs_["C3P2"] = C3P2;

    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
    
    qs::quantum_channel chn = chns[3]; // 3S1-3D1
    // WPC potential
    // ------------- 
    gsl_matrix* pot_V_mtx = Pot.get_matrix_no_onshell(chn, true);
    ph::eigen_t diag_res = ph::solve_SE(p_grid, w_grid, number_of_p_points, chn, pot_V_mtx);
    if (TEST) {
        std::cout << "The eigenvalues in 3S1-3D1 LO MWPC" << quantum_channel_to_string(chn) << " is: "
            << std::endl;
        ph::print_v(diag_res.eigenvalues);    
    }
    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_complex_free(diag_res.eigenvalues); 
    delete[] pot_V_mtx;
    // ------------
    // Nijmegen1  E = -2.224575 (nn-on-line)
    // ------------
    double E_nijmegen = -2.224575;
    Lambda = 10000.0;
    double mu,q_on_shell;
    LS_Solver::get_mu_q_on_shell(0.0, chn, &mu, &q_on_shell);
    
    //Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, Lambda, &nijm_correct_arg);
    Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, Lambda, &cdbonn_correct_arg);
    pot_V_mtx = nijmegen.get_matrix(10.0,chn);
    pot_V_mtx = nijmegen.get_matrix_no_onshell(chn);
    diag_res = ph::solve_SE(p_grid, w_grid, number_of_p_points, chn, pot_V_mtx);
    if (TEST) {
        std::cout << "The eigenvalues in 3S1-3D1 LO nijmegen1" << quantum_channel_to_string(chn) << " is: "
            << std::endl;
        ph::print_v(diag_res.eigenvalues);    
    }
    std::cout << "Nijmegen1 binding energy: -2.224575 MeV" << std::endl;
    
    for (int i = 0; i < (int)diag_res.eigenvalues->size; i++)
    {
        double E = GSL_REAL(gsl_vector_complex_get(diag_res.eigenvalues,i));
        if (E < 0.0)
        {
            std::cout << "This code: " << E << " MeV" << std::endl;
            std::cout << "Relative error: " << (E-E_nijmegen)/E_nijmegen << std::endl; 
        }
    }
    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_complex_free(diag_res.eigenvalues); 


    delete[] pot_V_mtx;
    //delete[] p_grid;
    //delete[] w_grid;
}

void check_MWPC(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string)
{
    std::cout << "Testing phase shifts of LO MWPC code" << std::endl << std::endl;
    
    //std::cout << " These tests are done with the finite p-mesh" << std::endl;
    //std::clock_t start, end;   
    
    // Make grid
    double* p_grid;
    double* w_grid;
    double Lambda = 500.0;
    bool FINITE_GRID = false;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    
    int cut_pow = 6;
    double C1S0	= -0.1/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.13/100.0; // contact term C3S1 for lambda = 450 [MeV]
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    terms.push_back("C3P0");
    terms.push_back("C3P2");


    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false);
    
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }

    // Set correct LECs
    Pot.LECs_["gA2"]  = 1.29*1.29;
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;
    Pot.LECs_["C3P0"] = C3P0;
    Pot.LECs_["C3P2"] = C3P2;


    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    //double rho_T;
    qs::quantum_channel chn = chns[0];
    for (int i = 0; i < (int)chns.size(); i++)
    {
        chn = chns[i];
        if (quantum_channel_to_string(chn) == chn_string) {
            break;
        }
    }
    std::cout << "Computing in chn: " << quantum_channel_to_string(chn) << std::endl;
    if (chn.coupled) {
        std::cout  << "E (MeV) |  uncoup | dm | dp | epsilon (all in deg)" << std::endl;
    } else {
        std::cout  << "E (MeV) |  uncoup (deg)" << std::endl;
    }
    for (int i = 0; i < 350; i++) 
    {
        double Tl = (double)(i+1);
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);
        Phase_shifts_chn phases = solver.solve_in_chn_T(Tl,chn,pot_V_mtx);
        
        double x = 180.0/M_PI;        
        if (chn.coupled) {
            std::cout << Tl << "   " <<  phases.delta_m*x << "   "
                << phases.delta_p*x << "   " << phases.epsilon*x << std::endl;
        } else {
            std::cout  << Tl << "   " << phases.delta_uncoupled*x << std::endl;
        } 
        gsl_matrix_free(pot_V_mtx);
    } 
}


void check_T_matrix(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string)
{
    std::cout << "Testing how the T matrix elements look" << std::endl << std::endl;
    
    //std::cout << " These tests are done with the finite p-mesh" << std::endl;
    //std::clock_t start, end;   
    
    // Make grid
    double* p_grid;
    double* w_grid;
    double Lambda = 500.0;
    bool FINITE_GRID = true;
    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,
                Lambda + 300.0,&p_grid,&w_grid);
    } else 
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }
    int cut_pow = 6;
    double C1S0	= -0.1/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.13/100.0; // contact term C3S1 for lambda = 450 [MeV]
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    terms.push_back("C3P0");
    terms.push_back("C3P2");


    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false);
    
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }

    // Set correct LECs
    Pot.LECs_["gA2"]  = 1.29*1.29;
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;
    Pot.LECs_["C3P0"] = C3P0;
    Pot.LECs_["C3P2"] = C3P2;


    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    //double rho_T
    qs::quantum_channel chn = chns[0];
    for (int i = 0; i < (int)chns.size(); i++)
    {
        chn = chns[i];
        if (quantum_channel_to_string(chn) == chn_string) {
            break;
        }
    }
    std::cout << "Computing in chn: " << quantum_channel_to_string(chn) << std::endl;
    if (chn.coupled) {
        std::cout  << "E (MeV) |  uncoup | dm | dp | epsilon (all in deg)" << std::endl;
    } else {
        std::cout  << "E (MeV) |  uncoup (deg)" << std::endl;
    }
    for (int i = 50; i < 100; i+=1) 
    {
        double Tl = (double)(i+1);
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn, true);
        
        // Solve for the T-matrix
        Phase_shifts_chn phases_T = solver.solve_in_chn_T(Tl,chn,pot_V_mtx);

        // Solve fot the R-matrix
        Phase_shifts_chn phases_R = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
        
        // Solve for the T matrix and take the matrix elemets directly
        std::complex<double>* T_elem = solver.solve_in_chn_T_Telem(Tl, chn, pot_V_mtx);
        double* R_elem = solver.solve_in_chn_R_Relem(Tl, chn, pot_V_mtx);
        std::complex<double>* Tel_from_Rel = solver.T_matrix_from_R_matrix
            (R_elem[0], R_elem[1], R_elem[2], mu, q_on_shell);

        std::vector<Phase_shifts_chn> phases;
        phases.push_back(phases_T);
        std::vector<qs::quantum_channel> chns_vec;
        chns_vec.push_back(chn);
        
        std::vector<std::complex<double>*> T_elem_vec = 
            sc::T_from_phase_shifts(phases, chns_vec,0.0);
        

        double x = 180.0/M_PI;        
        if (chn.coupled) {
            std::cout << "T: " << Tl << "   " <<  phases_T.delta_m*x << "   "
                << phases_T.delta_p*x << "   " << phases_T.epsilon*x << std::endl;
            std::cout << "R: " << Tl << "   " <<  phases_R.delta_m*x << "   "
                << phases_R.delta_p*x << "   " << phases_R.epsilon*x << std::endl;
        } else {
            std::cout  << "T: " << Tl << "   " << phases_T.delta_uncoupled*x << std::endl;
            std::cout  << "R: " << Tl << "   " << phases_R.delta_uncoupled*x << std::endl;
        } 


                
        std::complex<double>* S_R = sc::S_from_Stapp(phases_R.delta_m, phases_R.delta_p,
                phases_R.epsilon);
        std::complex<double>* S_T = sc::S_from_Stapp(phases_T.delta_m, phases_T.delta_p,
                phases_T.epsilon);
        const std::complex<double> imag_u(0.0,1.0);
        // Since I convert to the 2/pi convention earlier
        std::complex<double> fac = -imag_u*4.0*mu*q_on_shell; 
        for (int k=0; k < 3; k++)
        {
            std::cout << "S_T: " << S_T[k] << "  " << std::endl;
            std::cout << "S_R: " << S_R[k] << "  " << std::endl;
            std::cout << "T_dir: " << fac*T_elem[k] << "  " << std::endl;
            std::cout << "T_phase: " << T_elem_vec[0][k] << "  " << std::endl;
            std::cout << "-----" << std::endl;
        }

        std::cout << "Telem: " <<  fac*T_elem[0] << "   " << fac*T_elem[1] << 
            "   " << fac*T_elem[2] << std::endl;
        std::cout << "Tel from R: " <<  fac*Tel_from_Rel[0] << "   " << fac*Tel_from_Rel[1] << 
            "   " << fac*Tel_from_Rel[2] << std::endl;
        std::cout << "Relem: " <<  R_elem[0] << "   " << R_elem[1] << 
            "   " << R_elem[2] << std::endl;
        delete[] S_T;
        delete[] S_R;
        delete[] T_elem;
        gsl_matrix_free(pot_V_mtx);
    } 
}


void check_born(std::string chn_string)
{
    dwba::make_tests(chn_string);
}
void check_DWBA(std::string chn_string)
{
    dwba::make_tests_DWBA(chn_string);
}

void check_NPOT(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot, std::string chn_string)
{
    std::cout << "Testing phase shifts of LO MWPC new pot" << std::endl << std::endl;
    
    //std::cout << " These tests are done with the finite p-mesh" << std::endl;
    //std::clock_t start, end;   
    
    // Make grid
    double* p_grid;
    double* w_grid;
    double Lambda = 500.0;
    bool FINITE_GRID = false;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    
    int cut_pow = 6;
    double C1S0	= -0.1/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.13/100.0; // contact term C3S1 for lambda = 450 [MeV]
    double C3P0 = 0.0;
    double C3P2 = 0.0;    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    terms.push_back("C3P0");
    terms.push_back("C3P2");

    Pot_mwpc<gsl_matrix_complex> Pot_complex = Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false);
    
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot_complex.populate_saved_mtx(chn,true); // Realtivistic factor on
    }

    // Set correct LECs
    Pot_complex.LECs_["gA2"]  = 1.29*1.29;
    Pot_complex.LECs_["C1S0"] = C1S0;
    Pot_complex.LECs_["C3S1"] = C3S1;
    Pot_complex.LECs_["C3P0"] = C3P0;
    Pot_complex.LECs_["C3P2"] = C3P2;

    qs::quantum_channel chn = chns[0];
    for (int i = 0; i < (int)chns.size(); i++)
    {
        chn = chns[i];
        double q_on_shell = 10.0;
        gsl_matrix_complex* pot_V_mtx = Pot_complex.get_saved_matrix(q_on_shell, chn, true);
        ph::print_matrix(pot_V_mtx);
        gsl_matrix_complex_free(pot_V_mtx);
    } 
}
