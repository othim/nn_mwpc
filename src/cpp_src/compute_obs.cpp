/*
 * File to compute observalbes and check the speed of the code.
 */


#include <iostream>
#include <fstream>
#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" 
#include <cstdio>
#include <ctime>
#include "scattering.h"
#include "physics_helpers.h"
#include "pot_ext.h"
#include <algorithm>
/*
 * This function can be called if this file is linked with 
 * the .o files from the fortran libray compiled.
 */
extern "C" {
    void nijmegen_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

// This function is not complete!!!
void nijm_correct_arg(double qi, double qo, bool coupled, int S, int J, int T, int Tz,  double* V_arr)
{
    int coup = (int)coupled;
    nijmegen_fort_interface(&qi, &qo, &coup, &S, &J, &T, &Tz, &V_arr[0]); 

}
/*
 * Function declarations
 */

void check_observable(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale, std::string obs_string, bool ope_J_geq9);

void create_ext_pot();

void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot);


void check_speed(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot);


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

int main(int argc, char** argv)
{

    using namespace sc;
    // test_f();
    //int a = 0;
    //std::cin >> a;
    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 96; // Number of points in angular integration
    unsigned int number_of_p_points = 100; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    
    // ----- JUST CHOOSE SOME VALUES TO REPRODUCE PHASE SHIFTS WITH -----
    //static double Lambda	= 450; 		  // cut-off for renormalization of LO  [MeV]
    //static double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    //static double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    int J_max =8;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;
    bool OPE_inclue = true;
    
    
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, true);   
    
    // Computing observables
    // compute_observables(chns,number_of_p_points,ang_int_points,J_max_in_pot,scale,Lambda,C1S0,C3S1);
    if (std::string(argv[1]) == "phase") {
        check_phase_shifts(chns, number_of_p_points,scale, ang_int_points, J_max_in_pot);
    }
    //std::complex<double> a =( -1.70140, 8.83681);
    //std::complex<double> e = (0.00069, 0.00306);
    //std::cout << std::real(std::conj(a)*e) << std::endl;
    
    // Check observables
    if (std::string(argv[1]) == "DSG") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"I 0000", OPE_inclue);
    } else if (std::string(argv[1]) == "PB") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"P n000", OPE_inclue);
    } else if (std::string(argv[1]) == "CKK") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"A 00mm", OPE_inclue);
    } else if (std::string(argv[1]) == "AYY") {
        check_observable(chns, number_of_p_points, scale, ang_int_points, J_max_in_pot,"C nn00", OPE_inclue);
    } else if (std::string(argv[1]) == "SPEED") {
        check_speed(chns, number_of_p_points, scale,ang_int_points, J_max_in_pot);
    }
    ph::physics_helpers_free();
    return 0;
}

void check_phase_shifts(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
   unsigned int J_max_in_pot)
{
    double* p_grid;
    double* w_grid;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

    double Lambda = 5000.0;

    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
  
    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
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

    LS_Solver solver = LS_Solver(chns,number_of_p_points,scale,true,Lambda,true);

    double mu;
    double q_on_shell;
    

    for (int c_i = 0; c_i < chns.size(); c_i++)
    {
    
        qs::quantum_channel chn = chns[c_i]; // 1S0 channel
        // Open a file
        std::ofstream myfile;
        std::string filename = "../../data/out_" + quantum_channel_to_string(chn) + ".txt"; 
        myfile.open(filename);
        std::cout << "File_name: " << filename << std::endl; 
        

        // Read in the correct file of data
        std::string data = "../../data/np_" + quantum_channel_to_string(chn) + "_nijm1.txt";   
        
        // Open file
        std::ifstream infile(data);
        if (infile.is_open())
        {
            std::cout << "OK" << std::endl;
        } else
        {
            std::cout << "Failed" << std::endl;
        }
        double D_energies[300];
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
                D_energies[i] = E;
                D_delta_uncoupled[i] = d;
                i++;
            }
        } else
        {
            double E, dp, dm, e;
            int i = 0;
            while(infile >> E >> dm >> dp >> e)
            {
                D_energies[i] = E;
                D_delta_m[i] = dm;
                D_delta_p[i] = dp;
                D_eps[i] = e;
                i++;
            }
        }


        std::cout << "E_lab d-uncoup \t dm \t dp \t eps " << std::endl; 
        double error = 0;
        double error_m = 0;
        double error_p = 0;
        double error_eps = 0;

        for (int E = 1; E < 301; E++)
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
                std::cout << T_lab << "   " << std::abs(phases.delta_uncoupled*180.0/M_PI - D_delta_uncoupled[E-1]) << std::endl;
                error += std::abs((phases.delta_uncoupled*180.0/M_PI - D_delta_uncoupled[E-1])/D_delta_uncoupled[E-1]); 
            } else 
            {
                double em = std::abs((phases.delta_m*180.0/M_PI - D_delta_m[E-1])/D_delta_m[E-1]);         
                double ep = std::abs((phases.delta_p*180.0/M_PI   - D_delta_p[E-1])/D_delta_p[E-1]);         
                double eps = std::abs((phases.epsilon*180.0/M_PI  - D_eps[E-1])/D_eps[E-1]);         
                std::cout << T_lab << "   -   " << em << "   " << ep << "   " << eps << std::endl;
                error_m += em;
                error_p += ep;
                error_eps += eps;
            } 
        
        }
        std::cout << "Errors: " << error/300 << "   " << error_m/300 << "   " << error_p/300 << "   " << error_eps  << std::endl;
        
        myfile.close();
        //double a;
        //std::cout << "One channel done" << std::endl;
        //std::cin >> a;
    }
}



void check_observable(std::string observable, double energy, Potential_ext& pot, LS_Solver& LS_Solver)
{


} 



void check_observable(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale,std::string obs_string, bool ope_J_geq_9)
{
    //std::clock_t start, end;   
    // Make grid
    double* p_grid;
    double* w_grid;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    // Choose terms in LO WPC potential
    
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");

    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
   
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
    }
  
    Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;

    // Construct potential and LS-Solver 

    // Choose terms in OPE potential
    
    double Lambda = 10000.0;
    int l_max = 50;
    std::vector<std::string> terms2;
    terms2.push_back("OPEP"); // To just test elements use just OPEP

    Potential_mwpc OPE = Potential_mwpc(terms2,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,Lambda);
  
    
    for (auto chn : chns)
    {
        if (chn.J > 9)
        {
            OPE.populate_saved_mtx(chn,true); // Realtivistic factor on
        }
    }
  
    //OPE.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
    OPE.LECs_["gA2"]  = 1.1*1.1; // Set correct LEC

    // Compute the observables
    std::string obs_string2;

    if (obs_string == "I 0000") {
        obs_string2 = "DSG";
    } else if (obs_string == "P n000") {
        obs_string2 = "PB";
    } else if (obs_string == "A 00mm") {
        obs_string2 = "CKK";
    } else if (obs_string == "C nn00") {
        obs_string2 = "AYY";
    }

    Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, Lambda, &nijm_correct_arg);
    LS_Solver solver = LS_Solver(chns,number_of_p_points,scale,true,Lambda,true);
   
    double q_on_shell;
    double mu;
    double rho_T;
   
    const int len = 3;
    //double energies[len] = {10.0, 50.0, 200.0};
    double energies[len] = {125.0, 200.0, 350.0};
    
    for (int i = 0; i < len; i++)
    {
        double Tl = energies[i];
        // Compute all the phase shifts in the channels
        std::vector<Phase_shifts_chn> phases_vec;
        //start = std::clock();
        std::cout << std::endl << "Testing " + obs_string2 + " with T_lab=" << Tl << " MeV" << std::endl << std::endl;  
    
        // Read in the correct file of data
        std::string data = "../../data/np_" + obs_string2 + "_" + std::to_string((int)Tl) + "_nijm1.txt";   
        
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
        double D_theta[180];
        double D_obs[180];
        double theta, obs;
        int k = 0;
        while(infile >> theta >> obs)
        {
            D_theta[k] = theta;
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
        std::string filename = "../../data/out_" + obs_string2 + "_" + std::to_string((int)Tl) + ".txt"; 
        myfile.open(filename);

        // Compute the observables
        std::cout << "Angle \t obs \t correct \t abs. rel. error" << std::endl;    

        double errors[180];
        double mean_error = 0;
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

            double obs = sc::compute_observable(saclay_amplitudes, obs_string);
            
            //std::cout << D_obs[ang-1] << " " << obs << std::endl; 
            if (D_obs[ang-1] != 0) {
                errors[ang-1] = std::abs((D_obs[ang-1] - obs)/D_obs[ang-1]);
            } else {
                errors[ang-1] = 0;
            }
            std::cout << angle << " a: " << saclay_amplitudes[0] << " b: " << saclay_amplitudes[1] <<
               " c: " << saclay_amplitudes[2] << " d: " << saclay_amplitudes[3] << " e: " << saclay_amplitudes[4] << std::endl; 
            
            std::cout << angle << "\t" << obs << "\t" << D_obs[ang-1] << "\t" << errors[ang-1]  << std::endl;
            myfile << angle << "\t" << obs << "\t" << D_obs[ang-1] << "\t" << errors[ang-1] << std::endl;
            mean_error += errors[ang-1];

            /*
            // With the M-matrix
            double obs_M;
            if (obs_string == "I 0000")
            {
                gsl_matrix_complex* M_matrix = get_M_matrix(chns, phases_vec, q_on_shell, angle*M_PI/180.0, rho_T, l_max);
                gsl_matrix_complex* eig = gsl_matrix_complex_alloc(2,2);
                gsl_matrix_complex_set_identity(eig);
                ph::print_m_complex(eig);
                ph::print_m_complex(M_matrix);
                obs_M = get_observables(eig,eig,eig,eig,M_matrix); 
            }            
            std::cout << angle << "   " << obs_M << " mb" << std::endl; */
        }
        std::cout << "Mean absolute relative error: " << mean_error/180.0  << std::endl;
        std::cout << "Maximum error: " << *(std::max_element(errors, errors + 180)) << std::endl;
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
}

void check_speed(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, 
        double scale,unsigned int ang_int_points, unsigned int J_max_in_pot)
{
    std::cout << "Testing speed of code with LO WPC potential and the observable DSG" << std::endl << std::endl;

    std::clock_t start, end;   
    
    // Make grid
    double* p_grid;
    double* w_grid;
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

    double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    
    // Choose terms in LO WPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");

    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
    
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
    double Lambda = 450.0;

    LS_Solver solver = LS_Solver(chns,number_of_p_points,scale,true,Lambda,true);
   
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
        double obs = sc::compute_observable(saclay_amplitudes, "I 0000");
        end = std::clock();
        if (ang == 1) {
            std::cout << "Compute OBS from Sac. amp.: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl << std::endl;
        }
    } 
    
    //end = std::clock();
    //std::cout << "Time calculate OBS for 180 angles: " << 1e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl;
}
