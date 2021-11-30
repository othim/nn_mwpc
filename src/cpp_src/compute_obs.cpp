/*
 * File to compute observalbes
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

    // Convention factor
    double factor = (M_PI/2.0);
    /*
    for (int i = 0; i < 6; i++)
    {
        std::cout << V_arr[i] << " ";        
        // V_arr[i] = factor*V_arr[i];
    }
    std::cout << std::endl;*/
}
/*
 * Function declarations
 */

void compute_observables(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale,
   double Lambda, double C1S0, double C3S1);

void create_ext_pot();

void compute_1S0(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
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
   // test_f();
    //int a = 0;
    //std::cin >> a;
    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 200.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 96; // Number of points in angular integration
    unsigned int number_of_p_points = 200; // Number of momentum-grid points
    unsigned int J_max_in_pot = 40; // Maximum J that is stored for L-polynomials
    
    // ----- JUST CHOOSE SOME VALUES TO REPRODUCE PHASE SHIFTS WITH -----
    static double Lambda	= 450; 		  // cut-off for renormalization of LO  [MeV]
    static double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
    static double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    int J_max = 5;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;

    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, true);   
    
    // Computing observables
    // compute_observables(chns,number_of_p_points,ang_int_points,J_max_in_pot,scale,Lambda,C1S0,C3S1);
    compute_1S0(chns, number_of_p_points,scale, ang_int_points, J_max_in_pot);
    

    ph::physics_helpers_free();
    return 0;
}


void compute_1S0(std::vector<qs::quantum_channel> chns, unsigned int number_of_p_points, double scale,unsigned int ang_int_points,
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
    qs::quantum_channel chn = chns[0]; // 1S0 channel
    
    // Open a file
    std::ofstream myfile;
    myfile.open("../../data/out_1S0.txt");
    for (int E = 1; E < 301; E++)
    {
        double T_lab = (double)E;

         LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
         
         //gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
         //ph::print_m(pot_V_mtx);
         gsl_matrix* pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn);
         
         Phase_shifts_chn phases = solver.solve_in_chn_R(T_lab,chn,pot_V_mtx);
         gsl_matrix_free(pot_V_mtx);
         std::cout << T_lab << "   " << phases.delta_uncoupled*180.0/M_PI << std::endl;
         myfile << T_lab << "   " <<  phases.delta_uncoupled*180.0/M_PI << "\n"; 
    }
    myfile.close();
}

void compute_observables(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale,
   double Lambda, double C1S0, double C3S1)
{
   std::clock_t start, end;   
   // Make grid
   double* p_grid;
   double* w_grid;
   ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   // Choose terms in the potential, LO WPC
   std::vector<std::string> terms;
   terms.push_back("OPEP"); // To just test elements use just OPEP
   terms.push_back("C1S0");
   terms.push_back("C3S1");

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
   
   Potential_ext nijmegen = Potential_ext(p_grid, number_of_p_points, 450.0, &nijm_correct_arg);

   gsl_matrix* m = nijmegen.get_matrix(10.0,chns[0]);
   std::cout << "Printing potential" << std::endl;
   ph::print_m(m);

   std::cout << "Populate saved matrices..." << std::endl;

   start = std::clock();
   for (auto chn : chns)
   {
      Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
   }
   end = std::clock();
    std::cout << "Time to compute save matrices: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
   LS_Solver solver = LS_Solver(chns,number_of_p_points,scale,true,Lambda,true);
   
   
   start = std::clock();
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   
   unsigned int len = 7;
   double energies[7] = {10.0,20.0,30.0,50.0,80.0,100.0,150.0};

  
   
   double q_on_shell;
   double mu;
   double rho_T;
   
   
   std::cout << "Solving LS equation..." << std::endl;
   std::cout << "Tlab (MeV)   |   cross section (mb)" << std::endl;   
   for (int i = 0; i < len; i++)
   {
      double Tl = energies[i];
      
      std::vector<Phase_shifts_chn> phases_vec;

      start = std::clock();
      for (auto chn : chns)
      {
         LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
         
         gsl_matrix* pot_V_mtx = Pot.get_saved_matrix(q_on_shell, chn,true);
         //gsl_matrix* pot_V_mtx = nijmegen.get_matrix(q_on_shell, chn);
         
         Phase_shifts_chn phases = solver.solve_in_chn_R(Tl,chn,pot_V_mtx);
         gsl_matrix_free(pot_V_mtx);
         phases_vec.push_back(phases);
      }

      end = std::clock();
      std::cout << "Time to compute phase shifts: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
      // Now all the pahse shifts in the relevent channels are known.
      // Now we can compute the total cross section for some on_shell
      // lab energy

      LS_Solver::get_mu_q_on_shell(Tl,chns[0], &mu,&q_on_shell);
   
      rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);

      start = std::clock();
      double cross_section = compute_total_cross_section(chns,phases_vec,q_on_shell,rho_T,30);
      
      // Compute defferential cross section 
      // ----------------------------------
      
      // Get the Saclay amplitudes
      
      // Get the M-matrix 
      
      // Compute observable from the M-matrix

      // ---------------------------------- 
      end = std::clock();
      std::cout << "Time to compute cross section: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
      std::cout << Tl << "\t \t" << cross_section << std::endl;
   } 
}

