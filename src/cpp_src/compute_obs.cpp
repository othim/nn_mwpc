/*
 * File to compute observalbes
 */


#include <iostream>
#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" 
#include <cstdio>
#include <ctime>
#include "scattering.h"
#include "physics_helpers.h"


int main(int argc, char** argv)
{
    // Define constants
    
    // Do precomputations

    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    unsigned int J_max = 30;
    unsigned int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;

    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, true);   
    compute_observables(chns,number_of_p_points,ang_int_points,J_max_in_pot,scale,Lambda,C1S0,C3S1);

    return 0;
}


void compute_observables(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double scale,
   double Lambda, double C1S0, double C3S1)
{
   std::clock_t start, end;   
   // Make grid
   double* p_grid;
   double* w_grid;
   gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   // Choose terms in the potential, LO WPC
   std::vector<std::string> terms;
   terms.push_back("OPEP"); // To just test elements use just OPEP
   terms.push_back("C1S0");
   terms.push_back("C3S1");

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
   
   std::cout << "Populate saves matrices..." << std::endl;

   start = std::clock();
   for (auto chn : chns)
   {
      Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
   }
   end = std::clock();
    std::cout << "Time to compute save matrices: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
   LS_Solver solver = LS_Solver(chns,&Pot,number_of_p_points,scale,true,Lambda,true);
   
   
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
         Phase_shifts_chn phases = solver.solve_in_chn_R(Tl,chn,true,true);
         phases_vec.push_back(phases);
      }

      /*std::cout << "Phase shifts (deg) (Stapp):" << std::endl;
      for (int i = 0; i < phases_vec.size(); i++)
      {
         std::cout << "Channel " << i << ": " << rad_to_deg(phases_vec[i].delta_m) << " " << rad_to_deg(phases_vec[i].delta_p) << " " << rad_to_deg(phases_vec[i].epsilon) << 
            " " << rad_to_deg(phases_vec[i].delta_uncoupled) << std::endl;
      }*/
      end = std::clock();
      std::cout << "Time to compute phase shifts: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
      // Now all the pahse shifts in the relevent channels are known.
      // Now we can compute the total cross section for some on_shell
      // lab energy
      

      get_mu_q_on_shell(Tl,chns[0], &mu,&q_on_shell);
   
      rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);

      start = std::clock();
      double cross_section = compute_total_cross_section(chns,phases_vec,q_on_shell,rho_T,30);
      end = std::clock();
      std::cout << "Time to compute cross section: " << 1e3*(end-start)/(double)CLOCKS_PER_SEC << " ms" << std::endl;
   
      std::cout << Tl << "\t \t" << cross_section << std::endl;
   } 
}
