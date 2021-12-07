/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include <iomanip> 
#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include <cstdio>
#include <ctime>
#include "scattering.h"
#include "physics_helpers.h"



// Returns true if all tests are passed
/*
bool test_potential_elements(qs::quantum_channel chn,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double T_lab,double scale,double* V_arr_correct,double tol)
{
   // Calculate mass and on-shell momentum
   

   // Make grid
   double* p_grid;
   double* w_grid;
   gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   // Choose terms in the potential
   std::vector<std::string> terms;
   terms.push_back("OPEP"); // To just test elements use just OPEP

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot);
   Pot.LECs_["gA2"] = constants::gA*constants::gA; // Set correct LEC

   double V_arr[6];
   Pot.calc_element_V_arr(q_on_shell,q_on_shell,chn.coupled,chn.J,&V_arr[0]);
   //Pot.calc_element_V_arr(60,60,true,1,&V_arr[0]);
   //for (int i = 0; i < 6; i++)
   //{
   //   std::cout << V_arr[i] << " ";
   //}
   std::cout << std::endl;
   bool passed = true;
   for (int i= 0; i < 6; i++)
   {
      std::cout << fabs((V_arr[i]*pow(2.0*M_PI,3) - V_arr_correct[i])/V_arr_correct[i]) << " ";
      if (fabs((V_arr[i]*pow(2.0*M_PI,3) - V_arr_correct[i])/V_arr_correct[i]) > tol)
      {
         passed = false;
      }
   }
   std::cout << std::endl;
   return passed;
}*/

bool test_potential_matrix(qs::quantum_channel chn)
{
   return false;
}

bool test_phase_shifts(qs::quantum_channel chn,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double T_lab,double scale,double* V_arr_correct,double tol,
   double Lambda, double C1S0, double C3S1)
{
    // Make grid
   double* p_grid;
   double* w_grid;
   ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   /*
   for (int i = 0; i < number_of_p_points; i++)
   {
      std::cout << p_grid[i] << " " << w_grid[i] << std::endl;
   }
   */
   // Choose terms in the potential, LO WPC
   std::vector<std::string> terms;
   terms.push_back("OPEP"); // To just test elements use just OPEP
   terms.push_back("C1S0");
   terms.push_back("C3S1");

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
   Pot.populate_saved_mtx(chn,true); // Realtivistic factor on

   std::vector<qs::quantum_channel> chns; // Do not do anythin
   chns.push_back(chn); // Do not do anything
   LS_Solver solver = LS_Solver(chns,&Pot,number_of_p_points,scale,true,Lambda,true);
   
   std::clock_t start, end;   
   start = std::clock();
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   std::cout << "Getting phase shifts..." << std::endl; 
   Phase_shifts_chn phases = solver.solve_in_chn_R(T_lab,chn,true,true);
   end = std::clock();
   std::cout << "Time taken to compute phase shifts (with R) is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 
   
   
   std::cout << std::setprecision(16) << std::endl << "Phases in Stapp convection (deg): \n" << " delta_m = " << ph::rad_to_deg(phases.delta_m) << "\n delta_p = " << ph::rad_to_deg(phases.delta_p) << 
      "\n epsilon = " << ph::rad_to_deg(phases.epsilon) << "\n delta_uncoupled = " << ph::rad_to_deg(phases.delta_uncoupled) << std::endl << std::endl;
  
   
   // Make same test but with T matrix 
   
   start = std::clock();
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   phases = solver.solve_in_chn_T(T_lab,chn,true,true);
   end = std::clock();
   std::cout << "Time taken to compute phase shifts (with T) is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 
  
   std::cout << std::setprecision(16) << std::endl << "Phases in Stapp convection (deg): \n" << " delta_m = " << ph::rad_to_deg(phases.delta_m) << "\n delta_p = " << ph::rad_to_deg(phases.delta_p) << 
      "\n epsilon = " << ph::rad_to_deg(phases.epsilon) << "\n delta_uncoupled = " << ph::rad_to_deg(phases.delta_uncoupled) << std::endl << std::endl;
   


   if (chn.coupled==false) {
      return (abs(ph::rad_to_deg(phases.delta_uncoupled)+3.061426389773196) < tol);
   } else {
      return (abs(ph::rad_to_deg(phases.delta_m)+75.27581878915144) < tol && 
         abs(ph::rad_to_deg(phases.delta_p)+0.787703968368532) < tol && 
         abs(ph::rad_to_deg(phases.epsilon)-1.637312987120185) < tol);
   }

}

void test_deutron_binding_energy(double Lambda, double C1S0, double C3S1,unsigned int number_of_p_points,double scale,
   unsigned int J_max_in_pot,unsigned int ang_int_points)
{
   // S=1, J=1 coupled=true is the deuteron channel
   qs::quantum_channel deutron_chn = {.J=1, .S=1,.tz=0,.coupled=true};

   // Make grid
   double* p_grid;
   double* w_grid;
   ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   
   // Choose terms in the potential, LO WPC
   std::vector<std::string> terms;
   terms.push_back("OPEP"); 
   terms.push_back("C1S0");
   terms.push_back("C3S1");

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot,450.0);
   
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   gsl_matrix* pot = Pot.get_matrix_no_onshell(deutron_chn,true);
   
   //std::cout << "Printing potential" << std::endl;
   //ph::print_m(pot);
   std::clock_t start, end;   
   start = std::clock();
   ph::eigen_t eig_data = ph::solve_SE(p_grid,w_grid,number_of_p_points,deutron_chn,pot);
   end = std::clock();
   std::cout << "Time to solve_SE: " << (double)(end-start)*1.0e6/(double)CLOCKS_PER_SEC<< " us" << std::endl;
}

void compute_observables(std::vector<qs::quantum_channel> chns,unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double T_lab,double scale,double* V_arr_correct,double tol,
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
         std::cout << "Channel " << i << ": " << ph::rad_to_deg(phases_vec[i].delta_m) << " " << ph::rad_to_deg(phases_vec[i].delta_p) << " " << ph::rad_to_deg(phases_vec[i].epsilon) << 
            " " << ph::rad_to_deg(phases_vec[i].delta_uncoupled) << std::endl;
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
   
   /*
   unsigned int s = 1;
   unsigned int mo = 1;
   unsigned int mi = -1;
   double cos_theta = 1.0;
   unsigned int l_max = 30;


   get_mu_q_on_shell(10.0,chns[0], &mu,&q_on_shell);
   std::cout << chns[0].tz << std::endl;
   std::cout << "q_on_shell: " << q_on_shell << std::endl;
   rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/(constants::Mn+constants::Mp);

   std::cout << "Computing M-matrix" << std::endl;
   std::vector<Phase_shifts_chn> phases_vec;
   for (auto chn : chns)
   {
         Phase_shifts_chn phases = solver.solve_in_chn_R(10.0,chn,true,true);
         phases_vec.push_back(phases);
   }

   std::cout << "Phase shifts (deg) (Stapp):" << std::endl;
   for (int i = 0; i < phases_vec.size(); i++)
   {
      std::cout << "Channel " << i << ": " << ph::rad_to_deg(phases_vec[i].delta_m) << " " << ph::rad_to_deg(phases_vec[i].delta_p) << " " << ph::rad_to_deg(phases_vec[i].epsilon) << 
         " " << ph::rad_to_deg(phases_vec[i].delta_uncoupled) << std::endl;
   }

   std::complex<double> m = get_M_matrix_p(chns,phases_vec,s,mo,mi,cos_theta,q_on_shell,rho_T,l_max);
   std::cout << "M-matrix emement: " << m << std::endl; */
}

void run_tests(qs::quantum_channel chn, unsigned int number_of_p_points,unsigned int ang_int_points,
   unsigned int J_max_in_pot,double T_lab,double scale,double* V_arr_correct,double tol,
   double Lambda, double C1S0, double C3S1)
{
   std::cout << "Testing in channel: " << "J=" << chn.J << " S=" << chn.S << " tz=" << chn.tz << 
      " coupled=" << chn.coupled << std::endl << "-------------------------------------------" << std::endl << std::endl;

   //std::cout << "Testing potential elements" << std::endl;
   //bool pot_test = test_potential_elements(chn,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,V_arr_correct,tol);
   //std::cout << "Test passed: " << pot_test << std::endl << "---------" << std::endl;

   //std::cout << "Testing potential matrix" << std::endl;
   //bool pot_test_mtx = test_potential_matrix(chn);
   //std::cout << "Test passed: " << pot_test_mtx << std::endl << "---------" << std::endl;

   std::cout << "Testing phase shifts:" << std::endl;
   bool test_phase = test_phase_shifts(chn,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,V_arr_correct,tol,Lambda,C1S0,C3S1);
   std::cout << "Test passed: " << test_phase << std::endl << "-------------------------------------------" << std::endl;

   //std::cout << "TEST PASSED: " << (pot_test && pot_test_mtx && test_phase) << std::endl;
   //std::cout << "-----------" << std::endl << "----END----" << std::endl << std::endl;
}

void run_speed_tests(qs::quantum_channel chn)
{

}



int main(int argc, char** argv)
{
    ph::physics_helpers_init();
   // ------ CONSTANTS TO CHANGE ------
   // ---------------------------------
   double scale = 100.0; // Scale of momenutm grid MeV
   unsigned int ang_int_points = 96; // Number of points in angular integration
   unsigned int number_of_p_points = 100; // Number of momentum-grid points
   unsigned int J_max_in_pot = 40; // Maximum J that is stored for L-polynomials
   double T_lab = 10.0; // Lab energy in MeV
   double rel_tol_pot_elements = 1e-4;

   double tol_ps = 1e-11;
   // ----- JUST CHOOSE SOME VALUES TO REPRODUCE PHASE SHIFTS WITH -----
   static double Lambda	= 450; 		  // cut-off for renormalization of LO  [MeV]
   static double C1S0	= -0.112927/100.0; // contact term C1S0 for lambda = 450 [MeV]
   static double C3S1	= -0.087340/100.0; // contact term C3S1 for lambda = 450 [MeV]
   
   // ---------------------------------
   double q_on_shell = sqrt(constants::Mp*constants::Mp*T_lab*(T_lab + 2.0*constants::Mn)/
         ((constants::Mp + constants::Mn)*(constants::Mp + constants::Mn) + 2.0*T_lab*constants::Mp));
   std::cout << "On-shell momentum: " << q_on_shell << std::endl;

   // Test the code in the following channels
   
   qs::quantum_channel chn_uncoup = {.J=1, .S=1,.tz=0,.coupled=false};
   double V_arr_correct1[] = {1.453716658589179581973e-04, 1.175326792502285721664e-04, 0.000000000000000000000e+00, 0.000000000000000000000e+00, -0.000000000000000000000e+00, -0.000000000000000000000e+00};   
   run_tests(chn_uncoup,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,&V_arr_correct1[0],rel_tol_pot_elements,Lambda,C1S0,C3S1);

   qs::quantum_channel chn_coup = {.J=1, .S=0,.tz=0,.coupled=true};
   double V_arr_correct2[] = {  0.000000000000000000000e+00, 0.000000000000000000000e+00, -1.114898705446752692805e-10, -7.129564851709590453523e-10, -3.897949168142152748363e-09, -3.897949168142152748363e-09};   
   run_tests(chn_coup,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,&V_arr_correct2[0],tol_ps,Lambda,C1S0,C3S1);



   
   // Test deuteron binding energy
   std::cout << "\n\nTesting deuteron binding energy" << std::endl;
   test_deutron_binding_energy(Lambda,C1S0,C3S1,number_of_p_points,scale,J_max_in_pot,ang_int_points);

   
   // Test wigxjpf
   // ------------

   double val6j;

   wig_table_init(2*100, 9);
   wig_temp_init(2*100);

   std::clock_t start, end;   
   start = std::clock();
   val6j = wig3jj(2*  1 , 2*  2 , 2*  1 ,
                  2*  0 , 2*  1 , 2*  -1 );
   end = std::clock();
   std::cout << "Time to compute wiegner symbol: " << 1.0e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl ;
 
   wig_temp_free();
   wig_table_free();

   std::cout << val6j << std::endl;


   int size = 100;

   gsl_matrix* gsl_m = gsl_matrix_alloc(size,size);
   for (int i=0; i < size; i++)
   {
      for (int j=0; j < size; j++)
      {
         double r = rand();
         gsl_matrix_set(gsl_m,i,j,r);
      }
   }

    
   // Diagonalize the matrix
   start = std::clock();

   gsl_vector_complex* eval = gsl_vector_complex_alloc(size);
   gsl_eigen_nonsymm_workspace* ws = gsl_eigen_nonsymm_alloc(size);
   gsl_eigen_nonsymm(gsl_m,eval,ws);
   
   end = std::clock();
   std::cout << "Time to diag matrix: " << 1.0e6*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl ;
   
   // Compute observables
   // -------------------

   // Construct the quantum states
   std::cout << "Constructing quantum states:" << std::endl;
   unsigned int J_max = 30;
   unsigned int J_min = 0;
   int Tz_min = 0;
   int Tz_max = 0;
   bool print = true;

   // Construct the quantum states
   std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
   
   // Construct the quantum scattering channels from the states
   std::vector<qs::quantum_channel> chns = get_channels(states, true);   
   compute_observables(chns,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,&V_arr_correct1[0],1.0e-6,Lambda,C1S0,C3S1);
    
    /* Test kronecker product
     */
    gsl_matrix_complex* sigma_x = ph::get_Pauli_matrix('z');
    gsl_matrix_complex* sigma_y = ph::get_Pauli_matrix('z');
    gsl_matrix_complex* prod = ph::kronecker_product(sigma_x,sigma_y);
    
    gsl_complex t = ph::trace(prod);
    std::cout << "trace=(" << GSL_REAL(t) << "," << GSL_IMAG(t) << ")" << std::endl; 
    ph::print_m_complex(prod);
   // -----------------
   // ------OLD--------
   
   
   /*
   // Testing LS_Solver
   // -----------------
   std::cout << "TESTING LS-SOLVER" << std::endl;
     
   Potential_mwpc pot2 = Potential_mwpc(terms2,96,p_grid_2,w_grid_2,number_of_points,40);
   pot2.populate_saved_mtx(chn2,true,true);

   Pot.LECs_["gA2"] = constants::gA*constants::gA;

   LS_Solver solver = LS_Solver(chns,&pot2,number_of_points);
   start = std::clock();

   pot2.LECs_["gA2"] = constants::gA*constants::gA;
   Phase_shifts_chn phases = solver.solve_in_chn(68.54,chn2,true,true);

   end = std::clock();
   std::cout << "Time taken to compute phase shifts is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
   std::cout << " mu sec " << std::endl;
   
   std::cout << "Phases: " << "delta_p=" << ph::rad_to_deg(phases.delta_p) << " delta_m=" << ph::rad_to_deg(phases.delta_m) << 
      " epsilon=" << ph::rad_to_deg(phases.epsilon) << " delta_uncoupled=" << ph::rad_to_deg(phases.delta_uncoupled) << " deg" << std::endl;



   //LS_Solver solver = LS_Solver(100,&Pot);
   
   //solver.solve_in_chn(10,chn,false,false);*/
   ph::physics_helpers_free();
   return 0;
}


