/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include <iomanip> 
#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" 
#include <cstdio>
#include <ctime>
#include "gsl_eigen.h"

// This data type contains information about eigenvalues and
// eigenvectors of a given matrix.
typedef struct 
{
   gsl_vector_complex* eigenvalues;
   gsl_matrix_complex* eigenvectors;
} eigen_t;
eigen_t solve_SE(double* p, double* w, unsigned int numer_of_grid_points,qs::quantum_channel chn, const gsl_matrix* V);


double rad_to_deg(double in)
{
   return in*180.0/M_PI;
}

void gauss_legendre_inf_mesh(unsigned int Number_of_points, double scale,double** p,double** w)
{
    // Make grid from -1 to 1
   
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, Number_of_points, -1.0, 1.0, 0, 0);
   
   double* p_grid = gsl_integration_fixed_nodes(int_ang_);
   double* w_grid = gsl_integration_fixed_weights(int_ang_);
   // Make transformation
   double pi_4 = M_PI/4.0;

   double* pp = (double*)malloc(Number_of_points*sizeof(double));
   double* ww = (double*)malloc(Number_of_points*sizeof(double));
  
   for (int i = 0; i < Number_of_points; i++)
   {
      double x = p_grid[i];
      pp[i] = scale*tan(pi_4*(x+1));
      ww[i] = (scale*pi_4/(cos(pi_4*(x+1))*cos(pi_4*(x+1))))*w_grid[i];
   }
   *p = pp;
   *w = ww;
}

void print_m(gsl_matrix* matrix)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < matrix->size1; i++)
   {
      for (std::size_t j = 0; j < matrix->size1; j++)
      {
         std::cout << gsl_matrix_get(matrix,i,j) << " ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}



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
   gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

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
   Phase_shifts_chn phases = solver.solve_in_chn_R(T_lab,chn,true,true);
   end = std::clock();
   std::cout << "Time taken to compute phase shifts (with R) is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 
   
   
   std::cout << std::setprecision(16) << std::endl << "Phases in Stapp convection (deg): \n" << " delta_m = " << rad_to_deg(phases.delta_m) << "\n delta_p = " << rad_to_deg(phases.delta_p) << 
      "\n epsilon = " << rad_to_deg(phases.epsilon) << "\n delta_uncoupled = " << rad_to_deg(phases.delta_uncoupled) << std::endl << std::endl;
  
   
   // Make same test but with T matrix 
   
   start = std::clock();
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   phases = solver.solve_in_chn_T(T_lab,chn,true,true);
   end = std::clock();
   std::cout << "Time taken to compute phase shifts (with T) is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC << " us" << std::endl; 
  
   std::cout << std::setprecision(16) << std::endl << "Phases in Stapp convection (deg): \n" << " delta_m = " << rad_to_deg(phases.delta_m) << "\n delta_p = " << rad_to_deg(phases.delta_p) << 
      "\n epsilon = " << rad_to_deg(phases.epsilon) << "\n delta_uncoupled = " << rad_to_deg(phases.delta_uncoupled) << std::endl << std::endl;
   


   if (chn.coupled==false) {
      return (abs(rad_to_deg(phases.delta_uncoupled)+3.061426389773196) < tol);
   } else {
      return (abs(rad_to_deg(phases.delta_m)+75.27581878915144) < tol && 
         abs(rad_to_deg(phases.delta_p)+0.787703968368532) < tol && 
         abs(rad_to_deg(phases.epsilon)-1.637312987120185) < tol);
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
   gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);

   
   // Choose terms in the potential, LO WPC
   std::vector<std::string> terms;
   terms.push_back("OPEP"); 
   terms.push_back("C1S0");
   terms.push_back("C3S1");

   Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,number_of_p_points,J_max_in_pot);
   
   Pot.LECs_["gA2"]  = constants::gA*constants::gA; // Set correct LEC
   Pot.LECs_["C1S0"] = C1S0;
   Pot.LECs_["C3S1"] = C3S1;
   gsl_matrix* pot = Pot.get_matrix_no_onshell(deutron_chn,true); // Should it be true here? probably not
  
   eigen_t eig_data = solve_SE(p_grid,w_grid,number_of_p_points,deutron_chn,pot);

   // Print eigenvalues



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
   
   qs::quantum_channel chn_uncoup = {.J=1, .S=0,.tz=0,.coupled=false};
   double V_arr_correct1[] = {1.453716658589179581973e-04, 1.175326792502285721664e-04, 0.000000000000000000000e+00, 0.000000000000000000000e+00, -0.000000000000000000000e+00, -0.000000000000000000000e+00};   
   run_tests(chn_uncoup,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,&V_arr_correct1[0],rel_tol_pot_elements,Lambda,C1S0,C3S1);

   qs::quantum_channel chn_coup = {.J=1, .S=0,.tz=0,.coupled=true};
   double V_arr_correct2[] = {  0.000000000000000000000e+00, 0.000000000000000000000e+00, -1.114898705446752692805e-10, -7.129564851709590453523e-10, -3.897949168142152748363e-09, -3.897949168142152748363e-09};   
   run_tests(chn_coup,number_of_p_points,ang_int_points,J_max_in_pot,T_lab,scale,&V_arr_correct2[0],tol_ps,Lambda,C1S0,C3S1);

   
   // Test deuteron binding energy
   std::cout << "\n\nTesting deuteron binding energy" << std::endl;
   test_deutron_binding_energy(Lambda,C1S0,C3S1,number_of_p_points,scale,J_max_in_pot,ang_int_points);

   // Test the speed of some calculations

   // -----------------
   // ------OLD--------
   /*
   qs::quantum_channel chn2= {.J=1, .S=0,.tz=0,.coupled=true};
   std::vector<qs::quantum_channel> chns;
   chns.push_back(chn2);
   unsigned int number_of_points = 100;
   double scale = 100.0;
   double* p_grid_2;
   double* w_grid_2;
   std::cerr << "making mesh" << std::endl;
   gauss_legendre_inf_mesh(number_of_points,scale,&p_grid_2,&w_grid_2);
   

   std::cout << "Running test.cpp" << std::endl;


   // Make some tests getting a matrix
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, number_of_points, -1.0, 1.0, 0, 0);

   double* p_grid = gsl_integration_fixed_nodes(int_ang_);
   double* w_grid = gsl_integration_fixed_weights(int_ang_);
   

   // Make a gsl_gauss_legendre_grid
   std::vector<std::string> terms2;
   terms2.push_back("OPEP");

   qs::quantum_channel chn= {.J=1, .S=1,.tz=0,.coupled=true};
  
   Potential_mwpc Pot = Potential_mwpc(terms2,96,p_grid_2,w_grid_2,number_of_points,40);
   Pot.populate_saved_mtx(chn,true,true);
  
   double potential[6];
   gsl_matrix* potential_mtx;
  
   
   std::clock_t start, end;
   
   start = std::clock();

   Pot.LECs_["gA2"] = constants::gA*constants::gA;
   
   potential_mtx = Pot.get_saved_matrix(100,chn,true,true);
   
   end = std::clock();
   std::cout << "Time taken to get_saved_matrix() is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
   std::cout << " mu sec " << std::endl;
   for (int i= 0; i < 6; i++)
   {
      std::cout << potential[i] << " ";
   }
   std::cout << std::endl;

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
   
   std::cout << "Phases: " << "delta_p=" << rad_to_deg(phases.delta_p) << " delta_m=" << rad_to_deg(phases.delta_m) << 
      " epsilon=" << rad_to_deg(phases.epsilon) << " delta_uncoupled=" << rad_to_deg(phases.delta_uncoupled) << " deg" << std::endl;



   //LS_Solver solver = LS_Solver(100,&Pot);
   //solver.solve_in_chn(10,chn,false,false);*/
   return 0;
}


// Implementation


eigen_t solve_SE(double* p, double* w, unsigned int number_of_grid_points,qs::quantum_channel chn, const gsl_matrix* V)
{
   // The potential is assumed to be in a partial wave basis with normalization 
   // <p'|p> = (pi/2)*\delta(p'-p)/p^2 (as in Landau).
   if (chn.coupled)
   {
      if (V->size1 != 2*number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   } else 
   {
     if (V->size1 != number_of_grid_points)
      {
         std::cerr << "Error in solve_SE(): Number of grid points do not match potential dimensions" << std::endl;
      }
   }
   // Get reduced mass of system
   double mu;
   if (chn.tz == -1)
   {
      mu = constants::Mn/2.0; // nn
   } else if (chn.tz == 0)
   {
      mu = constants::Mp*constants::Mn/(constants::Mn+constants::Mp); // np
   } else if (chn.tz == 1)
   {
      mu = constants::Mp/2.0; // pp
   } else { 
      std::cout << "Error: unknown isospin" << std::endl;
   }
   
   gsl_matrix* H = gsl_matrix_alloc(V->size1,V->size2);
   
   // Construct Hamiltonian
   for (int i = 0; i < H->size1; i++)
   {
      for (int j=0; j < H->size2; j++)
      {
         // This is to still use the same momenta
         int l = j;
         if (!(j<number_of_grid_points))
         {
            l = j-number_of_grid_points;
         }
         double p2 = p[l]*p[l];
         if (i==0) {
            //std::cout << p2 << std::endl;
         }
         //double el = gsl_matrix_get(V,i,j)*p[l]*p[k]*sqrt(w[j]*w[k]);
         double el = gsl_matrix_get(V,i,j)*p2*w[l];
         if (i==j) {
            el += p2/(2.0*mu);
         }
         gsl_matrix_set(H,i,j,el);
      }
   }
   //print_m(H);

   // Diagonalize the matrix
   gsl_vector_complex* eval = gsl_vector_complex_alloc(V->size1);
   gsl_eigen_nonsymm_workspace* ws = gsl_eigen_nonsymm_alloc(V->size1);
   gsl_eigen_nonsymm(H,eval,ws);
   
   // Use a test matrix
   /* finds eigenvalues to this...
   gsl_matrix* H_test = gsl_matrix_alloc(2,2);
   gsl_matrix_set(H_test,0,0,1);
   gsl_matrix_set(H_test,0,1,0);
   gsl_matrix_set(H_test,1,0,0);
   gsl_matrix_set(H_test,1,1,1);

   eval = gsl_vector_complex_alloc(2);
   gsl_eigen_nonsymm_workspace* ws_test = gsl_eigen_nonsymm_alloc(2);
   gsl_eigen_nonsymm(H_test,eval,ws_test);
   */
   
   
   for (int i = 0; i < V->size1; i++)
   {
      if (GSL_REAL(gsl_vector_complex_get(eval,i)) < 0)
      {
         std::cout << "(" << GSL_REAL(gsl_vector_complex_get(eval,i)) << " " << GSL_IMAG(gsl_vector_complex_get(eval,i)) << ")" << std::endl;
      }
   }

   eigen_t e;
   return e;
}
