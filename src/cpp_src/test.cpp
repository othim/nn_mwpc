/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" 
#include <cstdio>
#include <ctime>

void gauss_legendre_inf_mesh(unsigned int Number_of_points, double scale,double** p,double** w)
{
    // Make grid from -1 to 1
   
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, Number_of_points, -1.0, 1.0, 0, 0);
   
   double* p_grid = gsl_integration_fixed_nodes(int_ang_);
   double* w_grid = gsl_integration_fixed_weights(int_ang_);
   // Make transformation
   double pi_4 = M_PI/4.0;


   for (int i = 0; i < Number_of_points; i++) 
   {
      std::cout << p_grid[i] << " " << w_grid[i] << std::endl;
   }
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


void test_computing_phase_shifts()
{
   // Test potential elements

   // Test potential matrix elements

   // Test phase shifts
}

double rad_to_deg(double in)
{
   return in*180.0/M_PI;
}

int main(int argc, char** argv)
{

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
   Pot.populate_saved_mtx(10.0,10.0,chn,true,true);
  
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
   pot2.populate_saved_mtx(10.0,10.0,chn2,true,true);

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
   //solver.solve_in_chn(10,chn,false,false);
   return 0;
}
