/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include "pot_nn_mwpc.h"
#include <cstdio>
#include <ctime>


int main(int argc, char** argv)
{
   std::cout << "Running test.cpp" << std::endl;


   // Make some tests getting a matrix
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, 100, -1.0, 1.0, 0, 0);

   double* p_grid = gsl_integration_fixed_nodes(int_ang_);
   double* w_grid = gsl_integration_fixed_weights(int_ang_);
   
   // Make a gsl_gauss_legendre_grid
   std::vector<std::string> terms2;
   terms2.push_back("OPEP");
   Potential_mwpc Pot = Potential_mwpc(terms2,96,p_grid,w_grid,100,20);

   double potential[6];
   gsl_matrix* potential_mtx;
   Pot.populate_saved_mtx(10.0,10.0,true,1,1,true,true);
 
   
   std::clock_t start, end;
   
   start = std::clock();

   Pot.LECs_["gA2"] = constants::gA*constants::gA;
   Pot.calc_element_V_arr(10,10,true,1,potential);
   //Pot.get_saved_matrix(10,true,1,1,true,true,potential_mtx);
   end = std::clock();
   std::cout << "Time taken to call calc_element_V_arr() is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
   std::cout << " mu sec " << std::endl;
   
   for (int i= 0; i < 6; i++)
   {
      std::cout << potential[i] << " ";
   }
   std::cout << std::endl;


   return 0;
}
