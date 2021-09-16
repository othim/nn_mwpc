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

   // Constructing potential with only OPEP-term
   std::vector<std::string> terms;
   terms.push_back("OPEP");
   //terms.push_back("C1S0");

   Potential_mwpc V = Potential_mwpc(terms);
   std::cout << "Size of V: " << sizeof(V) * 8 << " bytes."<< std::endl;
   
   double potential[6];
   V.calc_element_V_arr(10,10,true,2,potential);

   for (int i = 0; i<6; i++)
   {
        printf(" %.15f",potential[i]);
   }
   std::cout << std::endl;

   // Make some tests getting a matrix;
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   gsl_integration_fixed_workspace* int_ang_ = gsl_integration_fixed_alloc(T, 100, -1.0, 1.0, 0, 0);

   double* p_grid = gsl_integration_fixed_nodes(int_ang_);
   double* w_grid = gsl_integration_fixed_weights(int_ang_);
   
   // Make a gsl_gauss_legendre_grid
   std::vector<std::string> terms2;
   terms2.push_back("OPEP");
   Potential_mwpc Pot = Potential_mwpc(terms2,96,p_grid,w_grid,100);

   gsl_matrix* potential_mtx;

   //Pot.populate_saved_mtx(10.0,10.0,true,30,1,true,true);

   std::clock_t start, end;

   start = std::clock();

   Pot.calc_element_V_arr(10,10,true,1,potential);
   
   end = std::clock();
   std::cout << "Time taken to call calc_element_V_arr() is : " << (double)(end-start)/(double)CLOCKS_PER_SEC; 
   std::cout << " sec " << std::endl;
   


   return 0;
}
