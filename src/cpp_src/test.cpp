/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include "pot_nn_mwpc.h"

int main(int argc, char** argv)
{
   std::cout << "Running test.cpp" << std::endl;

   Potential_mwpc* V = new Potential_mwpc();

   double potential[6];
   V->opep_get_el(100,100,true,3,potential);
   
   for (int i = 0; i<6; i++)
   {
        printf(" %.15f",potential[i]);
   }
   std::cout << std::endl;
   
   delete V;

   return 0;
}
