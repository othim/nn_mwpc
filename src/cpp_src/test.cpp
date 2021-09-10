/*
   Test file to be able to run the C++ code
*/
#include <iostream>
#include "pot_nn_mwpc.h"

int main(int argc, char** argv)
{
   std::cout << "Running test.cpp" << std::endl;

   Potential_mwpc* V = new Potential_mwpc();

   delete V;

   return 0;
}
