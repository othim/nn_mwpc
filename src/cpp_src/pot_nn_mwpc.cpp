#include "pot_nn_mwpc.h"
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" // GL integration
#include "Constants.h"
#include <iostream>
// Constructor
Potential_mwpc::Potential_mwpc()
{
   // Make grid for angular integration
   
   std::cout << "V object created" << std::endl;
}

// Destructor
Potential_mwpc::~Potential_mwpc()
{
   // Free all memory allocations
   std::cout << "V object deleted" << std::endl;
}

/* 
   ** Some helper functions to opep_get_el **
   ******************************************
*/


/*
   V_\alpha from Erkelenz without the isospin factor
*/
double pot_OPEP_mom(double qi, double qo, double z)
{
   double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	return -(gA*gA/(4.0*fpi*fpi))*(1.0/(q2+mpi*mpi));
}

/*
   This function returns A^{J,l}(q',q) in Erkelenz, where V_\alpha is for ope.
*/
double ang_Integral(double qi, double qo, int J,int l)
{
   // Make grid
   const gsl_integration_fixed_type * T = gsl_integration_fixed_hermite;
   gsl_integration_fixed_workspace* integration = 
   gsl_integration_fixed_alloc(T, 100, 1, 2,1,2);
   
   return 0.0;
}

/*
   Computes the total isospin factor from \tau_1 \cdot \tau_2 from
   the constraint J+L+T = odd 
*/
int isoFac(int L,int S)
{
   return 0;
}

void Potential_mwpc::opep_get_el(double qi,double qo, bool coupled, int J, double* output)
{

   // Precompute some integrals

   // Check which elements that are non-zero by checking if the Channel is
   // coupled or not

   // Make the output point to the values

}
