#include "pot_nn_mwpc.h"
#include "Constants.h"
#include <iostream>
#include <cmath>

// Constructor
Potential_mwpc::Potential_mwpc()
{
   // Make grid for angular integration
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   int_ang = gsl_integration_fixed_alloc(T, N_GLI_PWA, -1.0, 1.0, 0, 0);
   
   double* d = gsl_integration_fixed_nodes(int_ang);
   for (int i = 0; i < 96; i++)
   {
      std::cout << d[i] << " ";
   }
   std::cout << "V object created" << std::endl;
}

// Destructor
Potential_mwpc::~Potential_mwpc()
{
   // Free all memory allocations
   gsl_integration_fixed_free(int_ang);
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

// This struct is to define what set of parameters the function f_int want
struct my_f_params { double qi; double qo; int J; int l; };

double f_int(double z, void* p)
{
   // Decode the void* to the parameters
   struct my_f_params * params = (struct my_f_params *)p;
   double qi = (params->qi);
   double qo = (params->qo);
   int J = (params->J);
   int l = (params->l);
   
   return pot_OPEP_mom(qi,qo,z)*gsl_pow_int(z,l)*gsl_sf_legendre_Pl(J, z);
}

/*
   This function returns A^{J,l}(q',q) in Erkelenz, where V_\alpha is for ope.
*/
double ang_Integral(double qi, double qo, int J,int l,gsl_integration_fixed_workspace* int_ang)
{ 
   // Define function to integrate
   gsl_function F;
   F.function = &f_int;
   struct my_f_params params = {qi, qo, J, l};
   F.params = &params;

   double result;
   // Perform the integration from -1 to 1
   gsl_integration_fixed(&F,&result,int_ang);

   return result*M_PI;
}

/*
   Computes the total isospin factor from \tau_1 \cdot \tau_2 from
   the constraint J+L+T = odd 
*/
int isoFac(int L,int S)
{
   int T = (1-L-S & 1); // L+S+T=odd
   return -3*(1-T) + T; // Return factor from \tau_i \cdot \tau_2 in terms of T
}


void Potential_mwpc::opep_get_el(double qi,double qo, bool coupled, int J,double* output)
{

   // Precompute some integrals
   double integral_0 = ang_Integral(qi,qo,J,0,int_ang);
	double integral_P = ang_Integral(qi,qo,J+1,0,int_ang);
   double integral_M = 0;
   if (J!=0) {
      integral_M = ang_Integral(qi,qo,J-1,0,int_ang);
   }

    // Define some variables
   double V_uncoupled_S0 = 0;
	double V_uncoupled_S1 = 0;
	double V_coupled_mm   = 0;
	double V_coupled_pm   = 0;
	double V_coupled_mp   = 0;
	double V_coupled_pp   = 0;

   // Check which elements that are non-zero by checking if the Channel is
   // coupled or not

   if (!coupled)
   {
      double integral_1 = ang_Integral(qi,qo,J,1,int_ang);
      // OPEP uncoupled interactions
      V_uncoupled_S0 = 2 * isoFac(J,0) * (-(qo*qo+qi*qi)*integral_0 + 2*qo*qi*integral_1);
      if (J!=0) {
		   V_uncoupled_S1 = 2 * isoFac(J,1) * ((qo*qi+qi*qi)*integral_0 - 2*qo*qi*(1.0/(2.0*J+1.0))*(J*integral_P + (J+1)*integral_M));
      }
   } else
   {
      V_coupled_pp    = isoFac(J+1,1)*(2.0/(2.0*J+1.0)) * (-(qo*qo+qi*qi)*integral_P + 2*qo*qi*integral_0);
      if (J!= 0)
      {
         V_coupled_mm    = isoFac(J-1,1) * (2./(2*J+1)) * ((qo*qo+qi*qi)*integral_M - 2*qo*qi*integral_0);
			V_coupled_mp    = isoFac(J-1,1) * (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*integral_P + qo*qo*integral_M - 2*qo*qi*integral_0);
			V_coupled_pm    = isoFac(J+1,1) * (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*integral_M + qo*qo*integral_P - 2*qo*qi*integral_0);
      }
   }
   // Fill the output array with the correct values
   output[0] = V_uncoupled_S0;
   output[1] = V_uncoupled_S1;
   output[2] = V_coupled_pp;
   output[3] = V_coupled_mm;
   output[4] = V_coupled_pm;
   output[5] = V_coupled_mp;
}
