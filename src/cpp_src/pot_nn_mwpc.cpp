#include "pot_nn_mwpc.h"
#include "Constants.h"
#include "Term.h"
#include <iostream>
#include <cmath>

// Constructor
Potential_mwpc::Potential_mwpc(std::vector<string> terms)
{
   // Construct terms and append them to terms_in_pot
   for (std::size_t i = 0; i < terms.size(); i++)
   {
      terms_in_pot.push_back(Term(terms[i]));
      std::cout << "Added " << terms[i] << " to terms_in_pot" << std::endl;
   }

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

double Potential_mwpc::calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz)
{

}

// This struct is to define what set of parameters the function f_int want
struct my_f_params { double qi; double qo; int J; int l; Term* term};

double Potential_mwpc::f_int_helper(double z, void* p)
{
   // Decode the void* to the parameters
   struct my_f_params * params = (struct my_f_params *)p;
   double qi = (params->qi);
   double qo = (params->qo);
   int J = (params->J);
   int l = (params->l);
   
   return *(params->term).get_v_alpha(qi,qo,z)*gsl_pow_int(z,l)*gsl_sf_legendre_Pl(J, z);
}

double Potential_mwpc::compute_A_integral(double qi, double qo, int J,int l,Term* term)
{
   // Define function to integrate
   gsl_function F;
   F.function = &f_int_helper;
   struct my_f_params params = {qi, qo, J, l,term};
   F.params = &params;

   double result;
   // Perform the integration from -1 to 1
   gsl_integration_fixed(&F,&result,int_ang);

   return result*M_PI;
}

void Potential_mwpc::calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr)
{
   for (std::size_t i = 0; i < terms.size(); i++)
   {  
      if (!term.is_lec())
      {
         // Calculate the A_x integrals (don't calculate all)
         A_0 = compute_A_integral(qi,qo,J,0,&terms_in_pot[i]);
         A_P = compute_A_integral(qi,qo,J,0,&terms_in_pot[i]);
         A_M = 0;
         if (J!=0) {
            compute_A_integral(qi,qo,J,0,&terms_in_pot[i]);
         }

         // call pwa with the correct spin structure from this term
      } else 
      {
         // Get JLS from the Term to know where to apend the lec value 

         // Find which ouput to append (if any)

         // Call Term function get_element() to get the 
      }
   }
}

void Potential_mwpc::pwa(double qi,double qo, bool coupled, int J,A_m,A_p,A_0,A_1,std::string spin_struct,bool isovector,double* V_arr)
{
   // Define some variables
   double V_uncoupled_S0 = 0;
   double V_uncoupled_S1 = 0;
   double V_coupled_mm   = 0;
   double V_coupled_pm   = 0;
   double V_coupled_mp   = 0;
   double V_coupled_pp   = 0;

   if (spin_struct == 'tensor')
   {
      // Check which elements that are non-zero by checking if the Channel is
      // coupled or not
      if (!coupled)
      {
         // OPEP uncoupled interactions
         V_uncoupled_S0    = 2 * (-(qo*qo+qi*qi)*A_0 + 2*qo*qi*A_1);
         if (J!=0) {
            V_uncoupled_S1 = 2 * ((qo*qi+qi*qi)*A_0 - 2*qo*qi*(1.0/(2.0*J+1.0))*(J*A_P + (J+1)*A_M));
         }
      } else
      {
         V_coupled_pp = (2.0/(2.0*J+1.0)) * (-(qo*qo+qi*qi)*A_P + 2*qo*qi*A_0);
         if (J!= 0)
         {
            V_coupled_mm = (2./(2*J+1)) * ((qo*qo+qi*qi)*A_M - 2*qo*qi*A_0);
            V_coupled_mp = 4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_P + qo*qo*A_M - 2*qo*qi*A_0);
            V_coupled_pm = (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_M + qo*qo*A_P - 2*qo*qi*A_0);
         }
      }

      // Add isospin factor if the term is an isovector
      if (isovector)
      {
         V_uncoupled_S0 *= isoFac(J,0);
         V_uncoupled_S1 *= isoFac(J,1);
         V_coupled_mm   *= isoFac(J-1,1);
         V_coupled_pm   *= isoFac(J+1,1);
         V_coupled_mp   *= isoFac(J-1,1);
         V_coupled_pp   *= isoFac(J+1,1);

      }
    } else
    {
       std::cerr << "Unknown <spin_struct> passed to function <pwa> in Potential_mwpc object." << std::endl;
    }
      // Fill the output array with the correct values
      V_arr[0] = V_uncoupled_S0;
      V_arr[1] = V_uncoupled_S1;
      V_arr[2] = V_coupled_pp;
      V_arr[3] = V_coupled_mm;
      V_arr[4] = V_coupled_pm;
      V_arr[5] = V_coupled_mp;
}