#include "pot_nn_mwpc.h"
#include "Constants.h"
#include <iostream>
#include <cmath>

// This struct is to define what set of parameters the function f_int want
struct my_f_params { double qi; double qo; int J; int l; Term* term; Potential_mwpc* this_pot;};


// Constructor
Potential_mwpc::Potential_mwpc(std::vector<std::string> terms)
{
   // Construct terms and append them to terms_in_pot
   for (std::size_t i = 0; i < terms.size(); i++)
   {
      terms_in_pot_.push_back(Term(terms[i]));
      std::cout << "Added " << terms[i] << " to terms_in_pot_" << std::endl;
   }

   // Make grid for angular integration
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   int_ang_ = gsl_integration_fixed_alloc(T, N_GLI_PWA_, -1.0, 1.0, 0, 0);
   
   double* d = gsl_integration_fixed_nodes(int_ang_);
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
   gsl_integration_fixed_free(int_ang_);
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
	return -(constants::gA*constants::gA/(4.0*constants::fpi*constants::fpi))*(1.0/(q2+constants::mpi*constants::mpi));
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


double Potential_mwpc::calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz)
{
   // TODO implement
   return 0;
}

double f_int_helper(double z, void* p)
{
   // Decode the void* to the parameters
   struct my_f_params * params = (struct my_f_params *)p;
   double qi = (params->qi);
   double qo = (params->qo);
   int J = (params->J);
   int l = (params->l);
   
   return (params->term)->get_v_alpha(qi,qo,z,params->this_pot->LECs_)*gsl_pow_int(z,l)*gsl_sf_legendre_Pl(J, z);
}

double Potential_mwpc::compute_A_integral(double qi, double qo, int J,int l,Term* term)
{
   // Define function to integrate
   gsl_function F;
   F.function = &f_int_helper;
   struct my_f_params params = {qi, qo, J, l, term, this};
   F.params = &params;

   double result;
   // Perform the integration from -1 to 1
   gsl_integration_fixed(&F,&result,int_ang_);

   return result*M_PI;
}

void Potential_mwpc::calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr)
{
   for (std::size_t i = 0; i < terms_in_pot_.size(); i++)
   {  
      if (!terms_in_pot_[i].is_lec())
      {
         // Calculate the A_x integrals (don't calculate all)
         double A_0 = compute_A_integral(qi,qo,J,0,&terms_in_pot_[i]);
         double A_1 = 0;
         double A_P = compute_A_integral(qi,qo,J+1,0,&terms_in_pot_[i]);
         double A_M = 0;

         if (J!=0) {
            A_M = compute_A_integral(qi,qo,J-1,0,&terms_in_pot_[i]);
         }
         if (!coupled) {
            A_1 = compute_A_integral(qi,qo,J,1,&terms_in_pot_[i]);
         }

         // Call pwa with the correct spin structure from this term
         // This will fill up the array V_arr with the correct potential elements
         pwa(qi,qo,coupled,J,A_M,A_P,A_0,A_1,terms_in_pot_[i].get_spin_structure(),terms_in_pot_[i].get_isovector(),V_arr);

      } else 
      {
         // Define some variables and initialize V_arr to all zeros.
         double V_uncoupled_S0 = 0;
         double V_uncoupled_S1 = 0;
         double V_coupled_mm   = 0;
         double V_coupled_pm   = 0;
         double V_coupled_mp   = 0;
         double V_coupled_pp   = 0;

         // Get JLS from the Term to know where to apend the lec value 
         LS_term LS_term = terms_in_pot_[i].get_LS_term();
         // Find which ouput to append (if any)
         if (J==LS_term.J)
         {
            if (coupled)
            {
               if (LS_term.Li == J-1 && LS_term.Lo == J-1) // --
               {
                  V_coupled_mm = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J+1 && LS_term.Lo == J+1) // ++
               {
                  V_coupled_pp = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // -+
               {
                  V_coupled_mp = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // +-
               {
                  V_coupled_pm = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               }
            } else
            {
               if (LS_term.S == 0) // S0
               {
                  V_uncoupled_S0 = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.S == 1) // S1
               {
                  V_uncoupled_S1 = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               }
            }
         }
         V_arr[0] = V_uncoupled_S0;
         V_arr[1] = V_uncoupled_S1;
         V_arr[2] = V_coupled_pp;
         V_arr[3] = V_coupled_mm;
         V_arr[4] = V_coupled_pm;
         V_arr[5] = V_coupled_mp;
      }
   }
}

void Potential_mwpc::pwa(double qi,double qo, bool coupled, int J,double A_M,double A_P,double A_0,double A_1,std::string spin_struct,bool isovector,double* V_arr)
{
   // Define some variables
   double V_uncoupled_S0 = 0;
   double V_uncoupled_S1 = 0;
   double V_coupled_mm   = 0;
   double V_coupled_pm   = 0;
   double V_coupled_mp   = 0;
   double V_coupled_pp   = 0;

   if (spin_struct == "tensor")
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
            V_coupled_mp = (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_P + qo*qo*A_M - 2*qo*qi*A_0);
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