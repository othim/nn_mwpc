#include "pot_nn_mwpc.h"

// Constructor
Potential_mwpc::Potential_mwpc()
{
   // Init everything
}

/* Some helper functions to opep_get_el
******************************************
*/

double pot_OPEP_mom(double qi,double qo, double z)
{
}

// Do with gsl
double find_root(int d, double x)
{
}

double ang_Integral(double qi,double qo,int J,int l)
{
}

int isoFac(int L,int S)
{
}

void Potential_mwpc::opep_get_el(double qi,double qo, bool coupled, int J, double* output)
{

   // Precompute some integrals

   // Check which elements that are non-zero by checking if the Channel is
   // coupled or not

   // Make the output point to the values

}
