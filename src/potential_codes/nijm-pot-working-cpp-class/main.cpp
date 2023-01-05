#include "nijmegen.h"

/*
// External C-function to call fortran script 
extern "C" {
    void nijmegen_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}
*/

int main(int argc, char** argv)
{
   

    double V_arr[6];
    double qo = 1;
    double qi = 1;
    bool coupled = false;
    int S = 1;
    int J = 0;
    int T = 0;
    int Tz = 0; // np

    potential_model* potential = potential_model::fetch_potential_ptr("nijmegen","np"); 
    potential->V(qi,qo,coupled,S,J,T,&V_arr[0]);
    
    for (int i = 0; i < 6; i++)
    {
        std::cout << V_arr[i] << " "; 
    }
    std::cout << std::endl;
    
    // Make potential pointer to the potential
    int coup = 0;
    nijmegen_fort_interface(&qi, &qo, &coup, &S, &J, &T, &Tz, &V_arr[0]); 
     
    std::cout << "my_f" << std::endl;

    for (int i = 0; i < 6; i++)
    {
        std::cout << V_arr[i] << " "; 
    }
    std::cout << std::endl;
    
    return 0;
}
