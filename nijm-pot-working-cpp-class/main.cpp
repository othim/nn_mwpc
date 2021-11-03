#include "nijmegen.h"


int main(int argc, char** argv)
{
   

    double V_arr[6];
    double qo = 1;
    double qi = 1;
    bool coupled = false;
    int S = 1;
    int J = 0;
    int T = 0;

    potential_model* potential = potential_model::fetch_potential_ptr("nijmegen","np"); 
    potential->V(qi,qo,coupled,S,J,T,&V_arr[0]);
    
    for (int i = 0; i < 6; i++)
    {
        std::cout << V_arr[i] << " "; 
    }
    std::cout << std::endl;
    return 0;
}