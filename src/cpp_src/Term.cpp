#include "Term.h"
#include "Constants.h"
#include <iostream>


Term::Term(std::string name)
{
    if (name == 'OPEP')
    {
        term_name = name;
        spin_structure = 'tensor';
        is_lec = false;
        my_v_alpha = &Term::v_alpha_OPEP; // Make my_v_alpha point to the correct function for OPEP
    } else if (name == 'C1S0')
    {
        term_name = name;
        spin_structure = 'none';
        is_lec = true;
        LS_term_lec.J = 0;
        LS_term_lec.L = 0;
        LS_term_lec.S = 1;
        my_v_alpha= &Term::mom_C1S0;
    } else
    {
        std::cerr << "Invalid input to Term constructor: term_name does not match any known term" << std::endl;
    }
}

Term::~Term()
{
    
}

/* Some getters */

std::string Term::get_term_name()
{
    return term_name;
}

std:string Term::get_spin_structure()
{
    return spin_structure;
}

std::vector<string> Term::get_lecs_in_term()
{
    return lecs_in_term;
}


double get_v_alpha(double qi, double qo, double z,std::unordered_map<string,double> LECs)
{
    return my_v_alpha(qi,qo,z,LECs);
}

/*
    These functions are the V_\alpha in Erkelenz for the 
    specific potential terms. Naming convention
    v_alpha_<term_name>.

    In general v_alpha_<term_name> depend on some LECs specified in the lecs_in_term variable. 
    To sucessfully call v_alpha_<term_name> std::unordered_map<string,double> LECs MUST 
    contain the LECs neseccary to be able to compute that v_alpha, else the program 
    will crash.
*/

// OPEP
static double Term::v_alpha_OPEP(double qi, double qo, double z,std::unordered_map<string,double> LECs)
{
    double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	return -(gA*gA/(4.0*fpi*fpi))*(1.0/(q2+mpi*mpi));
}

// C1S0
static double Term::mom_C1S0(double qi, double qo, double z, std::unordered_map<string,double> LECs)
{
    return LECs['C1S0']; // If momentum dependent, multiply here
}