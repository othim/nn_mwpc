#include "Term.h"
#include "Constants.h"
#include <iostream>


Term::Term(std::string name)
{
    if (name == "OPEP")
    {
        term_name_ = name;
        spin_structure_ = "tensor";
        is_lec_ = false;
        my_v_alpha = &Term::v_alpha_OPEP; // Make my_v_alpha point to the correct function for OPEP
        lecs_in_term_.push_back("gA2");
    } else if (name == "C1S0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        is_lec_ = true;
        LS_term_lec_.J = 0;
        LS_term_lec_.Li = 0;
        LS_term_lec_.Lo = 0; 
        LS_term_lec_.S = 1;
        my_LEC_term= &Term::mom_C1S0;
        lecs_in_term_.push_back("C1S0");
    } else
    {
        std::cerr << "Invalid input to Term constructor: term_name does not match any known term" << std::endl;
    }
}

Term::~Term()
{
    
}

/* Some getters */

bool Term::is_lec()
{
    return is_lec_;
}

std::string Term::get_term_name()
{
    return term_name_;
}

std::string Term::get_spin_structure()
{
    return spin_structure_;
}

bool Term::get_isovector()
{
    return isovector_;
}

LS_term Term::get_LS_term()
{
    return LS_term_lec_;
}

std::vector<std::string> Term::get_lecs_in_term()
{
    return lecs_in_term_;
}


double Term::get_v_alpha(double qi, double qo, double z,std::unordered_map<std::string,double> LECs)
{
    return my_v_alpha(qi,qo,z,LECs);
}

double Term::get_LEC_element(double qi, double qo,std::unordered_map<std::string,double> LECs)
{
    return my_LEC_term(qi,qo,LECs);
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
double Term::v_alpha_OPEP(double qi, double qo, double z,std::unordered_map<std::string,double> LECs)
{
    double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	return -(constants::gA*constants::gA/(4.0*constants::fpi*constants::fpi))*(1.0/(q2+constants::mpi*constants::mpi));
}

// C1S0
double Term::mom_C1S0(double qi, double qo, std::unordered_map<std::string,double> LECs)
{
    return LECs["C1S0"]; // If momentum dependent, multiply here
}