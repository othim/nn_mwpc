#include "Term.h"
#include "Constants.h"
#include <iostream>


Term::Term(std::string name)
{
    if (name == "OPEP")
    {
        term_name_ = name;
        spin_structure_ = "tensor";
        well_def_pw_ = false;
        my_v_alpha = &Term::v_alpha_OPEP; // Make my_v_alpha point to the correct function for OPEP
        lecs_in_term_.push_back("gA2");
        isovector_ = true;
    } else if (name == "C1S0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0; 
        LS_well_def_pw_.S = 0;
        my_v_alpha_well_def_pw = &Term::mom_C1S0;
        isovector_ = false;
        lecs_in_term_.push_back("C1S0");
    } else if (name == "C3S1")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0; 
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw= &Term::mom_C3S1;
        lecs_in_term_.push_back("C3S1");
    } else 
    {
        std::cerr << "Invalid input to Term constructor: term_name does not match any known term" << std::endl;
    }
}

Term::~Term()
{
    
}

/* Some getters */

bool Term::well_def_pw()
{
    return well_def_pw_;
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
    return LS_well_def_pw_;
}

std::vector<std::string> Term::get_lecs_in_term()
{
    return lecs_in_term_;
}


std::vector<double> Term::get_v_alpha(double qi, double qo, double* z,unsigned int z_len,std::unordered_map<std::string,double> LECs)
{
    return my_v_alpha(qi,qo,z,z_len,LECs);
}

double Term::get_v_alpha_well_def_pw(double qi, double qo,std::unordered_map<std::string,double> LECs)
{
    return my_v_alpha_well_def_pw(qi,qo,LECs);
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
std::vector<double> Term::v_alpha_OPEP(double qi, double qo, double* z,unsigned int z_len, std::unordered_map<std::string,double> LECs)
{
    double lec = LECs["gA2"];
    std::vector<double> tmp(z_len);
    double q2;
    
    for (int i = 0; i < z_len; i++)
    {
        q2 = qi*qi + qo*qo - 2*qi*qo*z[i];
        tmp[i] = -(lec/(4.0*constants::fpi*constants::fpi))*(1.0/(q2+constants::mpi*constants::mpi));
    }
	return tmp;
}

// C1S0
double Term::mom_C1S0(double qi, double qo, std::unordered_map<std::string,double> LECs)
{
    return LECs["C1S0"]; // If momentum dependent, multiply here
}
// C3S1
double Term::mom_C3S1(double qi, double qo, std::unordered_map<std::string,double> LECs)
{
    return LECs["C3S1"]; // If momentum dependent, multiply here
}