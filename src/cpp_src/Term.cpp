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
    } else if (name == "C3P0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_C3P0;
        lecs_in_term_.push_back("C3P0");
    } else if (name == "C3P2")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_C3P2;
        lecs_in_term_.push_back("C3P2");
    } else if (name == "C3D2")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 2;
        LS_well_def_pw_.Lo = 2;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_C3D2;
        lecs_in_term_.push_back("C3D2");
    } else if (name == "Yamaguchi_1S0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0;
        LS_well_def_pw_.S = 0;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_Yamaguchi_1S0;
        lecs_in_term_.push_back("Yamaguchi_1S0");
        params_in_term_.push_back("Yamaguchi_beta");
    }
    else if (name == "Yamaguchi_3S1") 
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J  = 1;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0;
        LS_well_def_pw_.S  = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_Yamaguchi_3S1;
        params_in_term_.push_back("Yamaguchi_lambda_00");
    }
    else if (name == "Yamaguchi_3S-D1") 
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J  = 1;
        LS_well_def_pw_.Li = 2;
        LS_well_def_pw_.Lo = 0;
        LS_well_def_pw_.S  = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_Yamaguchi_3S_D1;
        params_in_term_.push_back("Yamaguchi_lambda_01");
    }
    else if (name == "Yamaguchi_3D-S1") 
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J  = 1;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 2;
        LS_well_def_pw_.S  = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_Yamaguchi_3D_S1;
        params_in_term_.push_back("Yamaguchi_lambda_10");
    }
    else if (name == "Yamaguchi_3D1") 
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J  = 1;
        LS_well_def_pw_.Li = 2;
        LS_well_def_pw_.Lo = 2;
        LS_well_def_pw_.S  = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_Yamaguchi_3D1;
        params_in_term_.push_back("Yamaguchi_lambda_11");
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

std::vector<std::string> Term::get_params_in_term()
{
    return params_in_term_;
}

std::vector<double> Term::get_v_alpha(double qi, double qo, double* z,unsigned int z_len,
        std::unordered_map<std::string,double> LECs,
        std::unordered_map<std::string,double>& params)
{
    return my_v_alpha(qi,qo,z,z_len,LECs,params);
}

double Term::get_v_alpha_well_def_pw(double qi, double qo,
        std::unordered_map<std::string,double> LECs,
        std::unordered_map<std::string,double>& params)
{
    return my_v_alpha_well_def_pw(qi,qo,LECs,params);
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
std::vector<double> Term::v_alpha_OPEP(double qi, double qo, double* z,unsigned int z_len, 
        std::unordered_map<std::string,double> LECs,
        std::unordered_map<std::string,double>& params)
{
    double lec = LECs["gA2"];
    std::vector<double> tmp(z_len);
    double q2;
    
    for (int i = 0; i < (int)z_len; i++)
    {
        q2 = qi*qi + qo*qo - 2*qi*qo*z[i];
        tmp[i] = -(lec/(4.0*constants::fpi*constants::fpi))*(1.0/(q2+constants::mpi*constants::mpi));
    }
	return tmp;
}

/*
 * Leading order contact terms in chiral EFT
 */

double Term::mom_C1S0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["C1S0"]; // If momentum dependent, multiply here
}

double Term::mom_C3S1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["C3S1"]; // If momentum dependent, multiply here
}

double Term::mom_C3P0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["C3P0"]*qi*qo; 
}

double Term::mom_C3P2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["C3P2"]*qi*qo; 
}

double Term::mom_C3D2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["C3D2"]*qi*qi*qo*qo; // If momentum dependent, multiply here
}

/*
 * Here is the function that defines the 1S0 separable Yamaguchi potential
 */

double Term::mom_Yamaguchi_1S0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    //double beta   = 20.0; // MeV
    double beta   = params["Yamaguchi_beta"];
    double lambda = LECs["Yamaguchi_1S0"];

    return lambda*(1.0/(qi*qi + beta*beta))*(1.0/(qo*qo + beta*beta));

}

/*
 *
 * Here are functions defining the 3S_1 - 3D_1 coupled separable 
 * Yamaguchi potential.
 *
 */

double g2(double p, double beta)
{
    return 1/((p*p + beta*beta)*(p*p + beta*beta));
}

double Term::mom_Yamaguchi_3S1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = params["Yamaguchi_lambda_00"];
    
    return lambda*g2(qo,beta)*g2(qi,beta);
}

double Term::mom_Yamaguchi_3S_D1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = params["Yamaguchi_lambda_01"];
    
    return lambda*g2(qo,beta)*g2(qi,beta)*qi*qi;
}

double Term::mom_Yamaguchi_3D_S1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = params["Yamaguchi_lambda_10"];
    
    return lambda*g2(qo,beta)*g2(qi,beta)*qo*qo;
}

double Term::mom_Yamaguchi_3D1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = params["Yamaguchi_lambda_11"];
    
    return lambda*g2(qo,beta)*g2(qi,beta)*qo*qo*qi*qi;
}
