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
        lecs_in_term_.push_back("Yamaguchi_lambda_00");
        params_in_term_.push_back("Yamaguchi_beta");
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
        lecs_in_term_.push_back("Yamaguchi_lambda_01");
        params_in_term_.push_back("Yamaguchi_beta");
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
        lecs_in_term_.push_back("Yamaguchi_lambda_10");
        params_in_term_.push_back("Yamaguchi_beta");
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
        lecs_in_term_.push_back("Yamaguchi_lambda_11");
        params_in_term_.push_back("Yamaguchi_beta");
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
    double lambda = LECs["Yamaguchi_lambda_00"];
    
    return lambda*g2(qo,beta)*g2(qi,beta);
}

double Term::mom_Yamaguchi_3S_D1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = LECs["Yamaguchi_lambda_01"];

    double pot = lambda*g2(qo,beta)*g2(qi,beta)*qi*qi;

    //std::cout << "beta=" << beta << ", lambda=" << lambda << std::endl;
    //std::cout << "SD," << "qi=" << qi << ", qo=" << qo << ", pot=" << pot << std::endl;
    
    return pot;
}

double Term::mom_Yamaguchi_3D_S1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = LECs["Yamaguchi_lambda_10"];
    double pot = lambda*g2(qo,beta)*g2(qi,beta)*qo*qo;
    //std::cout << "beta=" << beta << ", lambda=" << lambda << std::endl;
    //std::cout << "DS," << "qi=" << qi << ", qo=" << qo << ", pot=" << pot << std::endl;
    return pot;
}

double Term::mom_Yamaguchi_3D1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double beta   = params["Yamaguchi_beta"];
    double lambda = LECs["Yamaguchi_lambda_11"];
    
    return lambda*g2(qo,beta)*g2(qi,beta)*qo*qo*qi*qi;
}

/*
 * nu=2 (N2LO) Chiral two-pion exchange contributions in dimensional
 * regularization. As eq. 4.9-4.12 in M&E Phys. Rept 503 (2011)
 */

double w_f(double q, double mpi)
{
    return std::sqrt(4*mpi*mpi+q*q);
}


double L_DR(double q, double mpi)
{
    double w = w_f(q,mpi);

    return (w/q)*std::log((w+q)/(2*mpi));
}

static double V_T_2pi_nu_2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double q = qo-qi; // momentum transfer

    return (-3.0*std::pow(gA,4)*L_DR(q,mpi))/(64.0*M_PI*M_PI*std::pow(fpi,4));
}

static double V_S_2pi_nu_2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double q = qo-qi; // momentum transfer
    return -q*q*V_T_2pi_nu_2(qi,qo,LECs,params);
}

static double W_C_2pi_nu_2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double q = qo-qi; // momentum transfer

    double gA2 = gA*gA;
    double gA4 = std::pow(gA,4);
    
    double w = w_f(q,mpi);
    
    double fac = 4.0*mpi*mpi*(5.0*gA4-4.0*gA2-1)+q*q*(23.0*gA4-10*gA2-1)+
        (48.0*gA4*std::pow(mpi,4))/(w*w);

    return ((-L_DR(q,mpi))/(384.0*M_PI*M_PI*std::pow(fpi,4)))*fac;

}


/*
 * nu=3 (N3LO) chiral two-pion exchange contributions in dimensional
 * regularization. As eq. 4.13-4.20 in M&E Phys. Rep. 503 (2011).
 */

static double w_tilde_f(double q, double mpi)
{
    return std::sqrt(2.0*mpi*mpi+q*q);
}

static double A_DR(double q, double mpi)
{
    return (1.0/(2.0*q))*std::atan(q/(2.0*mpi));
}

static double V_C_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];
    double c1  = LECs["c1"];
    double c3  = LECs["c3"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA2 = gA*gA;
    
    double w = w_f(q,mpi);
    double w_t = w_tilde_f(q,mpi);

    
    double tmp1 = (gA2*std::pow(mpi,5))/(16.0*mN*w*w);

    double tmp2 = 2.0*mpi*mpi*(2.0*c1-c3) - q*q*(c3 + (3.0*gA2)/(16.0*mN));
    
    double tmp3 = (3.0*gA2)/(16.0*M_PI*std::pow(fpi,4));

    return tmp3*(tmp1-tmp2*w_t*w_t*A_DR(q,mpi));
}   


static double W_C_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA2 = gA*gA;
    
    double w = w_f(q,mpi);
    double w_t = w_tilde_f(q,mpi);

    double tmp1 = (3.0*gA2*std::pow(mpi,5))/(w*w);

    double tmp2 = 4.0*mpi*mpi+2.0*q*q-gA2*(4.0*mpi*mpi+3.0*q*q);

    double tmp3 = (gA2)/(128.0*M_PI*mN*std::pow(fpi,4));

    return tmp3*(tmp1-tmp2*w_t*w_t*A_DR(q,mpi));
}

static double V_T_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA4 = std::pow(gA,4);
    
    double w_t = w_tilde_f(q,mpi);

    return (9.0*gA4*w_t*w_t*A_DR(q,mpi))/(512.0*M_PI*mN*std::pow(fpi,4));
}

static double V_S_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double q = qo-qi;

    return -q*q*V_T_2pi_nu_3(qi,qo,LECs,params);
}

static double W_T_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];
    double c4  = params["c4"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA2 = std::pow(gA,2);
    
    double w = w_f(q,mpi);

    
    double tmp1 = (c4+1.0/(4.0*mN))*w*w - (gA2/(8.0*mN))*(10.0*mpi*mpi+3.0*q*q);

    double tmp2 = (-gA2*A_DR(q,mpi))/(32.0*M_PI*std::pow(fpi,4));
    
    return tmp2*tmp1;
}

static double W_S_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double q = qo-qi;

    return -q*q*W_T_2pi_nu_3(qi,qo,LECs,params);
}

static double V_LS_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA4 = std::pow(gA,4);
    
    double w_t = w_tilde_f(q,mpi);

    return (3.0*gA4*w_t*w_t*A_DR(q,mpi))/(32.0*M_PI*mN*std::pow(fpi,4));
}

static double W_LS_2pi_nu_3(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = 2.0*constants::Mp*constants::Mn/(constants::Mp+constants::Mn);
    double q = qo-qi; // momentum transfer

    double gA2 = std::pow(gA,2);
    
    double w = w_f(q,mpi);

    return (gA2*(1-gA2)*w*w*A_DR(q,mpi))/(32.0*M_PI*mN*std::pow(fpi,4));
}


