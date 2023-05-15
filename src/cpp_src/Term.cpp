#include "Term.h"
#include "Constants.h"
#include <iostream>


Term::Term(std::string name)
{
    if (name == "W_T_1pi_nu_0")
    {
        term_name_ = name;
        spin_structure_ = "T";
        well_def_pw_ = false;
        my_v_alpha = &Term::v_alpha_OPEP; // Make my_v_alpha point to the correct function for OPEP
        params_in_term_.push_back("gA");
        isovector_ = true;
    } else if (name == "V_T_2pi_nu_2")
    {
        term_name_ = name;
        spin_structure_ = "T";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_T_2pi_nu_2; 
        params_in_term_.push_back("gA");
        isovector_ = false;
    } else if (name == "V_S_2pi_nu_2")
    {
        term_name_ = name;
        spin_structure_ = "S";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_S_2pi_nu_2; 
        params_in_term_.push_back("gA");
        isovector_ = false;
    } else if (name == "W_C_2pi_nu_2")
    {
        term_name_ = name;
        spin_structure_ = "C";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_C_2pi_nu_2; 
        params_in_term_.push_back("gA");
        isovector_ = true;
    } else if (name == "W_T_1pi_nu_2")
    {
        term_name_ = name;
        spin_structure_ = "T";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_T_1pi_nu_2; 
        params_in_term_.push_back("gA");
        params_in_term_.push_back("d18");
        isovector_ = true;
    } else if (name == "V_C_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "C";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_C_2pi_nu_3; 
        params_in_term_.push_back("gA");
        lecs_in_term_.push_back("c1");
        lecs_in_term_.push_back("c3");
        isovector_ = false;
    } else if (name == "W_C_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "C";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_C_2pi_nu_3; 
        params_in_term_.push_back("gA");
        isovector_ = true;
    } else if (name == "V_T_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "T";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_T_2pi_nu_3; 
        params_in_term_.push_back("gA");
        isovector_ = false;
    } else if (name == "V_S_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "S";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_S_2pi_nu_3; 
        params_in_term_.push_back("gA");
        isovector_ = false;
    } else if (name == "W_T_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "T";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_T_2pi_nu_3; 
        params_in_term_.push_back("gA");
        lecs_in_term_.push_back("c4");
        isovector_ = true;
    } else if (name == "W_S_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "S";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_S_2pi_nu_3; 
        params_in_term_.push_back("gA");
        lecs_in_term_.push_back("c4");
        isovector_ = true;
    } else if (name == "V_LS_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "LS";
        well_def_pw_ = false;
        my_v_alpha = &Term::V_LS_2pi_nu_3; 
        params_in_term_.push_back("gA");
        isovector_ = false;
    } else if (name == "W_LS_2pi_nu_3")
    {
        term_name_ = name;
        spin_structure_ = "LS";
        well_def_pw_ = false;
        my_v_alpha = &Term::W_LS_2pi_nu_3; 
        params_in_term_.push_back("gA");
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
        
    } else if (name == "D1S0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0; 
        LS_well_def_pw_.S = 0;
        isovector_ = false;
        my_v_alpha_well_def_pw= &Term::mom_D1S0;
        lecs_in_term_.push_back("D1S0");
        
    } else if (name == "D3S1")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0; 
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw= &Term::mom_D3S1;
        lecs_in_term_.push_back("D3S1");
        
    } else if (name == "D1P1")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1; 
        LS_well_def_pw_.S = 0;
        isovector_ = false;
        my_v_alpha_well_def_pw= &Term::mom_D1P1;
        lecs_in_term_.push_back("D1P1");
        
    } else if (name == "D3P0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_D3P0;
        lecs_in_term_.push_back("D3P0");
    } else if (name == "D3P1")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_D3P1;
        lecs_in_term_.push_back("D3P1");
    } else if (name == "D3P2")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_D3P2;
        lecs_in_term_.push_back("D3P2");
    } else if (name == "D_SD")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 2;
        LS_well_def_pw_.Lo = 0;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_D_SD;
        lecs_in_term_.push_back("D_SD");
    } else if (name == "D_DS")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 1;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 2;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_D_DS;
        lecs_in_term_.push_back("D_SD"); // NOTE D_SD! 
    } else if (name == "E1S0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 0;
        LS_well_def_pw_.Lo = 0;
        LS_well_def_pw_.S = 0;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_E1S0;
        lecs_in_term_.push_back("E1S0");
    } else if (name == "E3P0")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 0;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_E3P0;
        lecs_in_term_.push_back("E3P0");
    } else if (name == "E3P2")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_E3P2;
        lecs_in_term_.push_back("E3P2");
    } else if (name == "E_PF")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 3;
        LS_well_def_pw_.Lo = 1;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_E_PF;
        lecs_in_term_.push_back("E_PF");
    } else if (name == "E_FP")
    {
        term_name_ = name;
        spin_structure_ = "none";
        well_def_pw_ = true;
        LS_well_def_pw_.J = 2;
        LS_well_def_pw_.Li = 1;
        LS_well_def_pw_.Lo = 3;
        LS_well_def_pw_.S = 1;
        isovector_ = false;
        my_v_alpha_well_def_pw = &Term::mom_E_FP;
        lecs_in_term_.push_back("E_PF");
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
        std::cerr << "Invalid input to Term constructor: term_name=" << name << 
            ", does not match any known term" << std::endl;
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

std::vector<double> Term::get_v_alpha(double qi, double qo, double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    return my_v_alpha(qi,qo,z,z_len,LECs,params,chn,loop_reg,lam_SFR);
}

double Term::get_v_alpha_well_def_pw(double qi, double qo,
        std::unordered_map<std::string,double>& LECs,
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
std::vector<double> Term::v_alpha_OPEP(double qi, double qo, double* z, int z_len, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double lec = params["gA"];
    std::vector<double> tmp(z_len);
    double q2;
    
    for (int i = 0; i < (int)z_len; i++)
    {
        q2 = qi*qi + qo*qo - 2*qi*qo*z[i];
        tmp[i] = -(lec*lec/(4.0*constants::fpi*constants::fpi))*(1.0/
                (q2+constants::mpi*constants::mpi));
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


double Term::mom_D1S0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D1S0"]*(qi*qi+qo*qo); 
}

double Term::mom_D3S1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D3S1"]*(qi*qi+qo*qo); 
}

double Term::mom_D1P1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D1P1"]*qi*qo; 
}

double Term::mom_D3P0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["D3P0"]*qi*qo; 
}

double Term::mom_D3P1(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)
{
    return LECs["D3P1"]*qi*qo; 
}

double Term::mom_D3P2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D3P2"]*qi*qo; 
}

double Term::mom_D_SD(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D_SD"]*qi*qi; 
}

double Term::mom_D_DS(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["D_SD"]*qo*qo; 
}

double Term::mom_E1S0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["E1S0"]*qo*qo*qi*qi; 
}

double Term::mom_E3P0(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["E3P0"]*qo*qi*(qi*qi+qo*qo); 
}

double Term::mom_E3P2(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["E3P2"]*qo*qi*(qi*qi+qo*qo); 
}

double Term::mom_E_PF(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["E_PF"]*qo*qi*(qi*qi); 
}

double Term::mom_E_FP(double qi, double qo, 
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params)

{
    return LECs["E_PF"]*qo*qi*(qo*qo); 
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
double Term::get_q(double qi, double qo, double z)
{
    return std::sqrt(qi*qi + qo*qo - 2*qi*qo*z);
}

double Term::w_f(double q, double mpi)
{
    return std::sqrt(4*mpi*mpi+q*q);
}


double Term::L_DR(double q, double mpi)
{
    double w = w_f(q,mpi);

    return (w/q)*std::log((w+q)/(2*mpi));
}

double Term::L_SFR(double q, double mpi, double lam_SFR)
{
    if (lam_SFR > 2.0*mpi)
    {
        double w = Term::w_f(q,mpi);
        
        double s = std::sqrt(lam_SFR*lam_SFR - 4.0*mpi*mpi);

        double ln_fac = (lam_SFR*lam_SFR*w*w + q*q*s*s + 2.0*lam_SFR*q*w*s)/(
                4.0*mpi*mpi*(lam_SFR*lam_SFR + q*q));

        return (w/(2.0*q))*std::log(ln_fac);
    } else 
    {
        return 0;
    }
}


double Term::A_SFR(double q, double mpi, double lam_SFR)
{
    if (lam_SFR > 2.0*mpi)
    {
        double tmp = (q*(lam_SFR - 2.0*mpi))/(q*q + 2.0*lam_SFR*mpi);

        return (1.0/(2.0*q))*std::atan(tmp);
    } else
    {
        return 0.0;
    }

}

double Term::w_tilde_f(double q, double mpi)
{
    return std::sqrt(2.0*mpi*mpi+q*q);
}

double Term::A_DR(double q, double mpi)
{
    return (1.0/(2.0*q))*std::atan(q/(2.0*mpi));
}


double Term::L_gen(double q, double mpi, std::string loop_reg, double lam_SFR)
{
    if (loop_reg == "DR")
    {
        return Term::L_DR(q,mpi);
    } else if (loop_reg == "SFR")
    {
        return Term::L_SFR(q, mpi, lam_SFR);
    } else
    {   
        std::cout << "Error, loop_reg_ is not 'DR' or 'SFR'" << std::endl;
        return -1.0;
    }
}

double Term::A_gen(double q, double mpi, std::string loop_reg, double lam_SFR)
{
    if (loop_reg == "DR")
    {
        return Term::A_DR(q, mpi);   
    } else if (loop_reg == "SFR")
    {
        return Term::A_SFR(q, mpi, lam_SFR);
    } else
    {   
        std::cout << "Error, loop_reg_ is not 'DR' or 'SFR'" << std::endl;
        return -1.0;
    }

}

//*****************************************************************************

std::vector<double> Term::V_T_2pi_nu_2(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    
    // Compute the function for all angles z = cos <qi,qo>
    std::vector<double> tmp(z_len);
    double q;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        tmp[i] = Term::V_T_2pi_nu_2(q, gA, mpi, fpi, loop_reg, lam_SFR);
    }
	return tmp;

}

double Term::V_T_2pi_nu_2(double q, double gA, double mpi, double fpi,
        std::string loop_reg, double lam_SFR)
{
    return (-3.0*std::pow(gA,4)*Term::L_gen(q,mpi,loop_reg,lam_SFR))/
        (64.0*M_PI*M_PI*std::pow(fpi,4));
}

//*****************************************************************************

std::vector<double> Term::V_S_2pi_nu_2(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    
    // Compute the function for all angles z = cos <qi,qo>
    std::vector<double> tmp(z_len);
    double q;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        tmp[i] = -q*q*Term::V_T_2pi_nu_2(q, gA, mpi, fpi, loop_reg, lam_SFR);
    }
	return tmp;
}

//*****************************************************************************

std::vector<double> Term::W_C_2pi_nu_2(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;

    
    std::vector<double> tmp(z_len);
    double q,w;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        tmp[i] = Term::W_C_2pi_nu_2(q, gA, mpi, fpi, w, loop_reg, lam_SFR);
    }
	return tmp;
}

double Term::W_C_2pi_nu_2(double q, double gA, double mpi, double fpi,
        double w, std::string loop_reg, double lam_SFR)
{
    double gA2 = gA*gA;
    double gA4 = std::pow(gA,4);

    double fac = 4.0*mpi*mpi*(5.0*gA4-4.0*gA2-1)+q*q*(23.0*gA4-10.0*gA2-1.0)+
        (48.0*gA4*std::pow(mpi,4))/(w*w);

    return ((-Term::L_gen(q,mpi,loop_reg, lam_SFR))/
            (384.0*M_PI*M_PI*std::pow(fpi,4)))*fac;
}

std::vector<double> Term::W_T_1pi_nu_2(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA = params["gA"];
    double d18 = params["d18"];
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    std::vector<double> tmp(z_len);
    double q;

    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        tmp[i] = ((gA*d18*mpi*mpi)/(fpi*fpi))*(1.0/(q*q + mpi*mpi));
    }
	return tmp;
}
//*****************************************************************************

/*
 * nu=3 (N3LO) chiral two-pion exchange contributions in dimensional
 * regularization. As eq. 4.13-4.20 in M&E Phys. Rep. 503 (2011).
 */


//*****************************************************************************

std::vector<double> Term::V_C_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    // Get the constants
    double gA  = params["gA"];
    double c1  = LECs["c1"];
    double c3  = LECs["c3"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);
    
    // Compute the function for all angles z = cos <qi,qo>
    std::vector<double> tmp(z_len);
    double q,w,w_t;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        w_t    = Term::w_tilde_f(q,mpi);
        tmp[i] = Term::V_C_2pi_nu_3(q, gA, c1, c3, mpi, fpi, mN, w, w_t, 
                loop_reg, lam_SFR);
    }
	return tmp;
}  

double Term::V_C_2pi_nu_3(double q, double gA, double c1, double c3, 
        double mpi, double fpi, double mN, double w, double w_t,
        std::string loop_reg, double lam_SFR)
{
    double gA2 = gA*gA;

    double tmp1 = (gA2*std::pow(mpi,5))/(16.0*mN*w*w);

    double tmp2 = 2.0*mpi*mpi*(2.0*c1-c3) - q*q*(c3 + (3.0*gA2)/(16.0*mN));
    
    double tmp3 = (3.0*gA2)/(16.0*M_PI*std::pow(fpi,4));

    return tmp3*(tmp1-tmp2*w_t*w_t*Term::A_gen(q,mpi,loop_reg,lam_SFR));

}

//*****************************************************************************

std::vector<double> Term::W_C_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    
    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);

    std::vector<double> tmp(z_len);
    double q,w,w_t;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        w_t    = Term::w_tilde_f(q,mpi);
        tmp[i] = Term::W_C_2pi_nu_3(q, gA, mpi, fpi, mN, w, w_t, loop_reg,
                lam_SFR);
    }
	return tmp;
}

double Term::W_C_2pi_nu_3(double q, double gA, double mpi, double fpi, 
        double mN, double w, double w_t, std::string loop_reg, double lam_SFR)
{
    double gA2 = gA*gA;
    
    double tmp1 = (3.0*gA2*std::pow(mpi,5))/(w*w);

    double tmp2 = 4.0*mpi*mpi+2.0*q*q-gA2*(4.0*mpi*mpi+3.0*q*q);

    double tmp3 = (gA2)/(128.0*M_PI*mN*std::pow(fpi,4));

    return tmp3*(tmp1-tmp2*w_t*w_t*Term::A_gen(q,mpi,loop_reg,lam_SFR));
}

//*****************************************************************************
        
std::vector<double> Term::V_T_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);

    
    std::vector<double> tmp(z_len);
    double q,w,w_t;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        w_t    = Term::w_tilde_f(q,mpi);
        tmp[i] = Term::V_T_2pi_nu_3(q, gA, mpi, fpi, mN, w, w_t, loop_reg,
                lam_SFR);
    }
	return tmp;
}

double Term::V_T_2pi_nu_3(double q, double gA, double mpi, double fpi, 
        double mN, double w, double w_t, std::string loop_reg, double lam_SFR)
{
    double gA4 = std::pow(gA,4);
    return (9.0*gA4*w_t*w_t*Term::A_gen(q,mpi,loop_reg,lam_SFR))/
        (512.0*M_PI*mN*std::pow(fpi,4));
}

//*****************************************************************************

std::vector<double> Term::V_S_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);

    
    std::vector<double> tmp(z_len);
    double q,w,w_t;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        w_t    = Term::w_tilde_f(q,mpi);
        tmp[i] = -q*q*Term::V_T_2pi_nu_3(q, gA, mpi, fpi, mN, w, w_t, loop_reg,
                lam_SFR);
    }
	return tmp;
}

//*****************************************************************************

std::vector<double> Term::W_T_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    double c4  = params["c4"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);
    

    std::vector<double> tmp(z_len);
    double q,w;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        tmp[i] = Term::W_T_2pi_nu_3(q, c4, gA, mpi, fpi, mN, w, loop_reg, 
                lam_SFR);
    }
	return tmp;
}

double Term::W_T_2pi_nu_3(double q, double gA, double c4, double mpi, double fpi,
        double mN, double w, std::string loop_reg, double lam_SFR)
{
    double gA2 = std::pow(gA,2);
    
    double tmp1 = (c4+1.0/(4.0*mN))*w*w - (gA2/(8.0*mN))*(10.0*mpi*mpi+3.0*q*q);

    double tmp2 = (-gA2*Term::A_gen(q,mpi,loop_reg,lam_SFR))/
        (32.0*M_PI*std::pow(fpi,4));
    
    return tmp2*tmp1;

}

//*****************************************************************************

std::vector<double> Term::W_S_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];
    double c4  = params["c4"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);
    

    std::vector<double> tmp(z_len);
    double q,w;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        tmp[i] = -q*q*Term::W_T_2pi_nu_3(q, c4, gA, mpi, fpi, mN, w, loop_reg,
                lam_SFR);
    }
	return tmp;
}

//*****************************************************************************

std::vector<double> Term::V_LS_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);

    
    std::vector<double> tmp(z_len);
    double q,w_t;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w_t    = Term::w_tilde_f(q,mpi);
        tmp[i] = Term::V_LS_2pi_nu_3(q, gA, mpi, fpi, mN, w_t, loop_reg, 
                lam_SFR);
    }
	return tmp;
}

double Term::V_LS_2pi_nu_3(double q, double gA, double mpi, double fpi,
        double mN, double w_t, std::string loop_reg, double lam_SFR)
{
    double gA4 = std::pow(gA,4);
    
    return (3.0*gA4*w_t*w_t*Term::A_gen(q,mpi,loop_reg,lam_SFR))/
        (32.0*M_PI*mN*std::pow(fpi,4));
}

//*****************************************************************************

std::vector<double> Term::W_LS_2pi_nu_3(double qi, double qo, 
        double* z, int z_len,
        std::unordered_map<std::string,double>& LECs,
        std::unordered_map<std::string,double>& params,
        qs::quantum_channel chn, std::string loop_reg, double lam_SFR)
{
    double gA  = params["gA"];

    double mpi = constants::mpi;
    double fpi = constants::fpi;
    double mN  = ph::get_mN(chn.Tz);
    
    std::vector<double> tmp(z_len);
    double q,w;
    for (int i = 0; i < (int)z_len; i++)
    {
        q      = Term::get_q(qi,qo,z[i]);
        w      = Term::w_f(q,mpi);
        tmp[i] = Term::W_LS_2pi_nu_3(q, gA, mpi, fpi, mN, w, loop_reg, lam_SFR);
    }
	return tmp;
}

double Term::W_LS_2pi_nu_3(double q, double gA, double mpi, double fpi, 
        double mN, double w, std::string loop_reg, double lam_SFR)
{
    double gA2 = std::pow(gA,2);
    
    return (gA2*(1-gA2)*w*w*Term::A_gen(q,mpi,loop_reg,lam_SFR))/
        (32.0*M_PI*mN*std::pow(fpi,4));
}
