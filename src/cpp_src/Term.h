#ifndef TERM_H
#define TERM_H

#include <string>
#include <vector>
#include <unordered_map>


/*
 * Momenta:
 * ----------------------------------------------------------------------------
 * q = p'-p, k = (1/2)(p'+p)
 * p' = qo, p  = qi, is a naming convention in the code.
 
 *
 * The spin structures of the potential are:
 * ----------------------------------------------------------------------------
 * (as eq. 4.7 in M&E Phys. Rep. 503 (2011))
 * (the last spin structure 'sigmak' is added from eq. B.1 in Ebelbaum et al.
 * Nuclear Physics A 747 (2005))
 * (Note that the definitions of the other spin structure agree between the two
 * equations.)
 * 
 * 'C' <-> no spin structure
 * 'S' <-> \sigma_1 \cdot \sigma_2
 * 'LS' <-> (i/2)(\sigma_1 + \sigma_2) (k \times q)
 * 'T' <-> (q \cdot \sigma_1) (q \cdot \sigma_2)
 * 'sigmaL' <-> (\sigma_1 \cdot (q \times k))(\sigma_2 \cdot (q \times k))
 * 'sigmak' <-> (k \cdot \sigma_1) (k \cdot \sigma_2)
 * 
 * Contanct potential and LEC naming conventions
 * ----------------------------------------------------------------------------
 * 'C' -> no momentum dependence
 * 'D' -> quadratic momentum dependence
 * 'E' -> fourth order momentum dependence
 * 'F' -> sixth order momentum dependence 
 * 
 * LEC name | Term in the potential 
 * --------------------------------
 * C1S0     | C1S0
 * C3S1     | C3S1
 * D1S0     | D1S0 (p'^2 + p^2)
 * D3S1     | D3S1 (p'^2 + p^2)
 * D1P1     | D1P1 p'p
 * D3P0     | D3P0 p'p
 * D3P2     | D3P2 p'p
 * D_SD     | D_SD p^2
 * D_DS     | D_SD p'^2 - NOTE, same LEC in both channels
 * E1S0     | E1S0 p'^2p^2
 * E3P0     | E3P0 p'p(p'^2 + p^2)
 * E3P2     | E3P2 p'p(p'^2 + p^2)
 * E_PF     | E_PF p'p(p^2)
 * E_FP     | E_PF p'p (p'^2) - NOTE, same LEC in both channels
 *
 *
 */




struct LS_term {int J; int Li; int Lo; int S;};

class Term
{
private:
    

    // Static member functions for f_int for all possible terms
    // Variable for spin structure
    std::string spin_structure_; // eg. spin-spin, tensor

    // Variable for which term is created 
    std::string term_name_; // eg. OPEP, C1S0
    
    // Which LECs are associated to the terms (independent of channel)
    std::vector<std::string> lecs_in_term_; // Will be filled in by the constructor
    
    // Which params are associated to the terms (independent of channel)
    std::vector<std::string> params_in_term_; // Will be filled in by the constructor

    // Boolean for knowing if the term is isovector or isoscalar
    bool isovector_; // \tau_1 \cdot \tau_2 factor if isovector

    // Function pointer
   
    double (*my_v_alpha_well_def_pw)(double qi, double qo,
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    // Boolean to indicate if the term is a lec-term
    bool well_def_pw_;
    LS_term LS_well_def_pw_;

public:
    std::vector<double> (*my_v_alpha)(double qi, double qo, double* z,unsigned int z_len,
            std::unordered_map<std::string,double> LECs,
            std::unordered_map<std::string,double>& params);

    // Constructor
    Term(std::string name);

    // Destructor
    ~Term();
    bool well_def_pw();
    std::string get_term_name();
    std::string get_spin_structure();
    bool get_isovector();
    LS_term get_LS_term();
    std::vector<std::string> get_lecs_in_term();
    std::vector<std::string> get_params_in_term();
    
    std::vector<double> get_v_alpha(double qi, double qo, double* z,unsigned int z_len,
            std::unordered_map<std::string,double> LECs,
            std::unordered_map<std::string,double>& params);

    double get_v_alpha_well_def_pw(double qi, double qo,
            std::unordered_map<std::string,double> LECs,
            std::unordered_map<std::string,double>& params);

    static std::vector<double> v_alpha_OPEP(double qi, double qo, double* z, unsigned int z_len, 
            std::unordered_map<std::string,double> LECs,
            std::unordered_map<std::string,double>& params);

    static double mom_C1S0(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    static double mom_C3S1(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    static double mom_C3P0(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    static double mom_C3P2(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    static double mom_C3D2(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    static double mom_Yamaguchi_1S0(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    static double mom_Yamaguchi_3S1(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double mom_Yamaguchi_3S_D1(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double mom_Yamaguchi_3D_S1(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double mom_Yamaguchi_3D1(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    /*
     * nu=2 (N2LO) Chiral two-pion exchange contributions in dimensional
     * regularization. As eq. 4.9-4.12 in M&E Phys. Rep. 503 (2011)
     *
     * This is also the same expressions as in eg.
     * E. Epelbaum et al. Nuc. Phys. A 671 (2000) which is also in DR.
     *
     * These potential expressions are using the naming conventions as defined
     * in the top of this file, and in addition 'W' means isovector and 'V' 
     * isoscalar. See the definitions of the spin structures in the top of this
     * file for more information about the naming conventions.
     *
     */
    
    static double w_f(double q, double mpi);   

    static double L_DR(double q, double mpi);

    static double V_T_2pi_nu_2(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double W_C_2pi_nu_2(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    
    /*
     * nu=3 (N3LO) chiral two-pion exchange contributions in dimensional
     * regularization. As eq. 4.13-4.20 in M&E Phys. Rep. 503 (2011).
     */

    static double w_tilde_f(double q, double mpi);

    static double A_DR(double q, double mpi);
    
    static double V_C_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double W_C_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double V_T_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);

    static double V_S_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double W_T_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double W_S_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double V_LS_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
    
    static double W_LS_2pi_nu_3(double qi, double qo, 
            std::unordered_map<std::string,double>& LECs,
            std::unordered_map<std::string,double>& params);
};

#endif
