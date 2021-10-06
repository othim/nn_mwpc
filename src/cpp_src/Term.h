#ifndef TERM_H
#define TERM_H

#include <string>
#include <vector>
#include <unordered_map>



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

    // Boolean for knowing if the term is isovector or isoscalar
    bool isovector_; // \tau_1 \cdot \tau_2 factor if isovector

    // Function pointer
   
    double (*my_v_alpha_well_def_pw)(double qi, double qo,std::unordered_map<std::string,double> LECs);

    // Boolean to indicate if the term is a lec-term
    bool well_def_pw_;
    LS_term LS_well_def_pw_;

public:
    std::vector<double> (*my_v_alpha)(double qi, double qo, double* z,unsigned int z_len,std::unordered_map<std::string,double> LECs);

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
    std::vector<double> get_v_alpha(double qi, double qo, double* z,unsigned int z_len,std::unordered_map<std::string,double> LECs);
    double get_v_alpha_well_def_pw(double qi, double qo,std::unordered_map<std::string,double> LECs);

    static std::vector<double> v_alpha_OPEP(double qi, double qo, double* z, unsigned int z_len, std::unordered_map<std::string,double> LECs);
    static double mom_C1S0(double qi, double qo, std::unordered_map<std::string,double> LECs);
    static double mom_C3S1(double qi, double qo, std::unordered_map<std::string,double> LECs);

};

#endif