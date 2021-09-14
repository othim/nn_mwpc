#include <string>

class Term
{
private:
    

    // Static member functions for f_int for all possible terms
    // Variable for spin structure
    std::string spin_structure_; // eg. spin-spin, tensor

    // Variable for which term is created 
    std::string term_name_; // eg. OPEP, C1S0
    
    // Which LECs are associated to the terms (independent of channel)
    std::vector<string> lecs_in_term_; // Will be filled in by the constructor

    // Boolean for knowing if the term is isovector or isoscalar
    bool isovector_; // \tau_1 \cdot \tau_2 factor if isovector

    // Function pointer
    double (*my_v_alpha)(double qi, double qo, double z,std::unordered_map<string,double> LECs);

    // Boolean to indicate if the term is a lec-term
    bool is_lec_;
    struct LS_term {int J; int L; int S} LS_term_lec_;

public:

    // Constructor
    Term(std::string name);

    // Destructor
    ~Term();

};