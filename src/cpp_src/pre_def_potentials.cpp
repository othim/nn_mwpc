#include "pre_def_potentials.h"
    

Potential_mwpc<gsl_matrix_complex>* pre_def_pot::create_pre_def_pot(
        std::string pre_def_name, int ang_int_points_, double* p_grid_, 
        double* w_grid_, int number_of_p_points_, int J_max_in_pot_,
        double cutoff_, int cut_pow_, bool sharp_cutoff_, 
        double sharp_cutoff_add_, bool inc_weights_in_pot_, bool cut_on_shell_,
        ph::constants_struct* program_const_)
{
    if (pre_def_name == "Yamaguchi_1S0")
    {
        /*
         * Yamaguchi_1S0
         *
         * This is a potential in the 1S0 channel
         */
        std::vector<std::string> terms;
        terms.push_back("Yamaguchi_1S0");
        
        bool inc_weights_in_pot = true; // This is always true
        
        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, 
                inc_weights_in_pot, cut_on_shell_,"DR",700.0,program_const_);

        return pot_complex_weights;
    
    } else if (pre_def_name == "Yamaguchi_3S-D1")
    {
        /*
         * Yamaguchi_3S-D1
         *
         * This is a potential in the 3S-D1 channel
         */
        std::vector<std::string> terms;
        terms.push_back("Yamaguchi_3S1");
        terms.push_back("Yamaguchi_3S-D1");
        terms.push_back("Yamaguchi_3D-S1");
        terms.push_back("Yamaguchi_3D1");
        
        bool inc_weights_in_pot = true; // This is always true

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,"DR",700.0,program_const_);

        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_LO_SP")
    {
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        
        bool inc_weights_in_pot = true; // This is always true

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,"DR",700.0,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_NLO_SP")
    {
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("C1S0");
        terms.push_back("D1S0");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_N2LO_SP")
    {
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
        
        // TODO: Maybe add 1pi correction here

        // LO pert corr
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        
        // NLO pert corr
        terms.push_back("D1S0");

        // N2LO contacts
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("E1S0");
        terms.push_back("E3P0");
        terms.push_back("E3P2");
        terms.push_back("E_PF");
        terms.push_back("E_FP");


        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_N3LO_SP")
    {
        // Choose terms in potential
        std::vector<std::string> terms;

        // Subleading TPE without reativistic corrections, i.e. all terms
        // proportional to 1/M_N are excluded
        terms.push_back("V_C_2pi_nu_3_no_rel");
        terms.push_back("W_T_2pi_nu_3_no_rel");
        terms.push_back("W_S_2pi_nu_3_no_rel");
        
        // TODO: Maybe add 1pi correction here

        // LO pert corr
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        
        // NLO pert corr
        terms.push_back("D1S0");

        // N2LO pert corr
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("E1S0");
        terms.push_back("E3P0");
        terms.push_back("E3P2");
        terms.push_back("E_PF");
        terms.push_back("E_FP");

        // N3LO contacts
        terms.push_back("F1S0");


        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_N3LO_SP_REL")
    {
        // Choose terms in potential
        std::vector<std::string> terms;

        // Subleading TPE without reativistic corrections, i.e. all terms
        // proportional to 1/M_N are excluded
        terms.push_back("V_C_2pi_nu_3");
        terms.push_back("W_C_2pi_nu_3");
        terms.push_back("V_T_2pi_nu_3");
        terms.push_back("V_S_2pi_nu_3");
        terms.push_back("W_T_2pi_nu_3");
        terms.push_back("V_LS_2pi_nu_3");
        terms.push_back("W_LS_2pi_nu_3");
        
        // TODO: Maybe add 1pi correction here

        // LO pert corr
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        
        // NLO pert corr
        terms.push_back("D1S0");

        // N2LO pert corr
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("E1S0");
        terms.push_back("E3P0");
        terms.push_back("E3P2");
        terms.push_back("E_PF");
        terms.push_back("E_FP");

        // N3LO contacts
        terms.push_back("F1S0");


        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_LO")
    {
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;
        
        std::cout << "Creating WPC_LO potential" << std::endl;
        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);

        std::cout << "Done creating WPC_LO potential" << std::endl;
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_NLO_DR")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_NLO_SFR")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "SFR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_N2LO_DR")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
        
        terms.push_back("V_C_2pi_nu_3");
        terms.push_back("W_C_2pi_nu_3");
        terms.push_back("V_T_2pi_nu_3");
        terms.push_back("V_S_2pi_nu_3");
        terms.push_back("W_T_2pi_nu_3");
        terms.push_back("W_S_2pi_nu_3");
        terms.push_back("V_LS_2pi_nu_3");
        terms.push_back("W_LS_2pi_nu_3");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_N2LO_SFR_K")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
        
        terms.push_back("V_C_2pi_nu_3");
        terms.push_back("W_C_2pi_nu_3");
        terms.push_back("V_T_2pi_nu_3");
        terms.push_back("V_S_2pi_nu_3");
        terms.push_back("W_T_2pi_nu_3");
        terms.push_back("W_S_2pi_nu_3");
        terms.push_back("V_LS_2pi_nu_3");
        terms.push_back("W_LS_2pi_nu_3");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "SFR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_N2LO_SFR_EM")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
        
        terms.push_back("V_C_2pi_nu_3");
        terms.push_back("W_C_2pi_nu_3");
        terms.push_back("V_T_2pi_nu_3");
        terms.push_back("V_S_2pi_nu_3");
        terms.push_back("W_T_2pi_nu_3");
        terms.push_back("W_S_2pi_nu_3");
        terms.push_back("V_LS_2pi_nu_3");
        terms.push_back("W_LS_2pi_nu_3");
        
        // To convert to the EM convention for the subtrantion of the
        // iterated 1PE
        terms.push_back("V_C_2pi_nu_3_to_EM");
        terms.push_back("W_C_2pi_nu_3_to_EM");
        terms.push_back("V_T_2pi_nu_3_to_EM");
        terms.push_back("V_S_2pi_nu_3_to_EM");
        terms.push_back("W_T_2pi_nu_3_to_EM");
        terms.push_back("W_S_2pi_nu_3_to_EM");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "SFR";
        double lam_SFR = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_N2LO_SFR_NO_REL")
    {
        std::vector<std::string> terms;
        // Pion terms
        terms.push_back("W_T_1pi_nu_0");
        
        terms.push_back("V_T_2pi_nu_2");
        terms.push_back("V_S_2pi_nu_2");
        terms.push_back("W_C_2pi_nu_2");
        
        terms.push_back("V_C_2pi_nu_3_no_rel");
        terms.push_back("W_T_2pi_nu_3_no_rel");
        terms.push_back("W_S_2pi_nu_3_no_rel");
         
        // Contact terms
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");
        terms.push_back("D1P1");
        terms.push_back("D3P1");
        terms.push_back("D1S0");
        terms.push_back("D3S1");
        terms.push_back("D_SD");
        terms.push_back("D_DS");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "SFR";
        double lam_SFR          = 700.0;

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, sharp_cutoff_add_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR,program_const_);
    
        return pot_complex_weights;
    } else 
    {
        std::cout << "Error, undefined potential name: " << pre_def_name 
            << std::endl;
        return nullptr;
    }
}
