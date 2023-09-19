
    

Pot_mwpc<gsl_matrix_complex>* load_pre_def_pot(
        std::string pre_def_name)
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
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);

        return pot_complex_weights;

    }
}
