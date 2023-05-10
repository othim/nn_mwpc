#include "pybind_interface_dwb.h"

#ifdef PYBIND

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"
    
#endif


nn_mwpc_dwb_interface::nn_mwpc_dwb_interface(double scale,
        int J_max_chn, double cutoff, int cut_pow, 
        bool sharp_cutoff, bool rel_corr,
        int number_of_p_points, bool finite_grid, bool cut_on_shell, bool print)
{
    // ************************************************************************
    // ****** CONSTANTS TO CHANGE *********************************************
    // ************************************************************************
    scale_              = scale; // Scale of momenutm grid MeV (100)
    number_of_p_points_ = number_of_p_points; // Number of momentum-grid points (60)
    ang_int_points_     = 76; // Number of points in angular integration
    J_max_in_pot_       = 50; // Maximum J that is stored for L-polynomials
    cutoff_             = cutoff; // Cutoff in LS-equation
    cut_pow_            = cut_pow;
    sharp_cutoff_       = sharp_cutoff;
    rel_corr_           = rel_corr;
    finite_grid_        = finite_grid;
    finite_grid_max_    = 0.0; // Just default value
    cut_on_shell_       = cut_on_shell;
    
    // For the quantum states
    int J_max           = J_max_chn;
    int J_min           = 0;
    int Tz_min          = 0;
    int Tz_max          = 0;
    // ************************************************************************
    // ************************************************************************
    
    // Initialize physics helpers
    ph::physics_helpers_init();
    
    // Make GL grid
    if (!finite_grid)
    {
        ph::gauss_legendre_inf_mesh(
                number_of_p_points_,scale_,&p_grid_,&w_grid_);
    } else 
    {
        // Make GL-grid that is finite. Since we use a sharp cutoff this grid 
        // is fine
        double sharp_cut_add = 300.0;
        finite_grid_max_ = cutoff_+sharp_cut_add;
        ph::gauss_legendre_finite_mesh(number_of_p_points_,0,
                finite_grid_max_,&p_grid_,&w_grid_);
    }
    
    // Construct the quantum states
    if (print) {
        std::cout << "Constructing quantum states..." << std::endl;
    }
    std::vector<qs::quantum_NN_state> states = 
        get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    if (print) {
        std::cout << "Contruction scattering channels..." << std::endl;
    }
    chns_ = get_channels(states, print);   
    

    // Construct a LS_Solver
    LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,finite_grid_);
}


nn_mwpc_dwb_interface::~nn_mwpc_dwb_interface()
{
    delete LS_Solver_;

    for (auto name : potential_names_)
    {
        // TODO
        delete potentials_[name];
    }
    
    // TODO: also delete other stuff

    delete[] p_grid_;
    delete[] w_grid_;
}
 

/*
 * ****************************************************************************
 * Public methods that are a part of the DWB-interface 
 * ****************************************************************************
 */
std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_DWBA_T_chn(
        double T_lab, qs::quantum_channel chn, int order, 
        const std::string& VI_name, const std::string& VII_name)
{
    // Solve for the full T-matrix
    gsl_matrix_complex* T_DWBA = solve_DWBA_full_T(T_lab, chn, order, VI_name,
            VII_name);
    
    // Get the on-shell values
    std::vector<std::complex<double>> T_arr = 
            get_on_shell_from_matrix(T_DWBA);
    
    // Return the on-shell values
    return T_arr;
}

std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_DWBA_T(
        double T_lab, int chn_index, int order, 
        const std::string& VI_name, const std::string& VII_name)
{
    qs::quantum_channel chn = chns_[0];
    if (chn_index>chns_.size()-1)
    {
        std::cout << "Error, too large channel index, returning 0." << std::endl;
        std::vector<std::complex<double>> T_arr;
        T_arr.push_back(std::complex<double>(0,0));
        return T_arr;
    } else
    {
        chn = chns_[chn_index];
        return solve_DWBA_T_chn(T_lab, chn, order, VI_name, VII_name);
    }
}

std::vector<std::complex<double>> 
        nn_mwpc_dwb_interface::solve_exact_pot_sum_T_chn(
        double T_lab, qs::quantum_channel chn, const std::string& VI_name,
        const std::string& VII_name)
{
    // Solve for the full T-matrix
    //std::cout << "solving T-sum\n";
    gsl_matrix_complex* T_sum = solve_exact_pot_sum_full_T(T_lab, chn, VI_name,
            VII_name);
    
    // Get the on-shell values
    std::vector<std::complex<double>> T_arr = 
            get_on_shell_from_matrix(T_sum);

    // Return the on-shell values
    return T_arr;
}

std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_exact_pot_sum_T(
        double T_lab, int chn_index, const std::string& VI_name,
        const std::string& VII_name)
{
    qs::quantum_channel chn = chns_[0];
    if (chn_index>chns_.size()-1)
    {
        std::cout << "Error, too large channel index, returning 0." << std::endl;
        std::vector<std::complex<double>> T_arr;
        T_arr.push_back(std::complex<double>(0,0));
        return T_arr;
    } else
    {
        chn = chns_[chn_index];
        return solve_exact_pot_sum_T_chn(T_lab, chn, VI_name, VII_name);
    }
}


std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_DWBA_T_PC(
        double T_lab, int chn_index, int order,
        const std::string& V_LO_name, const std::string& V_NLO_name="none",
        const std::string& V_N2LO_name="none", const std::string& V_N3LO_name="none")
{
    qs::quantum_channel chn = chns_[0];
    if (chn_index>chns_.size()-1)
    {
        std::cout << "Error, too large channel index, returning 0." << std::endl;
        std::vector<std::complex<double>> T_arr;
        T_arr.push_back(std::complex<double>(0,0));
        return T_arr;
    } else
    {
        chn = chns_[chn_index];
    }

    // Get on-shell momentum and reduced mass
    double q_on_shell,mu;
    LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
    
    // Get G0 vector matrix
    gsl_vector_complex* prop_vec = 
        LS_Solver_->setup_G0_vector_complex(q_on_shell,chn.coupled,mu);
    // Make it to a diagonal matrix
    gsl_matrix_complex* G0 = 
        gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    ph::matrix_from_vector(G0,prop_vec);
    gsl_vector_complex_free(prop_vec);
    
    // Load the potentials that are needed at this order
    
    // LO
    gsl_matrix_complex* V_LO = 
        potentials_[V_LO_name]->get_saved_matrix(q_on_shell, chn, rel_corr_);
    
    // Load the higher order corrections if neccessary
    gsl_matrix_complex* V_NLO  = nullptr; 
    gsl_matrix_complex* V_N2LO = nullptr;
    gsl_matrix_complex* V_N3LO = nullptr;
    if (order>0)
    {
        V_NLO = potentials_[V_NLO_name]->
            get_saved_matrix(q_on_shell, chn, rel_corr_);
    }
    if (order>1)
    {
        V_N2LO = potentials_[V_N2LO_name]->
            get_saved_matrix(q_on_shell, chn, rel_corr_);
    }
    if (order>2)
    {
        V_N3LO = potentials_[V_N3LO_name]->
            get_saved_matrix(q_on_shell, chn, rel_corr_);
    }
    
    // Solve for TI
    gsl_matrix_complex* TI = 
        LS_Solver_->solve_in_chn_T_fullT_weights(T_lab,chn,V_LO,G0);
    
    if (order == 1)
    {
        gsl_matrix_complex* T1 = nullptr;
        
        T1 = dwba::pw_T_DWBA_PC_NLO(TI, G0, V_NLO);
        ph::matrix_add(TI,T1);
        
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(V_NLO);
    }
    else if (order == 2)
    {
        gsl_matrix_complex* T1 = nullptr;
        gsl_matrix_complex* T2 = nullptr;
        
        T1 = dwba::pw_T_DWBA_PC_NLO(TI, G0, V_NLO);
        T2 = dwba::pw_T_DWBA_PC_N2LO(TI, G0, V_NLO, V_N2LO);
        ph::matrix_add(TI,T1);
        ph::matrix_add(TI,T2);
        
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(T2);
        gsl_matrix_complex_free(V_NLO);
        gsl_matrix_complex_free(V_N2LO);
    } 
    else if (order == 3)
    {
        gsl_matrix_complex* T1 = nullptr;
        gsl_matrix_complex* T2 = nullptr;
        gsl_matrix_complex* T3 = nullptr;
        
        T1 = dwba::pw_T_DWBA_PC_NLO (TI, G0, V_NLO);
        T2 = dwba::pw_T_DWBA_PC_N2LO(TI, G0, V_NLO, V_N2LO);
        T3 = dwba::pw_T_DWBA_PC_N3LO(TI, G0, V_NLO, V_N2LO, V_N3LO);
        
        ph::matrix_add(TI,T1);
        ph::matrix_add(TI,T2);
        ph::matrix_add(TI,T3);
        
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(T2);
        gsl_matrix_complex_free(T3);
        gsl_matrix_complex_free(V_NLO);
        gsl_matrix_complex_free(V_N2LO);
        gsl_matrix_complex_free(V_N3LO);
    }
    else
    {
        std::cout << "Error, 'order' must be <=3." << std::endl;
    }
    
    // Get the on-shell values
    std::vector<std::complex<double>> T_arr = 
            get_on_shell_from_matrix(TI);
    
    gsl_matrix_complex_free(TI);
    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(V_LO);


    return T_arr;
}
/*
 * ************************************************************
 * Functions to compute other things than scattering properties
 * ************************************************************
 */


/*
 * This function computes the lowest eigenvalue to the Hamiltonian in the 
 * specified channel. The unit is MeV
 */
std::vector<double> nn_mwpc_dwb_interface::compute_binding_energy(int chn_number, 
        bool rel_corr, const std::string& V_name)
{
    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    if (chn_number < (int)chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }
    gsl_matrix_complex* pot_V_mtx = 
        potentials_[V_name]->get_matrix_no_onshell(chn, rel_corr_);
    
    ph::eigen_t_herm diag_res = ph::solve_SE_complex_weights(p_grid_, w_grid_, 
            number_of_p_points_, chn, pot_V_mtx);
    
    std::vector<double> eigenvalues;
    for (int i = 0; i < (int)pot_V_mtx->size1; i++) {
        eigenvalues.push_back(gsl_vector_get(diag_res.eigenvalues,i));
    }

    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_free(diag_res.eigenvalues); 
    gsl_matrix_complex_free(pot_V_mtx);
    return eigenvalues;
}

/*
 * *********************************************
 * Functions to create and manipulate potentials
 * *********************************************
 */
void nn_mwpc_dwb_interface::create_new_potential(const std::string& potential_name, 
        std::string pre_def_name, double lam_SFR)
{
    // Make a new potential of this type
    Pot_mwpc<gsl_matrix_complex>* pot = load_pre_def_pot(pre_def_name, 
            lam_SFR);
    
    // Insert the potential in the list of potentials
    potentials_.insert( std::make_pair(potential_name,pot) );

    // Add the potential name to the list of potential names
    potential_names_.push_back(potential_name);
}

void nn_mwpc_dwb_interface::print_potential_info(const std::string& potential_name)
{   
    potentials_[potential_name]->print_LECs_and_params_info();
}

void nn_mwpc_dwb_interface::set_LECs_in_potential(const std::string& potential_name, 
            const std::vector<double>& LECs)
{
    // Set the LECs 
    int i=0;

    for (auto& it: potentials_[potential_name]->LECs_in_use_)
    {
         std::cout << "LEC: " <<  it << "=" << LECs[i] << std::endl;
         potentials_[potential_name]->LECs_[it] = LECs[i++];
    }
}

void nn_mwpc_dwb_interface::set_params_in_potential(const std::string& potential_name, 
            const std::vector<double>& params)
{
    // Set the params 
    int i=0;
    for (auto& it: potentials_[potential_name]->params_in_use_)
    {
         potentials_[potential_name]->params_[it] = params[i++];
    }
}

void nn_mwpc_dwb_interface::save_potential_decomposition(
        const std::string& potential_name)
{
    // Save potential in all channels
    for (auto chn : chns_)
    {
        potentials_[potential_name]->populate_saved_mtx(chn,rel_corr_);
    }
}

void nn_mwpc_dwb_interface::print_potential_names()
{
    std::cout << "Constructed potentials:" << std::endl;
    // Save potential in all channels
    for (auto name : potential_names_)
    {
        std::cout << "Potential name: " << name << "\n\n";
        print_potential_info(name);
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/*
 * **********************************
 * Functions to get various constants
 * **********************************
 */

double nn_mwpc_dwb_interface::get_on_shell_momentum(double T_lab)
{    
    double mu, q_on_shell;
    LS_Solver::get_mu_q_on_shell(T_lab, chns_[0], &mu, &q_on_shell);
    return q_on_shell;
}

double nn_mwpc_dwb_interface::get_scale()
{
    return scale_;
}

int nn_mwpc_dwb_interface::get_ang_int_points()
{
    return ang_int_points_;
}

int nn_mwpc_dwb_interface::get_momentum_grid_points()
{
    return number_of_p_points_;
}

int nn_mwpc_dwb_interface::get_chn_len()
{
    return (int)chns_.size();
} 

std::string nn_mwpc_dwb_interface::get_chn_LS_term(int chn_number)
{
    if (!(chn_number < (int)chns_.size()))
    {
        std::string s = "Invalid chhannel number";
        return s;
    }

    return quantum_channel_to_string(chns_[chn_number]);
}

double nn_mwpc_dwb_interface::get_gA()
{
    return constants::gA;
}

double nn_mwpc_dwb_interface::get_fpi()
{
    return constants::fpi;
}

double nn_mwpc_dwb_interface::get_mpi()
{
    return constants::mpi;
}

double nn_mwpc_dwb_interface::get_Mp()
{
    return constants::Mp;
}

double nn_mwpc_dwb_interface::get_Mn()
{
    return constants::Mn;
}

/*
 * ****************************************************************************
 * Public methods that are NOT a part of the DWB-interface 
 * ****************************************************************************
 */

gsl_matrix_complex* nn_mwpc_dwb_interface::solve_exact_pot_sum_full_T(
        double T_lab, qs::quantum_channel chn, std::string VI_name,
        std::string VII_name)
{
    // Get G0 and potentials for the correct on-shell point
    gsl_matrix_complex* G0, *VI, *VII;
    //std::cout << "Get G0 and potentials\n";
    get_G0_and_potentials(T_lab, chn, &G0, &VI, &VII, VI_name, VII_name);

    // Solve LS equation for SUM of potentials, NOTE: VII is now VI+VII!!!
    ph::matrix_add(VII,VI);
    
    // Solve for T
    //std::cout << "Solve LS_equation\n";
    gsl_matrix_complex* T_full = 
        LS_Solver_->solve_in_chn_T_fullT_weights(T_lab,chn,VII,G0);
    
    // Delete and return
    gsl_matrix_complex_free(VI);
    gsl_matrix_complex_free(VII);
    gsl_matrix_complex_free(G0);
    return T_full;
}

gsl_matrix_complex* nn_mwpc_dwb_interface::solve_DWBA_full_T(double T_lab, 
        qs::quantum_channel chn, int order, std::string VI_name,
        std::string VII_name)
{
    // Get G0 and potentials for the correct on-shell point
    gsl_matrix_complex* G0, *VI, *VII;
    get_G0_and_potentials(T_lab, chn, &G0, &VI, &VII, VI_name,VII_name);

    // Solve for TI
    gsl_matrix_complex* TI = 
        LS_Solver_->solve_in_chn_T_fullT_weights(T_lab,chn,VI,G0);
    
    // Solve DWB series
    gsl_matrix_complex* T_DWBA = dwba::pw_T_DWBA(order,TI,VI,VII,G0);
    
    // Delete and return
    gsl_matrix_complex_free(VI);
    gsl_matrix_complex_free(VII);
    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(TI);
    return T_DWBA;
}



/*void nn_mwpc_dwb_interface::solve_DWBA_PC_full_T(int order)
{

}
*/

/*
 * ****************************************************************************
 * Pivate helper methods
 * ****************************************************************************
 */

void nn_mwpc_dwb_interface::get_G0_and_potentials(double T_lab, 
        qs::quantum_channel chn, gsl_matrix_complex** G0, 
        gsl_matrix_complex** VI, gsl_matrix_complex** VII,
        std::string VI_name, std::string VII_name)
{
    // Get on-shell momentum and reduced mass
    double q_on_shell,mu;
    LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
    
    // Get G0 vector matrix
    gsl_vector_complex* prop_vec = 
        LS_Solver_->setup_G0_vector_complex(q_on_shell,chn.coupled,mu);
    // Make it to a diagonal matrix
    gsl_matrix_complex* G0_loc = 
        gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    ph::matrix_from_vector(G0_loc,prop_vec);
    gsl_vector_complex_free(prop_vec);

    // Get the sum of VI and VII with weights
    gsl_matrix_complex* VI_loc = 
        potentials_[VI_name]->get_saved_matrix(q_on_shell, chn, rel_corr_);

    gsl_matrix_complex* VII_loc = 
        potentials_[VII_name]->get_saved_matrix(q_on_shell, chn, rel_corr_);
    
    // Set the value of the external pointers
    *G0  = G0_loc;
    *VI  = VI_loc;
    *VII = VII_loc;
}

std::vector<std::complex<double>> nn_mwpc_dwb_interface::
    get_on_shell_from_matrix(gsl_matrix_complex* M)
{
    gsl_complex T_mm,T_pm,T_pp,T_uncoup;
    // If coupled channel
    if ((int)M->size1==(int)(2*number_of_p_points_+2))
    {
        T_mm = gsl_matrix_complex_get(M,number_of_p_points_,number_of_p_points_);
        T_pm = gsl_matrix_complex_get(M,2*number_of_p_points_+1,number_of_p_points_);
        T_pp = gsl_matrix_complex_get(M,2*number_of_p_points_+1,2*number_of_p_points_+1);
        T_uncoup = gsl_complex_rect(0.0,0.0);

    } else
    {
        T_uncoup = gsl_matrix_complex_get(M,number_of_p_points_,number_of_p_points_);
        T_mm = gsl_complex_rect(0.0,0.0);
        T_pm = gsl_complex_rect(0.0,0.0);
        T_pp = gsl_complex_rect(0.0,0.0);
    
    }
    std::vector<std::complex<double>> T_arr;
    T_arr.push_back(std::complex<double>(GSL_REAL(T_mm),GSL_IMAG(T_mm)));
    T_arr.push_back(std::complex<double>(GSL_REAL(T_pm),GSL_IMAG(T_pm)));
    T_arr.push_back(std::complex<double>(GSL_REAL(T_pp),GSL_IMAG(T_pp)));
    T_arr.push_back(std::complex<double>(GSL_REAL(T_uncoup),GSL_IMAG(T_uncoup)));
    return T_arr;
}

Pot_mwpc<gsl_matrix_complex>*  nn_mwpc_dwb_interface::
        load_pre_def_pot(std::string pre_def_name, double lam_SFR)
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
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);

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
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_NLO_SP")
    {
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("C1S0");
        terms.push_back("D1S0");
        
        bool inc_weights_in_pot = true; // This is always true

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);
    
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

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_LO")
    {
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        
        bool inc_weights_in_pot = true; // This is always true

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_);

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

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
    
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

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
    
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

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_N2LO_SFR")
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

        // Make the potential complex
        Pot_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Pot_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR);
    
        return pot_complex_weights;
    } else 
    {
        std::cout << "Error, undefined potential name: " << pre_def_name 
            << std::endl;
        return nullptr;
    }
}
