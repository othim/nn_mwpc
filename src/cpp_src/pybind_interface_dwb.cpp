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
    ang_int_points_     = 96; // Number of points in angular integration
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
    
    
    // Initialize all the saved T-matrix elements to zero
    std::vector<std::complex<double>> T_arr_0;
    T_arr_0.push_back(std::complex<double>(0.0,0.0));
    T_arr_0.push_back(std::complex<double>(0.0,0.0));
    T_arr_0.push_back(std::complex<double>(0.0,0.0));
    T_arr_0.push_back(std::complex<double>(0.0,0.0));
    
    for (auto chn : chns_)
    {
        saved_T_LO_.push_back  (T_arr_0);
        saved_T_NLO_.push_back (T_arr_0);
        saved_T_N2LO_.push_back(T_arr_0);
        saved_T_N3LO_.push_back(T_arr_0);
    }
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
            ph::get_on_shell_from_matrix(T_DWBA,number_of_p_points_);
    
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
            ph::get_on_shell_from_matrix(T_sum,number_of_p_points_);

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
    
    // LO
    gsl_matrix_complex* V_LO = 
        potentials_[V_LO_name]->get_saved_matrix(q_on_shell, chn, rel_corr_);
    gsl_matrix_complex* TI = 
        LS_Solver_->solve_in_chn_T_fullT_weights(T_lab,chn,V_LO,G0);

    
    // Solve for TI
    gsl_matrix_complex* T_res = nullptr; 
    if (order == 0)
    {
        T_res = TI;
    }
    else if (order == 1)
    {
        T_res = dwba::pw_T_DWBA_PC_NLO(TI, G0, V_NLO);
        gsl_matrix_complex_free(V_NLO);
    }
    else if (order == 2)
    {   
        T_res = dwba::pw_T_DWBA_PC_N2LO(TI, G0, V_NLO, V_N2LO);
        gsl_matrix_complex_free(V_NLO);
        gsl_matrix_complex_free(V_N2LO);
    } 
    else if (order == 3)
    {
        T_res = dwba::pw_T_DWBA_PC_N3LO(TI, G0, V_NLO, V_N2LO, V_N3LO);
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
            ph::get_on_shell_from_matrix(T_res,number_of_p_points_);
    
    if (order==0)
    {
        gsl_matrix_complex_free(T_res);
    } else {
        gsl_matrix_complex_free(T_res);
        gsl_matrix_complex_free(TI);
    }

    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(V_LO);

    return T_arr;
}

std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_BA_T_PC(
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
    
    // LO
    gsl_matrix_complex* V_LO = 
        potentials_[V_LO_name]->get_saved_matrix(q_on_shell, chn, rel_corr_);
    gsl_matrix_complex* TI = 
        LS_Solver_->solve_in_chn_T_fullT_weights(T_lab,chn,V_LO,G0);

    // ****
    // OBS
    // ****
    gsl_matrix_complex_set_zero(TI);
    
    // Solve for TI
    gsl_matrix_complex* T_res = nullptr; 
    if (order == 0)
    {
        T_res = TI;
    }
    else if (order == 1)
    {
        T_res = dwba::pw_T_DWBA_PC_NLO(TI, G0, V_NLO);
        gsl_matrix_complex_free(V_NLO);
    }
    else if (order == 2)
    {   
        T_res = dwba::pw_T_DWBA_PC_N2LO(TI, G0, V_NLO, V_N2LO);
        gsl_matrix_complex_free(V_NLO);
        gsl_matrix_complex_free(V_N2LO);
    } 
    else if (order == 3)
    {
        T_res = dwba::pw_T_DWBA_PC_N3LO(TI, G0, V_NLO, V_N2LO, V_N3LO);
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
            ph::get_on_shell_from_matrix(T_res,number_of_p_points_);
    
    if (order==0)
    {
        gsl_matrix_complex_free(T_res);
    } else {
        gsl_matrix_complex_free(T_res);
        gsl_matrix_complex_free(TI);
    }

    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(V_LO);

    return T_arr;
}

std::vector<std::complex<double>> nn_mwpc_dwb_interface::solve_DWBA_T_PC_full(
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
    
    //gsl_matrix_complex* V_LO = 
    //    potentials_[V_LO_name]->get_matrix(q_on_shell, chn, rel_corr_);
    
    
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
    if (order == 0)
    {
        // Do nothing, since the TI result is correct
    }
    else if (order == 1)
    {
        gsl_matrix_complex* T1 = nullptr;
        
        T1 = dwba::pw_T_DWBA_PC_NLO(TI, G0, V_NLO);
        
        std::vector<std::complex<double>> T_vec_I = 
                ph::get_on_shell_from_matrix(TI,number_of_p_points_);
        std::vector<std::complex<double>> T_vec_1 = 
                ph::get_on_shell_from_matrix(T1,number_of_p_points_);
        
        std::cout << "T from C++" << std::endl;
        std::cout << std::real(T_vec_I[3]) << std::endl;
        std::cout << std::real(T_vec_1[3]) << std::endl;
        
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
            ph::get_on_shell_from_matrix(TI,number_of_p_points_);
    
    gsl_matrix_complex_free(TI);
    gsl_matrix_complex_free(G0);
    gsl_matrix_complex_free(V_LO);


    return T_arr;
}


void nn_mwpc_dwb_interface::solve_save_T_chn_PC(double T_lab, 
        std::vector<int> chn_index_LO, std::vector<int> orders,
        const std::string& V_LO_name, const std::string& V_NLO_name,
        const std::string& V_N2LO_name, const std::string& V_N3LO_name)
{
    // Solve in specified channels and orders and save the result
    energy_saved_ = T_lab;

    #pragma omp parallel
    {
        for (auto order : orders)
        {
            //int tid = omp_get_thread_num();
            //std::cout << "Hello from thread: " << tid << std::endl;
            //std::cout << "order=" << order << std::endl;
            #pragma omp for
            for (int chn_index = 0; chn_index < chns_.size(); chn_index++)
            {
                //std::cout << "chn_index=" << chn_index << std::endl;
                // If a LO channel
                std::vector<std::complex<double>> t;
                if ( std::find(chn_index_LO.begin(), chn_index_LO.end(), chn_index) 
                        != chn_index_LO.end())
                {
                    t = solve_DWBA_T_PC(T_lab,chn_index,order,V_LO_name,V_NLO_name,
                            V_N2LO_name,V_N3LO_name);
                
                } else
                {
                    t = solve_BA_T_PC(T_lab,chn_index,order,V_LO_name,V_LO_name,
                            V_NLO_name,V_N2LO_name);
                
                }
                // Save the the saved on-shell T-arrays
                save_order(order,chn_index,t);
            }
        }
    }
}


void nn_mwpc_dwb_interface::save_DWBA_T_chn_PC(double Tlab,int order, int chn_index,
        std::complex<double> T11, std::complex<double> T12,
        std::complex<double> T22, std::complex<double> T_uncoup)
{
    energy_saved_ = Tlab;
    // Construct T-element
    std::vector<std::complex<double>> t;
    t.push_back(T11);
    t.push_back(T12);
    t.push_back(T22);
    t.push_back(T_uncoup);
    
    save_order(order,chn_index,t);
}

void nn_mwpc_dwb_interface::save_order(
        int order, int chn_index, std::vector<std::complex<double>> t)
{
    // Check if sizes are initialized correctly
    if (chns_.size() != saved_T_LO_.size() 
            || chns_.size() != saved_T_NLO_.size()
            || chns_.size() != saved_T_N2LO_.size() 
            || chns_.size() != saved_T_N3LO_.size())
    {
        std::cerr << "Error, something is wrong with one of the save vectors."
            << std::endl;
    }

    // Set the correct saved T-value
    if (order == 0)
    {
        saved_T_LO_[chn_index] = t;
    } else if (order == 1 )
    {
        saved_T_NLO_[chn_index] = t;
    } else if (order == 2)
    {
        saved_T_N2LO_[chn_index] = t;
    } else if (order == 3)
    {
        saved_T_N3LO_[chn_index] = t;
    } else
    {
        std::cerr << "Error, order is out of bounds." << std::endl;   
    }
}

std::complex<double> nn_mwpc_dwb_interface::observable_from_saved_T(
        const std::string& obs_name, double theta, int order)
{
    if (order > 3 || order < 0)
    {
        std::cerr << "Error: order out of bounds, returning 0." << std::endl;
        return (std::complex<double>)0.0;
    }

    double mu, q_on_shell;
    //std::cout << "Energy saved:" << energy_saved_ << std::endl;

    LS_Solver::get_mu_q_on_shell(energy_saved_,chns_[0], &mu,&q_on_shell);
    
    // Construct the vector of on-shell-T matrix
    std::vector<std::complex<double>*> T_vec;

    // Compute the factor
    const std::complex<double> imag_u(0.0,1.0);
    
    std::complex<double> fac = -2.0*mu*imag_u*M_PI*q_on_shell;

    int chn_index = 0;
    for (auto chn : chns_)
    {
        // Init values to zero
        std::complex<double>* T = new std::complex<double>[3];
        T[0] = (std::complex<double>)0.0;
        T[1] = (std::complex<double>)0.0;
        T[2] = (std::complex<double>)0.0;

        if (order>=0)
        {
            if (chn.coupled)
            {
                T[0] += fac*saved_T_LO_[chn_index][0];
                T[1] += fac*saved_T_LO_[chn_index][1];
                T[2] += fac*saved_T_LO_[chn_index][2];
            } else
            {   
                // To match the different saving format conventions
                T[0] += fac*saved_T_LO_[chn_index][3];
            }
        }
        if (order >= 1)
        {
            if (chn.coupled)
            {
                T[0] += fac*saved_T_NLO_[chn_index][0];
                T[1] += fac*saved_T_NLO_[chn_index][1];
                T[2] += fac*saved_T_NLO_[chn_index][2];
            } else
            {   
                // To match the different saving format conventions
                T[0] += fac*saved_T_NLO_[chn_index][3];
            }
        }
        if (order >= 2)
        {
            if (chn.coupled)
            {
                T[0] += fac*saved_T_N2LO_[chn_index][0];
                T[1] += fac*saved_T_N2LO_[chn_index][1];
                T[2] += fac*saved_T_N2LO_[chn_index][2];
            } else
            {   
                // To match the different saving format conventions
                T[0] += fac*saved_T_N2LO_[chn_index][3];
            }
        }
        if (order >= 3)
        {
            if (chn.coupled)
            {
                T[0] += fac*saved_T_N3LO_[chn_index][0];
                T[1] += fac*saved_T_N3LO_[chn_index][1];
                T[2] += fac*saved_T_N3LO_[chn_index][2];
            } else
            {   
                // To match the different saving format conventions
                T[0] += fac*saved_T_N3LO_[chn_index][3];
            }
        }
        
        //std::cout << "Setting T: " << T[0] << "," << T[1] << "," << T[2] << 
        //    std::endl;

        T_vec.push_back(&T[0]);
        chn_index++;
    }

    // Compute saclay amplitudes 
    std::vector<std::complex<double> > saclay_amplitudes;

    // S = 1-2*i*rho_T*T
    double rho_T = M_PI*q_on_shell*mu; // In the convention used
    
    saclay_amplitudes = sc::compute_Saclay_amplitudes(chns_, T_vec, 
            theta*M_PI/180.0, q_on_shell, rho_T, J_max_in_pot_);
    
    // Print saclay amplitudes
    /*
    std::cout << "Amplitudes:" << std::endl;
    for (auto s : saclay_amplitudes)
    {
        std::cout << std::real(s) << "," << std::imag(s) << std::endl;
    }
    */

    // Compute the observable from the amplitudes
    double obs_value;
    // The A 00kk is an observable that is defined in terms of other vectors
    // than the n,l,m. That is why the _lab function is used.
    if (obs_name == "A 00kk") {
        obs_value = sc::compute_observable_lab(saclay_amplitudes, q_on_shell, 
                obs_name, theta*M_PI/180.0);
    } else {
        obs_value = sc::compute_observable(saclay_amplitudes, q_on_shell, obs_name);
    }

    // Delete the vector of on-shell T
    while (!T_vec.empty())
    {
        delete T_vec.back();
        T_vec.pop_back();
    }

    return obs_value;
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
    Potential_mwpc<gsl_matrix_complex>* pot = load_pre_def_pot(pre_def_name, 
            lam_SFR);
    std::cout << "Loaded pre def potential" << std::endl;
    // Insert the potential in the list of potentials
    potentials_.insert( std::make_pair(potential_name,pot) );

    // Add the potential name to the list of potential names
    potential_names_.push_back(potential_name);
}

void nn_mwpc_dwb_interface::print_potential_info(const std::string& potential_name)
{   
    potentials_[potential_name]->print_LECs_and_params_info();

    // Added to check that the potential is pw decomposed correctly
    
    /*
    double qi = 200;
    double qo = 200;

    // 1S0
    qs::quantum_channel chn = chns_[0];
    double V_arr[6];
    potentials_[potential_name]->calc_element_V_arr(qi,qo,chn,&V_arr[0]);
    
    
    std::cout << "V_arr in " << quantum_channel_to_string(chn) << std::endl;
    std::cout << "qi = " << qi << std::endl;
    std::cout << "qo = " << qo << std::endl;

    std::cout << "[V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]" << std::endl;
    std::cout << std::setprecision(16);
    for (int i=0;i<6;i++)
    {
        double mu = ph::get_mN(0)/2.0;
        double rel_cut = potentials_[potential_name]->get_rel_cut(qi,qo,mu,true);
        std::cout << V_arr[i]*rel_cut << "   ";
    }
    std::cout << "\n\n";
    
    // 3S-D1
    chn = chns_[3];
    potentials_[potential_name]->calc_element_V_arr(qi,qo,chn,&V_arr[0]);
    
    std::cout << "V_arr in " << quantum_channel_to_string(chn) << std::endl;
    std::cout << "qi = " << qi << std::endl;
    std::cout << "qo = " << qo << std::endl;

    std::cout << "[V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]" << std::endl;
    std::cout << std::setprecision(16);
    for (int i=0;i<6;i++)
    {
        double mu = ph::get_mN(0)/2.0;
        double rel_cut = potentials_[potential_name]->get_rel_cut(qi,qo,mu,true);
        std::cout << V_arr[i]*rel_cut << "   ";
    }
    std::cout << "\n\n";
    */
}

void nn_mwpc_dwb_interface::set_LECs_in_potential(const std::string& potential_name, 
            const std::vector<double>& LECs)
{
    // Set the LECs 
    int i=0;

    for (auto& it: potentials_[potential_name]->LECs_in_use_)
    {
         //std::cout << "LEC: " <<  it << "=" << LECs[i] << std::endl;
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

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < chns_.size(); i++)
        {
            qs::quantum_channel chn = chns_[i];
            potentials_[potential_name]->populate_saved_mtx(chn,rel_corr_);
            std::cout << "Saved channel: " << i << std::endl;
        }
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

int nn_mwpc_dwb_interface::get_chn_coupled(int chn_number)
{
    if (!(chn_number < (int)chns_.size()))
    {
        std::string s = "Invalid chhannel number";
        return -1;
    }
    return (int)chns_[chn_number].coupled;
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

Potential_mwpc<gsl_matrix_complex>*  nn_mwpc_dwb_interface::
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
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
        std::string loop_reg    = "DR";

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
    
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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
    
        return pot_complex_weights;
    } else if (pre_def_name == "MWPC_N3LO_SP")
    {
        std::cout << "Adding MWPC_N3LO_SP" << std::endl;
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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
    
        return pot_complex_weights;
    } else if (pre_def_name == "WPC_LO")
    {
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        
        bool inc_weights_in_pot = true; // This is always true
        std::string loop_reg    = "DR";
        
        std::cout << "Creating WPC_LO potential" << std::endl;
        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg,lam_SFR);
;

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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR);
    
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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR);
    
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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
                w_grid_, number_of_p_points_,J_max_in_pot_,
                cutoff_, cut_pow_, sharp_cutoff_, inc_weights_in_pot, 
                cut_on_shell_,loop_reg, lam_SFR);
    
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

        // Make the potential complex
        Potential_mwpc<gsl_matrix_complex>* pot_complex_weights = 
                new Potential_mwpc<gsl_matrix_complex>(terms,ang_int_points_,p_grid_,
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
