#include "pybind_interface.h"
#include "pybind_interface_dwb.h"

#ifdef PYBIND

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"
    
namespace py = pybind11;

PYBIND11_MODULE(nn_mwpc, m) 
{
    py::class_<nn_mwpc_interface>(m,"nn_mwpc_interface")
        .def(py::init<const std::string&,int,double,int,bool,bool,bool,double,
                bool,bool,bool>())
        .def("solve_LS", &nn_mwpc_interface::solve_LS,
                py::return_value_policy::copy)
        .def("solve_LS_ext_pot", &nn_mwpc_interface::solve_LS_ext_pot,
                py::return_value_policy::copy)
        .def("compute_observable", &nn_mwpc_interface::compute_observable,
                py::return_value_policy::copy)
        .def("get_saved_phase_shifts", &nn_mwpc_interface::get_saved_phase_shifts,
                py::return_value_policy::copy)
        .def("compute_observable_l", &nn_mwpc_interface::compute_observable_l,
                py::return_value_policy::copy)
        .def("compute_phase_shift",&nn_mwpc_interface::compute_phase_shift,
                py::return_value_policy::copy)
        .def("compute_T_on_shell",&nn_mwpc_interface::compute_T_on_shell,
                py::return_value_policy::copy)
        .def("compute_M_element",&nn_mwpc_interface::compute_M_element,
                py::return_value_policy::copy)
        .def("compute_binding_energy",&nn_mwpc_interface::compute_binding_energy,
                py::return_value_policy::copy)
        .def("print_LECs_in_use", &nn_mwpc_interface::print_LECs_in_use,
                py::return_value_policy::copy)
        .def("print_LEC_values", &nn_mwpc_interface::print_LEC_values, 
                py::return_value_policy::copy)
        .def("get_on_shell_momentum",&nn_mwpc_interface::get_on_shell_momentum, 
                py::return_value_policy::copy)
        .def("get_scale", &nn_mwpc_interface::get_scale, 
                py::return_value_policy::copy)
        .def("get_ang_int_points", &nn_mwpc_interface::get_ang_int_points, 
                py::return_value_policy::copy)
        .def("get_momentum_grid_points", &nn_mwpc_interface::get_momentum_grid_points, 
                py::return_value_policy::copy)
        .def("get_chn_len", &nn_mwpc_interface::get_chn_len, 
                py::return_value_policy::copy)
        .def("get_chn_LS_term", &nn_mwpc_interface::get_chn_LS_term, 
                py::return_value_policy::copy)
        .def("get_gA", &nn_mwpc_interface::get_gA, py::return_value_policy::copy)
        .def("get_fpi", &nn_mwpc_interface::get_fpi, py::return_value_policy::copy)
        .def("get_mpi", &nn_mwpc_interface::get_mpi, py::return_value_policy::copy)
        .def("get_Mp", &nn_mwpc_interface::get_Mp, py::return_value_policy::copy)
        .def("get_Mn", &nn_mwpc_interface::get_Mn, py::return_value_policy::copy);
    
    /*
     * ****************
     * DWB module
     * ****************
     */
    py::class_<nn_mwpc_dwb_interface>(m,"nn_mwpc_dwb_interface")
        .def(py::init<double,int,double,int,bool,double,bool,int,
                bool,double,bool,bool,double,double,double,double,double>())
        
        .def("solve_DWBA_T", 
                &nn_mwpc_dwb_interface::solve_DWBA_T,
                py::return_value_policy::copy)
        .def("solve_exact_pot_sum_T", 
                &nn_mwpc_dwb_interface::solve_exact_pot_sum_T,
                py::return_value_policy::copy)
        .def("compute_binding_energy", 
                &nn_mwpc_dwb_interface::compute_binding_energy,
                py::return_value_policy::copy)
        
        .def("solve_DWBA_T_PC_full", 
                &nn_mwpc_dwb_interface::solve_DWBA_T_PC_full,
                py::return_value_policy::copy)
        
        .def("solve_DWBA_T_PC", 
                &nn_mwpc_dwb_interface::solve_DWBA_T_PC,
                py::return_value_policy::copy)
        
        .def("solve_BA_T_PC", 
                &nn_mwpc_dwb_interface::solve_BA_T_PC,
                py::return_value_policy::copy)
        
        .def("solve_save_T_chn_PC", 
                &nn_mwpc_dwb_interface::solve_save_T_chn_PC,
                py::return_value_policy::copy)
        
        .def("save_DWBA_T_chn_PC", 
                &nn_mwpc_dwb_interface::save_DWBA_T_chn_PC,
                py::return_value_policy::copy)
        
        .def("observable_from_saved_T_vec", 
                &nn_mwpc_dwb_interface::observable_from_saved_T_vec,
                py::return_value_policy::copy)
        .def("observable_arr_from_saved_T_vec", 
                &nn_mwpc_dwb_interface::observable_arr_from_saved_T_vec,
                py::return_value_policy::copy)

        .def("create_new_potential", 
                &nn_mwpc_dwb_interface::create_new_potential,
                py::return_value_policy::copy)
        .def("print_potential_info", &nn_mwpc_dwb_interface::print_potential_info,
                py::return_value_policy::copy)
        .def("set_LECs_in_potential", 
                &nn_mwpc_dwb_interface::set_LECs_in_potential,
                py::return_value_policy::copy)
        .def("set_params_in_potential", 
                &nn_mwpc_dwb_interface::set_params_in_potential,
                py::return_value_policy::copy)
        .def("save_potential_decomposition", 
                &nn_mwpc_dwb_interface::save_potential_decomposition,
                py::return_value_policy::copy)
        .def("save_ho_me_decomp", 
                &nn_mwpc_dwb_interface::save_ho_me_decomp,
                py::return_value_policy::copy)
        .def("save_ho_me", 
                &nn_mwpc_dwb_interface::save_ho_me,
                py::return_value_policy::copy)
        .def("save_ho_me_diff_chn", 
                &nn_mwpc_dwb_interface::save_ho_me_diff_chn,
                py::return_value_policy::copy)
        .def("print_potential_names", 
                &nn_mwpc_dwb_interface::print_potential_names,
                py::return_value_policy::copy)
        .def("get_V_matrix", 
                &nn_mwpc_dwb_interface::get_V_matrix,
                py::return_value_policy::copy)
        .def("get_G0_matrix", 
                &nn_mwpc_dwb_interface::get_G0_matrix,
                py::return_value_policy::copy)
        .def("get_DWBA_T_matrix", 
                &nn_mwpc_dwb_interface::get_DWBA_T_matrix,
                py::return_value_policy::copy)
        .def("save_decomp_T_chn_PC", 
                &nn_mwpc_dwb_interface::save_decomp_T_chn_PC,
                py::return_value_policy::copy)
        .def("set_saved_T_vec_from_saved_decomp", 
                &nn_mwpc_dwb_interface::set_saved_T_vec_from_saved_decomp,
                py::return_value_policy::copy)
        .def("set_saved_T_vec_from_saved_orders", 
                &nn_mwpc_dwb_interface::set_saved_T_vec_from_saved_orders,
                py::return_value_policy::copy)

        

        .def("get_on_shell_momentum",&nn_mwpc_dwb_interface::get_on_shell_momentum, 
                py::return_value_policy::copy)
        .def("get_scale", &nn_mwpc_dwb_interface::get_scale, 
                py::return_value_policy::copy)
        .def("get_ang_int_points", &nn_mwpc_dwb_interface::get_ang_int_points, 
                py::return_value_policy::copy)
        .def("get_momentum_grid_points", &nn_mwpc_dwb_interface::get_momentum_grid_points, 
                py::return_value_policy::copy)
        .def("get_chn_len", &nn_mwpc_dwb_interface::get_chn_len, 
                py::return_value_policy::copy)
        .def("get_chn_LS_term", &nn_mwpc_dwb_interface::get_chn_LS_term, 
                py::return_value_policy::copy)
        .def("get_chn_coupled", &nn_mwpc_dwb_interface::get_chn_coupled, 
                py::return_value_policy::copy)
        .def("get_p_points", &nn_mwpc_dwb_interface::get_p_points, 
                py::return_value_policy::copy)
        .def("get_w_points", &nn_mwpc_dwb_interface::get_w_points, 
                py::return_value_policy::copy)
        .def("get_gA", &nn_mwpc_dwb_interface::get_gA, py::return_value_policy::copy)
        .def("get_fpi", &nn_mwpc_dwb_interface::get_fpi, py::return_value_policy::copy)
        .def("get_mpi", &nn_mwpc_dwb_interface::get_mpi, py::return_value_policy::copy)
        .def("get_Mp", &nn_mwpc_dwb_interface::get_Mp, py::return_value_policy::copy)
        .def("get_Mn", &nn_mwpc_dwb_interface::get_Mn, py::return_value_policy::copy);
        
}
#endif

nn_mwpc_interface::nn_mwpc_interface(const std::string& model_name, 
        int J_max_chn, double cutoff, int cut_pow, 
        bool sharp_cutoff, bool pre_comp_pot, bool rel_corr,
        int number_of_p_points,bool finite_grid,bool inc_weights_in_pot,
        bool cut_on_shell, ph::constants_struct* program_const)
{

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    scale_ = 100.0; // Scale of momenutm grid MeV (100)
    number_of_p_points_ = number_of_p_points; // Number of momentum-grid points (60)
    ang_int_points_ = 76; // Number of points in angular integration
    J_max_in_pot_ = 50; // Maximum J that is stored for L-polynomials
    cutoff_ = cutoff; // Cutoff in LS-equation
    cut_pow_ = cut_pow;
    sharp_cutoff_ = sharp_cutoff;
    pre_comp_pot_ = pre_comp_pot; // If pre-computations should be made
    rel_corr_ = rel_corr;
    finite_grid_ = finite_grid;
    finite_grid_max_ = 0.0;
    inc_weights_in_pot_ = inc_weights_in_pot;
    cut_on_shell_ = cut_on_shell;
    
    program_const_ = program_const;
    
    J_pot_ext_cut_ = 5000;
    Pot_ext_aux_ = nullptr;
    
    // For the quantum states
    int J_max = J_max_chn;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = false;
    // ---------------------------------
    // ---------------------------------

    // Initialize physics helpers
    ph::physics_helpers_init();

    // Make GL grid
    if (!finite_grid)
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points_,scale_,&p_grid_,&w_grid_);
    } else 
    {
        // Make GL-grid that is finite. Since we use a sharp cutoff this grid is fine
        double sharp_cut_add = 300.0;
        finite_grid_max_ = cutoff_+sharp_cut_add;
        ph::gauss_legendre_finite_mesh(number_of_p_points_,0,finite_grid_max_,&p_grid_,&w_grid_);
    }

    // These are the pre-determined models
    if ("WPC_LO"==model_name)
    {
        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
                Tz_min, Tz_max, print);
     
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   

        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");

        // Construct potential
        //Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
        //        number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_);
        Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
            number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_,
            300.0,inc_weights_in_pot_,cut_on_shell_,"DR",700.0,program_const_);
        Pot_ext_ = nullptr;

        double ga = 1.29;
        Pot_->params_["gA"] = ga;
        if (pre_comp_pot_)
        {
            // Save potential
            for (auto chn : chns_)
            {
                Pot_->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_,program_const_);

    } else if("MWPC_LO_1"==model_name)
    { 
        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
                Tz_min, Tz_max, print);
         
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   
        
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");

        // Construct potential
        //Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
        //        number_of_p_points_,J_max_in_pot_,cutoff_, cut_pow_,sharp_cutoff_);
        Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
            number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_,
            300.0,inc_weights_in_pot_,cut_on_shell_,"DR",700.0,program_const_);
        Pot_ext_ = nullptr;

        if (pre_comp_pot_)
        {
            // Save potential
            for (auto chn : chns_)
            {
                Pot_->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_,program_const_);

    } else if("MWPC_LO_J"==model_name)
    { 
        /* This potential is the same as the one in 
         * PhysRevC.103.054304.
         * The channels that are included are:
         * 1s0, 3s1-3d1, 1p1, 3p0, 3p1, 3p2-3f2.
         * These are all channels with l \leq 1 and coupled channels are
         * included if some part of it has l \leq 1.
         * The regulator used is exp((p/\Lambda)^4 \theta(\Lambda + 300 - p)
         * where \theta(x) = 1 if x > 0, otherwise 0.
         */
        J_max = 2; // Not changable 
        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
                Tz_min, Tz_max, print);
         
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   
        
        // Delete channel 5 and 6 
        chns_.erase(chns_.begin() + 6);
        chns_.erase(chns_.begin() + 5);
        
        std::cout << "Channels in this (MWPC_LO_J) potential: " << std::endl;
        for (auto& chn : chns_) 
        {
            std::cout << quantum_channel_to_string(chn) << std::endl;
        }

        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("D3P0");
        terms.push_back("D3P2");

        // Construct potential
        //Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
        //        number_of_p_points_,J_max_in_pot_,cutoff_, cut_pow_,sharp_cutoff_);
        Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
            number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_,
            300.0,inc_weights_in_pot_,cut_on_shell_,"DR",700.0,program_const_);
        Pot_ext_ = nullptr;

        if (pre_comp_pot_)
        {
            // Save potential
            for (auto chn : chns_)
            {
                Pot_->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_, program_const_);

    } else if("nijmegen1"==model_name)
    {
        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
                Tz_min, Tz_max, print);
     
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   
        Pot_ = nullptr;
        
        Pot_ext_ = new Potential_ext<gsl_matrix>(p_grid_, number_of_p_points_, cutoff_, 
                &nijm_correct_arg);
        
        std::vector<std::string> terms;
        terms.push_back("W_T_1pi_nu_0");

        // Construct potential
        J_pot_ext_cut_ = 10;
        rel_corr_      = false;
        sharp_cutoff_  = false;

        Pot_ext_aux_ = 
            new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
            number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_,
            300.0,inc_weights_in_pot_,cut_on_shell_,"DR",700.0,program_const_);
        
        if (pre_comp_pot_)
        {
            // Save potential
            for (auto chn : chns_)
            {
                Pot_ext_aux_->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }
        double ga = 1.29;
        Pot_ext_aux_->params_["gA"] = ga;
        //std::cout << Pot_ext_aux_->LECs_in_use_[0] << std::endl;

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_,program_const_);

    } else if ("cdbonn"==model_name)
    {
        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
                Tz_min, Tz_max, print);
     
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   
        Pot_ = nullptr;
        Pot_ext_ = new Potential_ext<gsl_matrix>(p_grid_, number_of_p_points_, cutoff_, &cdbonn_correct_arg);
        
        // Can't precompute this potential (which is kind of stupid...)
        // Construct the LS_Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_,program_const_);
    } else if ("Yamaguchi_1S0" == model_name)
    {

        // Construct the quantum states and quantum channels
        std::vector<qs::quantum_NN_state> states = get_states_NN(0, 0, 
                0, 0, print);
         
        // Construct the quantum scattering channels from the states
        chns_ = get_channels(states, print);   
        
        chns_.erase(chns_.begin() + 1);
        std::cout << "Channels in this (Yamaguchi_1S0) potential: " << std::endl;
        for (auto& chn : chns_) 
        {
            std::cout << quantum_channel_to_string(chn) << std::endl;
        }
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("Yamaguchi_1S0");

        // Construct potential
        Pot_ = new Potential_mwpc<gsl_matrix>(terms,ang_int_points_,p_grid_,w_grid_,
            number_of_p_points_,J_max_in_pot_,cutoff_,cut_pow_,sharp_cutoff_,
            300.0,inc_weights_in_pot_,cut_on_shell_,"DR",700.0,program_const_);
        Pot_ext_ = nullptr;

        if (pre_comp_pot_)
        {
            // Save potential in 1S0
            for (auto chn : chns_)
            {
                Pot_->populate_saved_mtx(chn,rel_corr_); // Realtivistic factor on
            }
        }

        // Construct LS Solver
        LS_Solver_ = new LS_Solver(number_of_p_points_,p_grid_,w_grid_,
                finite_grid_,program_const_);
    } else
    {
        std::cout << "Error, not a valid potential model name" << std::endl;
    }
    //std::cout << J_max_in_pot_ << std::endl;
    std::cout << "Done!" << std::endl;
}
/*
nn_mwpc_interface::initialize()
{
}
*/
nn_mwpc_interface::~nn_mwpc_interface()
{
    // Delete all pointers
    delete LS_Solver_;
    delete Pot_;
    delete Pot_ext_;

    delete[] p_grid_;
    delete[] w_grid_;

    ph::physics_helpers_free();
}
    
void nn_mwpc_interface::solve_LS(double T_lab, std::vector<double> LECs)
{
    // Set the couplings by looping through LECs in use
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }
    // Compute the phase shifts
    phase_shifts_ = compute_phase_shifts(T_lab);
    energy_saved_ = T_lab;
}

void nn_mwpc_interface::solve_LS_ext_pot(double T_lab)
{
    // Compute the phase shifts
    phase_shifts_ = compute_phase_shifts(T_lab);
    energy_saved_ = T_lab;
}

double nn_mwpc_interface::compute_observable(const std::string& name, double angle)
{
    double mu, q_on_shell;
    LS_Solver_->get_mu_q_on_shell(energy_saved_,chns_[0], &mu,&q_on_shell);
    
    // If the observable is omcputed with the optical theorem
    if (name=="SGT" || name=="SGTL" || name=="SGTT") {
        angle = 0;
    }

    // To avoid numerical instability
    if (std::abs(angle-90.0) < 0.001) {
        angle = 90.001;
    }
    // S = 1-2*i*rho_T*T
    double rho_T = M_PI*q_on_shell*mu; // In the convention used
    
    std::vector<std::complex<double> > saclay_amplitudes;
    // Convert the angle to radians. This uses the pre-computed phase
    // shifts phase_shifts_ that are stored in the class.
    saclay_amplitudes = sc::compute_Saclay_amplitudes(chns_, phase_shifts_, 
            angle*M_PI/180.0, q_on_shell, rho_T, J_max_in_pot_,program_const_);
    
    //for (int i=0; i<saclay_amplitudes.size(); i++)
    //{
    //    std::cout << saclay_amplitudes[i] << " ";
    //}
    // Compute the observable from the amplitudes
    double obs;
    // The A 00kk is an observable that is defined in terms of other vectors
    // than the n,l,m. That is why the _lab function is used.
    if (name == "A 00kk") {
        obs = sc::compute_observable_lab(saclay_amplitudes, q_on_shell, name, 
                angle*M_PI/180.0,program_const_);
    } else {
        obs = sc::compute_observable(saclay_amplitudes, q_on_shell, name, 
                program_const_);
    }
    
    return obs;
}

std::vector<double> nn_mwpc_interface::get_saved_phase_shifts(int chn_number)
{
            
    if (!(chn_number < (int)chns_.size()))
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    } 

    std::vector<double> vec;
    vec.push_back(phase_shifts_[chn_number].delta_p);
    vec.push_back(phase_shifts_[chn_number].delta_m);
    vec.push_back(phase_shifts_[chn_number].epsilon);
    vec.push_back(phase_shifts_[chn_number].delta_uncoupled);
    return vec;
}

std::vector<double> nn_mwpc_interface::compute_observable_l(const std::string& name, 
        std::vector<double> angles, std::vector<double> T_lab, std::vector<double> LECs)
{
    // Set the couplings by looping through LECs in use
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }
    // Compute the phase shifts with those couplings
    std::vector<double> obs_vec;
    double mu, q_on_shell;

    for (auto Tl : T_lab)
    {
        std::vector<Phase_shifts_chn> phases_vec = compute_phase_shifts(Tl);
        LS_Solver_->get_mu_q_on_shell(Tl,chns_[0], &mu,&q_on_shell);
        
        for (auto angle : angles)
        {
            // If the observable is omcputed with the optical theorem
            if (name=="SGT" || name=="SGTL" || name=="SGTT") {
                angle = 0;
            }

            // To avoid numerical instability
            if (std::abs(angle-90.0) < 0.001) {
                angle = 90.001;
            }
            double rho_T = M_PI*q_on_shell*mu; // In the convention used
            
            std::vector<std::complex<double> > saclay_amplitudes;
            // Convert the angle to radians. This uses the pre-computed phase
            // shifts phase_shifts_ that are stored in the class.
            saclay_amplitudes = sc::compute_Saclay_amplitudes(chns_, phases_vec, 
                    angle*M_PI/180.0, q_on_shell, rho_T, J_max_in_pot_,
                    program_const_);
                
            // Compute the observable from the amplitudes
            double obs = sc::compute_observable(saclay_amplitudes, q_on_shell, 
                    name,program_const_);

            obs_vec.push_back(obs);
            // If angle independent observable: Ignore loop over angles 
            if (name=="SGT" || name=="SGTL" || name=="SGTT") {
                break;
            }
        }
    }
    /* 
     * This vector contains the observables in a long array where all the
     * angles for one energy is in a row.
     */
    return obs_vec;
}

void nn_mwpc_interface::print_LECs_in_use()
{
    // Print the LECs in the same order as they are set in the compute functions.
    std::cout << "The LECs in the same order as they must be set in the compute functions." << std::endl;
    for (auto s : Pot_->LECs_in_use_)
    {
        std::cout << s << ", ";
    }
    std::cout << std::endl;
}

void nn_mwpc_interface::print_LEC_values()
{
    std::cout << "Printing the LECs and their current values." << std::endl;
    for (auto& it: Pot_->LECs_)
    {
        std::cout << it.first << " = " << it.second << std::endl;
    }
}


std::vector<double> nn_mwpc_interface::compute_phase_shift(int chn_number, double T_lab, 
        std::vector<double> LECs)
{
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }

    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    
    if (chn_number < (int)chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    double mu, q_on_shell;
    LS_Solver_->get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
        
    gsl_matrix* pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);

    /*
    // Some print
    std::cout << "p_grid" << std::endl;
    for (int i=0; i<number_of_p_points_; i++) {
        std::cout << p_grid_[i] << std::endl;
    }
    std::cout << "pot" << std::endl;
    for (int i=0; i<number_of_p_points_; i++) {
        double pot = 0;
        pot = gsl_matrix_get(pot_V_mtx,i,i);
        std::cout << i << "," << pot << std::endl;
    }
    */
    Phase_shifts_chn phases = LS_Solver_->solve_in_chn_R(T_lab,chn,pot_V_mtx);
        
    gsl_matrix_free(pot_V_mtx);
    std::vector<double> phases_vec;
    phases_vec.push_back(phases.delta_p);
    phases_vec.push_back(phases.delta_m);
    phases_vec.push_back(phases.epsilon);
    phases_vec.push_back(phases.delta_uncoupled);

    return phases_vec;
}
/*
std::vector<double> compute_phase_shift_l(int chn_number, 
        std::vector<double> T_lab, std::vector<double> LECs)
{

    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }

    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    
    if (chn_number < chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    double mu, q_on_shell;
    std::vector<double> phases;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < T_lab.size(); i++)
        {
        
        }
    }

    LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
        
    gsl_matrix* pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);

    Phase_shifts_chn phases = LS_Solver_->solve_in_chn_R(T_lab,chn,pot_V_mtx);
        
    gsl_matrix_free(pot_V_mtx);
    std::vector<double> phases_vec;
    phases_vec.push_back(phases.delta_p);
    phases_vec.push_back(phases.delta_m);
    phases_vec.push_back(phases.epsilon);
    phases_vec.push_back(phases.delta_uncoupled);
}
*/

std::vector<double> nn_mwpc_interface::compute_binding_energy( 
        int chn_number, std::vector<double> LECs)
{
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }

    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    
    if (chn_number < (int)chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    gsl_matrix* pot_V_mtx = Pot_->get_matrix_no_onshell(chn, rel_corr_);
    
    ph::eigen_t diag_res = ph::solve_SE(p_grid_, w_grid_, number_of_p_points_, 
            chn, pot_V_mtx,program_const_->Mn,program_const_->Mp);
    
    std::vector<double> eigenvalues;
    for (int i = 0; i < (int)pot_V_mtx->size1; i++) {
        eigenvalues.push_back(GSL_REAL(gsl_vector_complex_get(diag_res.eigenvalues,i)));
    }

    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_complex_free(diag_res.eigenvalues); 
    gsl_matrix_free(pot_V_mtx);
    return eigenvalues;
}

std::vector<double> nn_mwpc_interface::compute_wave_function( 
        int chn_number, std::vector<double> LECs)
{
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }

    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    
    if (chn_number < (int)chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    gsl_matrix* pot_V_mtx = Pot_->get_matrix_no_onshell(chn, rel_corr_);
    
    ph::eigen_t diag_res = ph::solve_SE(p_grid_, w_grid_, number_of_p_points_, 
            chn, pot_V_mtx,program_const_->Mn,program_const_->Mp);
    

    std::vector<double> wave_function;
    
    for (int i = 0; i < (int)pot_V_mtx->size1; i++) {
        // Take the eigenvector corresponding to the lowest eigenvalue.
        //eigenvalues.push_back(GSL_REAL(gsl_vector_complex_get(diag_res.eigenvalues,i)));
    }
    
    // De-allocate the arrays
    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_complex_free(diag_res.eigenvalues); 
    gsl_matrix_free(pot_V_mtx);
    return wave_function;
}
/*
 * Helper functions
 */

std::vector<Phase_shifts_chn> nn_mwpc_interface::compute_phase_shifts(double Tl)
{
    // Allogacate a vector of the same length as numer of channels and 
    // fill with zeros.
    std::vector<Phase_shifts_chn> phases_vec(chns_.size());
    double mu, q_on_shell;
    
    // Make this loop parallel
    //#pragma omp parallel
    //{
        //#pragma omp for
        for (int i = 0; i < (int)chns_.size(); i++)
        {
            //int th_id = omp_get_thread_num();
            //std::cout << "Hello from thread: " << th_id << std::endl;
            qs::quantum_channel chn = chns_[i];

            LS_Solver_->get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
             
            gsl_matrix* pot_V_mtx;
            pot_V_mtx = get_my_potential_matrix(q_on_shell, chn);
            
            Phase_shifts_chn phases = LS_Solver_->solve_in_chn_R(Tl,chn,pot_V_mtx);
        
            gsl_matrix_free(pot_V_mtx);
            phases_vec[i] = phases;
        }
    //}
    return phases_vec;
}

gsl_matrix* nn_mwpc_interface::get_my_potential_matrix(double q_on_shell,
        qs::quantum_channel chn)
{
    gsl_matrix* pot_V_mtx;
    if (Pot_ext_ == nullptr) {
        pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);
    } else {
        if (chn.J < J_pot_ext_cut_) {
            pot_V_mtx = Pot_ext_->get_matrix(q_on_shell, chn,false);
        } else {
            pot_V_mtx = Pot_ext_aux_->get_matrix(q_on_shell, chn, rel_corr_);
        }
    }
    return pot_V_mtx;
}


std::vector<std::complex<double>> nn_mwpc_interface::compute_T_on_shell(
            int chn_number, double T_lab, std::vector<double> LECs)
{
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }

    qs::quantum_channel chn = chns_[0]; // Just to have it initialized
    
    if (chn_number < (int)chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    double mu, q_on_shell;
    LS_Solver_->get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
        
    gsl_matrix* pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);

    std::complex<double>* telem = LS_Solver_->solve_in_chn_T_Telem(T_lab,chn,pot_V_mtx);
        
    gsl_matrix_free(pot_V_mtx);
    std::vector<std::complex<double>> T_elem_vec;
    for (int i=0; i < 4;i++)
    {
        T_elem_vec.push_back(telem[i]);
    }
    delete[] telem;

    return T_elem_vec;

}

std::complex<double> nn_mwpc_interface::compute_M_element(double T_lab, 
        double theta_cm, int S, int mo, int mi)
{
    // Check that the input is OK
    if (S==0) {
        if (mo != 0 || mi != 0) {
            std::cout << "Warning! inconsistent spin indices... returning 0." << std::endl;
            return std::complex<double>(0.0,0.0);
        }
    } else if (S==1) {
        if (mo > 1 || mo < -1 || mi > 1 || mi < -1) {
            std::cout << "Warning! inconsistent spin indices... returning 0." << std::endl;
            return std::complex<double>(0.0,0.0);
        }  
    } else {
        std::cout << "Warning! inconsistent spin indices... returning 0." << std::endl;
        return std::complex<double>(0.0,0.0);
    }
    // Should be OK if program compes to here
    double theta_rad = theta_cm*M_PI/180.0;
    double mu, q_on_shell;
    LS_Solver_->get_mu_q_on_shell(energy_saved_,chns_[0], &mu,&q_on_shell);
    
    double rho_T = M_PI*q_on_shell*mu; // In the convention used
    
    std::complex<double> M_el =  sc::get_M_matrix_p(chns_,phase_shifts_,
            S,mo,mi,std::cos(theta_rad),q_on_shell,rho_T,J_max_in_pot_);

    return M_el;
}

double nn_mwpc_interface::get_on_shell_momentum(double T_lab)
{
    double mu, q_on_shell;
    LS_Solver_->get_mu_q_on_shell(T_lab, chns_[0], &mu, &q_on_shell);
    return q_on_shell;
}

double nn_mwpc_interface::get_scale()
{
    return scale_;
}

int nn_mwpc_interface::get_ang_int_points()
{
    return ang_int_points_;
}

int nn_mwpc_interface::get_momentum_grid_points()
{
    return number_of_p_points_;
}

int nn_mwpc_interface::get_chn_len()
{
    return (int)chns_.size();
} 

std::string nn_mwpc_interface::get_chn_LS_term(int chn_number)
{
    if (!(chn_number < (int)chns_.size()))
    {
        std::string s = "Invalid chhannel number";
        return s;
    }

    return quantum_channel_to_string(chns_[chn_number]);
}

double nn_mwpc_interface::get_gA()
{
    std::cout << "gA not a constant!" << std::endl;
    return -1;
}

double nn_mwpc_interface::get_fpi()
{
    return program_const_->fpi;
}

double nn_mwpc_interface::get_mpi()
{
    return program_const_->mpi;
}

double nn_mwpc_interface::get_Mp()
{
    return program_const_->Mp;
}

double nn_mwpc_interface::get_Mn()
{
    return program_const_->Mn;
}

