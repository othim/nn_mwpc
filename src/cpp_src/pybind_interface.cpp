#include "pybind_interface.h"

#ifdef PYBIND

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(nn_mwpc, m) 
{
    py::class_<nn_mwpc_interface>(m,"nn_mwpc_interface")
        .def(py::init<const std::string&,int,double,bool,bool>())
        .def("compute_observable", &nn_mwpc_interface::compute_observable,
                py::return_value_policy::copy)
        .def("compute_phase_shift",&nn_mwpc_interface::compute_phase_shift,
                py::return_value_policy::copy)
        .def("print_LECs_in_use", &nn_mwpc_interface::print_LECs_in_use,py::return_value_policy::copy)
        .def("print_LEC_values", &nn_mwpc_interface::print_LEC_values, py::return_value_policy::copy);
}
#endif

nn_mwpc_interface::nn_mwpc_interface(const std::string& model_name, 
        int J_max_chn, double cutoff, bool pre_comp_pot, bool rel_corr)
{

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    scale_ = 100.0; // Scale of momenutm grid MeV
    number_of_p_points_ = 100; // Number of momentum-grid points
    ang_int_points_ = 96; // Number of points in angular integration
    J_max_in_pot_ = 50; // Maximum J that is stored for L-polynomials
    cutoff_ = cutoff; // Cutoff in LS-equation
    pre_comp_pot_ = pre_comp_pot; // If pre-computations should be made
    rel_corr_ = rel_corr;

    // For the quantum states
    int J_max = J_max_chn;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = true;
    // ---------------------------------
    // ---------------------------------

    // Initialize physics helpers
    ph::physics_helpers_init();

    // Construct the quantum states and quantum channels
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, 
            Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    chns_ = get_channels(states, print);   

    // Make GL grid
    ph::gauss_legendre_inf_mesh(number_of_p_points_,scale_,&p_grid_,&w_grid_);
    
    // These are the pre-determined models
    if ("WPC_LO"==model_name)
    {
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("OPEP");
        terms.push_back("C1S0");
        terms.push_back("C3S1");

        // Construct potential
        Pot_ = new Potential_mwpc(terms,ang_int_points_,p_grid_,w_grid_,
                number_of_p_points_,J_max_in_pot_,cutoff_);
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
        LS_Solver_ = new LS_Solver(chns_,number_of_p_points_, scale_,
                cutoff_,rel_corr_);

    } else if("MWPC_LO_1"==model_name)
    { 
    
        // Choose terms in potential
        std::vector<std::string> terms;
        terms.push_back("OPEP");
        terms.push_back("C1S0");
        terms.push_back("C3S1");
        terms.push_back("C3P0");
        terms.push_back("C3P2");

        // Construct potential
        Pot_ = new Potential_mwpc(terms,ang_int_points_,p_grid_,w_grid_,
                number_of_p_points_,J_max_in_pot_,cutoff_);
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
        LS_Solver_ = new LS_Solver(chns_,number_of_p_points_, scale_,
                cutoff_,rel_corr_);

    } else if("nijmegen1"==model_name)
    {
        // TODO: implement this 
        // Construct potential
        Pot_ = nullptr;
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
        LS_Solver_ = new LS_Solver(chns_,number_of_p_points_, scale_,
                cutoff_,rel_corr_);

    } else
    {
        std::cout << "Error, not a valid model_name" << std::endl;
    }
    //std::cout << J_max_in_pot_ << std::endl;
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
    
std::vector<double> nn_mwpc_interface::compute_observable(const std::string& name, 
        std::vector<double> angles, std::vector<double> T_lab, std::vector<double> LECs)
{
    // Compute phase shifts in all channels for the given energy
    //std::cout << J_max_in_pot_ << std::endl; 
    //std::cout << scale_ << std::endl;
    //std::cout << ang_int_points_ << std::endl;
    // Set the couplings by looping through LECs in use
    int i = 0;
    for (auto& it: Pot_->LECs_in_use_)
    {
         Pot_->LECs_[it] = LECs[i++];
    }
    // Compute the phase shifts with those couplings
    std::vector<double> obs_vec;
    double mu, q_on_shell, rho_T;

    for (auto Tl : T_lab)
    {
        std::vector<Phase_shifts_chn> phases_vec = compute_phase_shifts(Tl);
        LS_Solver::get_mu_q_on_shell(Tl,chns_[0], &mu,&q_on_shell);
        double obs;
        if (name != "SGT")
        {
            for (auto ang : angles)
            {
                if (std::abs(ang-90.0) < 0.001) {
                    ang = 90.001;
                }
                rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/
                        (constants::Mn+constants::Mp); // In the convention used in the code
                //std::cout << J_max_in_pot_ << std::endl;            
                std::vector<std::complex<double> > saclay_amplitudes;
                saclay_amplitudes = sc::compute_Saclay_amplitudes(chns_, phases_vec, 
                        ang*M_PI/180.0, q_on_shell, rho_T, J_max_in_pot_);
            
                // Compute the observable from the amplitudes
                obs = sc::compute_observable(saclay_amplitudes, name);
                obs_vec.push_back(obs);
            }
        } else 
        {
            rho_T = M_PI*q_on_shell*constants::Mn*constants::Mp/
                    (constants::Mn+constants::Mp); // In the convention used in the code
            
            // Compute the observable from the amplitudes
            obs = sc::compute_total_cross_section(chns_, phases_vec, q_on_shell, rho_T, J_max_in_pot_);
            obs_vec.push_back(obs);
        } 
    }
    /* 
     * This vector contains the observables in a long array where all the
     * angles for one energy is in a row.
     */
    return obs_vec;
}

std::string nn_mwpc_interface::print_LECs_in_use()
{
    // Print the LECs in the same order as they are set in the compute functions.
    std::string output;
    output.append(" Print the LECs in the same order as they are set in the compute functions. \n");
    for (auto s : Pot_->LECs_in_use_)
    {
        output.append(s);
        output.append("\n");
    }
    return output;
}

std::string nn_mwpc_interface::print_LEC_values()
{
    std::string output;
    output.append("Printing the LECs and their current values. \n");
    for (auto& it: Pot_->LECs_)
    {
         output.append(it.first + " = " +std::to_string(it.second) +  "\n");
    }
    return output;
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
    
    if (chn_number < chns_.size())
    {
        chn = chns_[chn_number];
    } else 
    {
        std::cout << "Error: chn_number out of range" << std::endl;
    }

    double mu, q_on_shell;
    LS_Solver::get_mu_q_on_shell(T_lab, chn, &mu, &q_on_shell);
        
    gsl_matrix* pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);

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
 * Helper functions
 */

std::vector<Phase_shifts_chn> nn_mwpc_interface::compute_phase_shifts(double Tl)
{
    std::vector<Phase_shifts_chn> phases_vec;
    double mu, q_on_shell;
    for (auto chn : chns_)
    {
        LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
        
        gsl_matrix* pot_V_mtx = Pot_->get_saved_matrix(q_on_shell, chn, rel_corr_);

        Phase_shifts_chn phases = LS_Solver_->solve_in_chn_R(Tl,chn,pot_V_mtx);
        
        gsl_matrix_free(pot_V_mtx);
        phases_vec.push_back(phases);
    }
    return phases_vec;
}
