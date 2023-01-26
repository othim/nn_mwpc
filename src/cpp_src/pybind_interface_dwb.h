/* pybind_interface_dwb.h
 * This file is a part of the nn_mwpc project.
 *
 * This file defines a class and some functions that will be acessable from
 * python. The class will provide functionality to solve the NN scattering 
 * problem in distorted wave perturbation theory.
 *
 * Oliver Thim 2023-01-
 * Department of Physics, Chalmers
 */

#ifndef PYBIND_INTERFACE_DWB
#define PYBIND_INTERFACE_DWB

#include <iostream>
#include <string>
#include <omp.h>

#include "pot_nn_mwpc.h"
#include "pot_ext.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "scattering.h"
#include "physics_helpers.h"

/* This class will be acessed from python through the bindings in pybind11.
 * The pybind code will be written in sucha a way that C++ will always have
 * ownership over the objects that are created even if it is created in 
 * python.
 *
 * By compiling into the package nn-mwpc-dwb and importing the code in
 * python with 'import nn-mwpc-dwb' all the functionality will be reached
 * by eg 'myobj = nn-mwpc-dwb.nn_mwpc_interface(<args>)' and an object is 
 * created from the class. Then one can call 
 * 'phases = myobj.compute_phase_shifts(<args>)' in python and get phase 
 * shifts.
 *
 * You can't access all the functionality from the code without modifying the
 * code in this class. The main idea is to be able to easily acess key 
 * pre-determined features from python, that are explicitly implemented here.
 * 
 * The initialization of the object can take some time since pre-computations 
 * are made. The idea is that the functions that will be called multiple times
 * should be as fast as possible.
 */
class nn_mwpc_dwb_interface
{
private:
/* Here variables that needs to be initialized are stored
 */
    // Constants
    double scale_;
    int ang_int_points_;
    int number_of_p_points_;
    int J_max_in_pot_;   
    bool rel_corr_;
    double cutoff_;
    int cut_pow_;
    bool sharp_cutoff_;
    bool pre_comp_pot_;
    bool inc_weights_in_pot_;
    bool cut_on_shell_;
    double* p_grid_;
    double* w_grid_;

    bool finite_grid_;
    double finite_grid_max_;
    // Objects 
    
    Potential_mwpc* Pot_;
    Potential_ext* Pot_ext_;

    // This potential is taking over the Pot_ext_ in prtial waves with
    // J >= J_pot_ext_cut_.
    Potential_mwpc* Pot_ext_aux_;
    int J_pot_ext_cut_;

    LS_Solver* LS_Solver_;
    
    // Saved data
    std::vector<qs::quantum_channel> chns_;
    std::vector<Phase_shifts_chn> phase_shifts_; 
    double energy_saved_;
    
    // Helper functions
    std::vector<Phase_shifts_chn> compute_phase_shifts(double Tl);
public:

    nn_mwpc_dwb_interface(const std::string& model_name, int J_max_chn, 
            double cutoff, int cut_pow, bool sharp_cutoff, bool pre_comp_pot, 
            bool rel_corr, int number_of_p_points, bool finite_grid,
            bool inc_weights_in_pot_ = false, bool cut_on_shell = true);
    ~nn_mwpc_dwb_interface();
    

    


    void print_LEC_values();
    void print_LECs_in_use();

    double get_on_shell_momentum(double T_lab);
    double get_scale();
    int    get_ang_int_points();
    int    get_momentum_grid_points();
    int    get_chn_len();
    std::string get_chn_LS_term(int chn_number);

    double get_gA();
    double get_fpi();
    double get_mpi();
    double get_Mp();
    double get_Mn();

};
#endif
