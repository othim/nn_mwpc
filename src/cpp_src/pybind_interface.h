/* pybind_interface.cpp
 * This file is a part of the nn_mwpc project.
 *
 * This file defines a class and some functions that will be acessable from
 * python. The class will serves as a simplified interface to the code.
 *
 * Oliver Thim 2021-12-
 * Department of Physics, Chalmers
 */
#ifndef PYBIND11_INTERFACE
#define PYBIND11_INTERFACE

#include <algorithm>
#include <iostream>
#include <ctime>
#include <string>

#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "scattering.h"
#include "physics_helpers.h"
#include "pot_ext.h"

/* This class will be acessed from python through the bindings in pybind11.
 * The pybind code will be written in sucha a way that C++ will always have
 * ownership over the objects that are created even if it is created in 
 * python.
 *
 * By compiling into the package nn-mwpc and importing the code in
 * python with 'import nn-mwpc' all the functionality will be reached
 * by eg 'myobj = nn-mwpc.nn_mwpc_interface(<args>)' and an object is 
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
class nn_mwpc_interface
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
    bool pre_comp_pot_;
    double* p_grid_;
    double* w_grid_;
    // Objects 
    Potential_mwpc* Pot_;
    Potential_ext* Pot_ext_;

    LS_Solver* LS_Solver_;

    std::vector<qs::quantum_channel> chns_;

    std::vector<Phase_shifts_chn> compute_phase_shifts(double Tl);
public:
    nn_mwpc_interface(const std::string& model_name, int J_max_chn, double cutoff,
            bool pre_comp_pot, bool rel_corr);
    ~nn_mwpc_interface();
    
    std::vector<double> compute_observable(const std::string& name, 
            std::vector<double> angles, std::vector<double> T_lab, std::vector<double> LECs);
    
    /*
     * chn_number is the channel number in the vector chns_
     * T_lab is the lab energy in MeV
     * LECs is the lecs in the correct order as in the potential,
     * the units are specified in the README file.
     *
     * The function return the phase shifts in the Stapp convention.
     */
    std::vector<double> compute_phase_shift(int chn_number, double T_lab, std::vector<double> LECs);
    /*
     * This function computes the lowest eigenvalue to the Hamiltonian in the 
     * specified channel. The unit is MeV
     */
    std::vector<double> nn_mwpc_interface::compute_binding( 
            int chn_number, std::vector<double> LECs);
    
    std::string print_LEC_values();
    std::string print_LECs_in_use();
};

/* Here is the code for the pybind11 interface
 */


//PYBIND11_MODULE(nn_mwpc, m);

/*
 */
#endif
