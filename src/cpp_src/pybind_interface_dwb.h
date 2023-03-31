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
#include "potential_mwpc.h"
#include "pot_ext.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "scattering.h"
#include "physics_helpers.h"
#include "born_approx.h"

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
/* 
 * Here are private class variables that need to be given as arguments to the 
 * constructor.
 *
 */
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

private:
/*
 * Other private variables that are set and used by memberfunctions
 *
 */        

    LS_Solver* LS_Solver_;

    // This map stores the constructed potentials
    std::unordered_map<std::string, Pot_mwpc<gsl_matrix_complex>> potentials_;

    // This vector contains all the namses of the defined potentials
    std::vector<std::string> potential_names_;


    // Saved data
    std::vector<qs::quantum_channel> chns_;
    double energy_saved_;


private:
    /*
     * Private helper methods
     */
    void get_G0_and_potentials(double T_lab, qs::quantum_channel chn,
            gsl_matrix_complex** G0, gsl_matrix_complex** VI,
            gsl_matrix_complex** VII, std::string VI_name,
            std::string VII_name);

    /*
     *
     * Method that takes out the on-shell values
     *
     * Returns:
     * std::vector with four elements: [T[0,0],T[1,0],T[1,1],Tuncoup], where the 
     * indices are over the 2x2 on-shell T-matrix. Tuncoup is the uncoupled 
     * on-shell T-matrix element.
     */
    std::vector<std::complex<double>> 
            get_on_shell_from_matrix(gsl_matrix_complex* M);


    /*
     * Function that defines potentials and populate potentials_ and
     * potential_names_
     *
     */
    void load_predefined_potentials();

public:

    nn_mwpc_dwb_interface(const std::string& model_name, int J_max_chn, 
            double cutoff, int cut_pow, bool sharp_cutoff, bool pre_comp_pot, 
            bool rel_corr, int number_of_p_points, bool finite_grid,
            bool inc_weights_in_pot_ = false, bool cut_on_shell = true);
    ~nn_mwpc_dwb_interface();
    
    /*
     *
     * This function solves the DWB series and returns the on-shell T-matrix
     *
     */
    std::vector<std::complex<double>>   
            solve_DWBA_T(double T_lab, qs::quantum_channel chn, int order,
            std::string VI_name, std::string VII_name);
    /*
     * This function solves for the full T-matrix for the 
     * the potential V_I + V_II
     *
     * Args:
     * -----
     *
     * T_lab : laboratory kinetic energy of projectile in MeV
     * chn   : object that stores the quantum channel that is to be solved for
     *
     * Returns:
     * --------
     *
     * coupled and uncoupled channels:
     *
     * std::vector with four elements: [T[0,0],T[1,0],T[1,1],Tuncoup], where the 
     * indices are over the 2x2 on-shell T-matrix. Tuncoup is the uncoupled 
     * on-shell T-matrix element.
     *
     * Units of the T-matrix is MeV^{-2} and the normalization conventions 
     * are that there is as given for the LS-Solver in the README.
     */
    std::vector<std::complex<double>> solve_exact_pot_sum_T(
            double T_lab, qs::quantum_channel chn, int order, 
            std::string VI_name, std::string VII_name);

    /*
     * This function is the same as 'solve_exact_pot_sum_T(...)' with the 
     * difference that the whole T-matrix is returned
     *
     * Args:
     * -----
     * T_lab : laboratory kinetic energy of projectile in MeV
     * chn   : object that stores the quantum channel that is to be solved for
     *
     * Returns:
     *
     * T-matrix is row-major format with the same units and conventions as
     * 'solve_exact_pot_sum_T()'.
     */
    gsl_matrix_complex* solve_exact_pot_sum_full_T(
            double T_lab, qs::quantum_channel chn,
            std::string VI_name, std::string VII_name);

    /*
     * TODO
     */
    gsl_matrix_complex* solve_DWBA_full_T(double T_lab, 
            qs::quantum_channel chn, int order,
            std::string VI_name, std::string VII_name);

    /*
     * TODO
     */
    void solve_DWBA_PC_full_T(int order);

    void print_LEC_values(std::string potential_name);
    void print_param_values(std::string potential_name);

    void print_LECs_in_use(std::string potential_name);
    void print_params_in_use(std::string potential_name);

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
