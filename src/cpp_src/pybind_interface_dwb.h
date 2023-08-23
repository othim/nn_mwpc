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
#include <algorithm>

#include "potential_mwpc.h"
#include "potential_ext.h"
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
    std::unordered_map<std::string, Potential_mwpc<gsl_matrix_complex>*> potentials_;

    // This vector contains all the namses of the defined potentials
    std::vector<std::string> potential_names_;


    // Saved data
    std::vector<qs::quantum_channel> chns_;
    double energy_saved_;

    // Saves T-matrix elements. This vector of values conrrespond to the 
    // vector of channels.
    std::vector<std::vector<std::complex<double>> > saved_T_LO_;
    std::vector<std::vector<std::complex<double>> > saved_T_NLO_;
    std::vector<std::vector<std::complex<double>> > saved_T_N2LO_;
    std::vector<std::vector<std::complex<double>> > saved_T_N3LO_;


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
    Potential_mwpc<gsl_matrix_complex>* load_pre_def_pot(
            std::string pre_def_name, double lam_SFR);

public:

    nn_mwpc_dwb_interface(double scale,
            int J_max_chn, double cutoff, int cut_pow, 
            bool sharp_cutoff, bool rel_corr,
            int number_of_p_points, bool finite_grid, bool cut_on_shell,
            bool print);
    ~nn_mwpc_dwb_interface();
    
    /*
     *
     * This function solves the DWB series and returns the on-shell T-matrix
     *
     */
    std::vector<std::complex<double>>   
            solve_DWBA_T_chn(double T_lab, qs::quantum_channel chn, int order,
            const std::string& VI_name, const std::string& VII_name);
    std::vector<std::complex<double>>   
            solve_DWBA_T(double T_lab, int chn_index, int order,
            const std::string& VI_name, const std::string& VII_name);
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
    std::vector<std::complex<double>> solve_exact_pot_sum_T_chn(
            double T_lab, qs::quantum_channel chn,
            const std::string& VI_name, const std::string& VII_name);
    std::vector<std::complex<double>> solve_exact_pot_sum_T(
            double T_lab, int chn_index, 
            const std::string& VI_name, const std::string& VII_name);

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
     * *********************************************************
     * Functions to compute amplitudes for chiral EFT potentials
     * *********************************************************
     */

    /*
     * Function to compute the on-shell T-matrix amplitudes in the chiral PC.
     * Note that if you choose LO (order = 0) this just solved the LS equation 
     * for the potential V_LO and returns the T-matrix.
     *
     * Args:
     * -----
     * T_lab     : laboratory kinetic energy for the projectile in MeV
     * chn_index : index of the quantum channel to compute in
     * order     : the chiral order (0 = LO, 1 = NLO, 2=N2LO, 3=N3LO)
     * V_NxLO    : potential name for each order. Note that you just need to 
     *             define the orders that are used. The others are set by
     *             default arguments
     *
     * Returns:
     * --------
     * std::vector with four elements: [T[0,0],T[1,0],T[1,1],Tuncoup], where the 
     * indices are over the 2x2 on-shell T-matrix. Tuncoup is the uncoupled 
     * on-shell T-matrix element.
     *
     * Units of the T-matrix is MeV^{-2} and the normalization conventions 
     * are that there is as given for the LS-Solver in the README.
     */
    std::vector<std::complex<double>> solve_DWBA_T_PC_full(
            double T_lab, int chn_index, int order,
            const std::string& V_LO_name, const std::string& V_NLO_name,
            const std::string& V_N2LO_name, const std::string& V_N3LO_name);
    
    /*
     * Same as above but just the T-matrix at order=order is returned, not
     * the whole sum.
     */
    std::vector<std::complex<double>> solve_DWBA_T_PC(
            double T_lab, int chn_index, int order,
            const std::string& V_LO_name, const std::string& V_NLO_name,
            const std::string& V_N2LO_name, const std::string& V_N3LO_name);
    
    /*
     * Same as above the LO T-matrix is zero since the LO potential is zero.
     * Note that this means that the V_LO_name potential does not do anything 
     * in this function.
     */
    std::vector<std::complex<double>> solve_BA_T_PC(
            double T_lab, int chn_index, int order,
            const std::string& V_LO_name, const std::string& V_NLO_name,
            const std::string& V_N2LO_name, const std::string& V_N3LO_name);
    
    /*
     * Function that does the same as 'solve_DWBA_T_PC()' but saves the 
     * T-matrix in the saved data.
     */
    void solve_save_T_chn_PC(double T_lab, std::vector<int> chn_index_LO,
            std::vector<int> orders, 
            const std::string& V_LO_name, const std::string& V_NLO_name,
            const std::string& V_N2LO_name, const std::string& V_N3LO_name);
    /*
     * Function that saves a T-matrix that is given. This can be used if one 
     * wants to compute the T-matrix in python and then send it back to C++
     * for further computation.
     */
    void save_DWBA_T_chn_PC(double Tlab,int order, int chn_index,
            std::complex<double> T11, std::complex<double> T12,
            std::complex<double> T22, std::complex<double> T_uncoup);
    /*
     * Computes observables by adding the contributions from all channels 
     * that have a saved T-matrix.
     */
    std::complex<double> observable_from_saved_T(const std::string& obs_name, 
            double theta, int order);
    
    /*
     * Function that saves the on-shell T-matrix element t in the correct saved
     * array.
     */
    void save_order(int order,int chn_index, 
            std::vector<std::complex<double>> t);

    /*
     * ************************************************************
     * Functions to compute other things than scattering properties
     * ************************************************************
     */


    /*
     * This function computes the lowest eigenvalue to the Hamiltonian in the 
     * specified channel. The unit is MeV
     */
    std::vector<double> compute_binding_energy(int chn_number, 
            bool rel_corr, const std::string& V_name);
    
    
    /*
     * *********************************************
     * Functions to create and manipulate potentials
     * *********************************************
     */
    void create_new_potential(const std::string& potential_name, 
            std::string pre_def_name, double lam_SFR);
    
    void print_potential_info(const std::string& potential_name);

    void set_LECs_in_potential(const std::string& potential_name, 
            const std::vector<double>& LECs);
    void set_params_in_potential(const std::string& potential_name, 
            const std::vector<double>& params);

    void save_potential_decomposition(const std::string& potential_name);

    void print_potential_names();
    
    /*
     * **********************************
     * Functions to get various constants
     * **********************************
     */
    double get_on_shell_momentum(double T_lab);
    double get_scale();
    int    get_ang_int_points();
    int    get_momentum_grid_points();
    int    get_chn_len();
    std::string get_chn_LS_term(int chn_number);
    int    get_chn_coupled(int chn_number);

    double get_gA();
    double get_fpi();
    double get_mpi();
    double get_Mp();
    double get_Mn();
};
#endif
