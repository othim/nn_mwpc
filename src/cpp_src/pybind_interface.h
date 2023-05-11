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
#include <omp.h>

#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "scattering.h"
#include "physics_helpers.h"
#include "potential_ext.h"

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

    /*
     * Constructor where the models are defined. cut_pow it the exponent in 
     * the regularization exp( (p/\Lambda)^cut_pow). 4 or 6 is a good choice
     */
    nn_mwpc_interface(const std::string& model_name, int J_max_chn, 
            double cutoff, int cut_pow, bool sharp_cutoff, bool pre_comp_pot, 
            bool rel_corr, int number_of_p_points, bool finite_grid,
            bool inc_weights_in_pot_ = false, bool cut_on_shell = true);
    ~nn_mwpc_interface();
    
    /*
     * Method that solves the LS equation and saves the phase shifts in
     * this class. The phase shifts are then accessed when an observable 
     * is computed. Thus, it is necessary to call this method before computing
     * observables. This 
     * method solves the LS in every channel which J less than J_max_chn enered
     * in the constructor. The energy is in MeV in the lab-frame  and the LECs 
     * are in MeV^{-n} depending on which LEC it is.
     */
    void solve_LS(double T_lab, std::vector<double> LECs);
    /*
     * Same as above but without LEC argument if the potential in 
     * the model does not have any LECs
     */
    void solve_LS_ext_pot(double T_lab);

    /*
     * Methods that computes observables from the saved phase shifts. The 
     * angles are in the center-of-mass-frame. 
     */
    double compute_observable(const std::string& name, 
            double angle);
    
    std::vector<double> get_saved_phase_shifts(int chn_number); 
    
    /*
     * Methods that take lists of angles and energies as arguments. The
     * angles are in degrees and the energy is in MeV. The LECs are in 
     * MeV^{-n} to some power n dependeing on the LEC. The angles are in the 
     * center-of-mass frame and the energy is in the lab-frame.
     */
    std::vector<double> compute_observable_l(const std::string& name, 
            std::vector<double> angles, std::vector<double> T_lab, 
            std::vector<double> LECs);
   
    /*
     * chn_number is the channel number in the vector chns_
     * T_lab is the lab energy in MeV
     * LECs is the lecs in the correct order as in the potential,
     * the units are specified in the README file.
     *
     * The function return the phase shifts in the Stapp convention.
     */
    std::vector<double> compute_phase_shift(int chn_number, double T_lab, 
            std::vector<double> LECs);
    
    /*
     * Same as 'compute_phase_shift' but accepts a list of energies. This
     * function is parallelized over the loop over energies.
     */ 
    std::vector<double> compute_phase_shift_l(int chn_number, 
            std::vector<double> T_lab, std::vector<double> LECs);
    
    /*
     * This function computes the lowest eigenvalue to the Hamiltonian in the 
     * specified channel. The unit is MeV
     */
    std::vector<double> compute_binding_energy( 
            int chn_number, std::vector<double> LECs);

    /*
     * This function computes the T-matrix elements.
     */
    std::vector<std::complex<double>> compute_T_on_shell(
            int chn_number, double T_lab, std::vector<double> LECs);
    
    /*
     * This function computes M-matrix elements in a coupled spin
     * basis. The M-matrix is computed in the convention that is 
     * described in the article.
     *
     * T_lab    : lab energy of the incomping nucleon in MeV
     * theta_cm : center-of-mass scattering angle in deg
     * S        : total spin of the NN-system, 0,1
     * mo       : spin projection of the outgoing NN-syatem, -1,0,1
     * mi       : spin projection of the incoming NN-system, -1,0,1
     */

    std::complex<double> compute_M_element(double T_lab, 
            double theta_cm, int S, int mo, int mi);

    /*
     * This function returns the wave function of the ground state. In case of
     * a coupled channel the wave funciton for the lowest l is the first
     * block.
     */
    std::vector<double> compute_wave_function(int chn_number, 
            std::vector<double> LECs);


    /*
     * You can choose if you have an external potential or an nn_mwpc potential.
     * There is an additional choice to have an external potential up to a
     * certain maximum J and then have a nn_mwpc potential. The logic for 
     * handling this is taken care of by this function.
     *
     * You should not use the 'Pot_->...' directly in to code, but instead use
     * this function to get the correct potential in each partial wave.
     */
    gsl_matrix* get_my_potential_matrix(double q_on_shell,qs::quantum_channel chn);
    
        
    void print_LEC_values();
    void print_LECs_in_use();

    
    /*
     * Get functions
     */

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
