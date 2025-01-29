/* potential_mwpc.h
   This file is part of the nn_mwpc project.

   The potential class creats nucleon-nucleon potential objects in a specific
   partial wave. The potential object is flexible in the sense that you can
   choose to create potentials that includes hand-picked contributions.
   By precomputing matrices the function call to compute the potential matrix
   for a given set of LECs should be quite efficient.

   The class ha functionality to be able to compute and return specific matrix
   elements.

   General information about definitions, convention and notation are
   specified in comments above the respective function decalarations.

   Oliver Thim 2021-09 --
   Department of Physics, Chalmers
*/

#ifndef POTENTIAL_MWPC_H
#define POTENTIAL_MWPC_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
//#include <filesystem>
#include "gsl_sf_legendre.h"
#include "gsl_integration.h" 
#include "gsl_matrix.h"
#include "gsl_matrix_complex_double.h"
#include "Term.h"
#include "quantum_states.h"
#include "physics_helpers.h"
#include "potential.h"

//#define ENABLE_DEBUG

template <class gsl_m>
class Potential_mwpc : public potential<gsl_m>
{
private:
    // The matrices are saved in a specific channel and the class
    // should be able to store a list of W's for different channels.
    // By calling populate_saved_matrices() the matrices are stored in memory.
    std::map<qs::quantum_channel, std::map<std::string, gsl_m*>, qs::comp> saved_matrices_;

    // Variables to store the Gauss-Legendre grid for the potential
    double* p_grid_;
    double* w_grid_;
    std::size_t mom_grid_size_;

    // Momentum cutoff for the potential. This is used to regulate 
    // the potential elements for solving the LS-equation in a more
    // stable way.
    double cutoff_Lambda_;

    // The power that Lambda and p is raised to
    int cut_pow_;
    bool sharp_cutoff_;
    double sharp_cutoff_add_;
    // This is the maximum J of the channels the potential class can
    // calculate matrix elements in, since Lagandre polynomials are just
    // precomputed and stored to this order.


    std::string loop_reg_; // == 'DR' or 'SFR'
    double lam_SFR_; // SFR cutoffs in MeV, used if loop_reg_ == 'SFR'


    unsigned int J_max_;

    // If the weights and momentum vectors should be included in the potential
    bool inc_grid_weights_in_pot_;

    // If the cutoff also should affect the on-shell part.
    bool cut_on_shell_;
    // Pointer to the array of the stored polynomials
    double** stored_Legendre_polynomials_;

    unsigned int N_GLI_PWA_; // Number of points in Gauss-Legandre angular integration

    gsl_integration_fixed_workspace* int_ang_; // Need to be saved to not delete pointers
    double* z_mesh; // GL integration points
    double* w_z_mesh; // GL integration weights
    unsigned int len_z_mesh; // GL integration number of points

    ph::constants_struct* program_const_;

private:
   /*
      ---------------------------
      Some private helper-methods
      ---------------------------
   */


   /*
      This method computes the \mathcal{A}-integrals (see Erkelenz). The method takes an
      array of v_alpha-values (see Erkelenz) to speed up the computation.
   */
   double compute_A_integral(double qi, double qo, int J,int l, 
           std::vector<double> v_alpha_arr);
   
   /*
      This method performs a partial-wave-decomposition (partial-wave-analysis) of the 
      potential terms v_alpha funciton. This is done by first computing the A-intagrals
      using compute_A_integral(), and then computing the matrix elemetnts of the form

      V_arr = [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
      where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc

      following the formulas in Erkelenz for a given spin structure. (NOTE some formulas are wrong
      is is pointed out in an appendix in one of Machleidts papers)
   */
   void pwa(double qi,double qo, bool coupled, int J_int, std::string spin_struct,
           bool isovector,std::vector<double>& v_alpha_arr, double* V_arr);
   
   
   double pot_OPEP_mom(double qo,double qi, double z);
    
   double get_total_rel_cut_weight_factor(double p_in, int j, 
           double p_out, int i, double mu, bool rel_correction);
   
   int isoFac(int L, int S);
public:
   /*
      This method clears the saved matrices from memory.
   */
   void clear_saved_matrices();
   
   double get_rel_cut(double p_in,
           double p_out, double mu, bool rel_correction);
   
    /*
      Some class variables that ideally should be private with getters and setters, but in the
      interest of computation time they are public for the time beeing. (The user of this 
      code is trusted!)

      The difference between params_ and LECs_ is that the code is assuming that
      all LEC dependence is linear in the potential, so it can be split.

      The params_ are more static variables of a potentials that is not intended
      to change after the potential is cunstructed. If some of these non-linear
      paramters are changed the storage of the precomputed matrices need to be 
      redone.
    */
    std::vector<std::string> LEC_names_; // List of ALL possible LEC names
    std::unordered_map<std::string, double> LECs_; // List of lecs  and their values
    std::vector<std::string> LECs_in_use_; // List of lec names of lecs in use in this potential
    
    std::vector<std::string> param_names_;
    std::unordered_map<std::string, double> params_; // List of lecs  and their values
    std::vector<std::string> params_in_use_; // List of lec names of lecs in use in this potential
    
    std::vector<Term> terms_in_pot_; // Terms in the potential

    


    /* 
      Constructor
    */
    Potential_mwpc(){} // Dummy default constructor

    Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,
           double* p_grid, double* w_grid, 
           std::size_t grid_size,unsigned int J_max, 
           double cutoff_Lambda, int cut_pow, bool sharp_cutoff,
           double sharp_cutoff_add,
           bool inc_grid_weights_in_pot, bool cut_on_shell,
           std::string loop_reg, double lam_SFR,
           ph::constants_struct* program_const);

    /* 
      Destructor
    */ 
    ~Potential_mwpc();
    /*
      This function computes matrix elements of the potential. NOTE that this is done with the
      CURRENT values of LECs_, so be sure to set them to the correct values prior
      to computation. There is no guarantee that the LECs_ values are laft untouched by
      routines in the class.

      The elements are on the form [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
      where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc
      NOTE: There is no minus sign on the off-diagonal elements as in some conventions!
    */
    void calc_element_V_arr(double qi,double qo, qs::quantum_channel chn, 
            double* V_arr);

    /*
     * Same as above bu includes the cutoff, relativistic and weights factors.
     * NOTE! Can not be used when cut_on_shell_ = false.
     */
    void calc_element_V_arr_full(double qi,double qo, qs::quantum_channel chn,
           bool rel_correction, bool inc_reg_cut_and_rel, double* V_arr);
    
    /*
     *
     * NOTE! Can not be used when cut_on_shell_ = false
     */
    double calc_element_LSJ_full(double p, double pp, int L, int Lp, int S, 
            int J, int T, int Tz);

    /* 
      This function returns a mom_grid_size_ + 1 x mom_grid_size_ + 1 matrix IF the channel is uncoupled
      and 2*mom_grid_size_ + 2 x 2*mom_grid_size_ +2 if the channel is coupled.

      The last momentum point (in the respective blocks in the coupled case) is the on-shell
      q_on_shell passed to the function.
    */
    gsl_m* get_matrix(double q_on_shell,qs::quantum_channel chn, bool rel_correction);

    /*
      This function returns a mom_grid_size_ x mom_grid_size_  matrix IF the channel is uncoupled
      and 2*mom_grid_size_ x 2*mom_grid_size_  if the channel is coupled. You probaly dont want 
      relativistic factors here!!

      This function is contrary to get_matrix() useful when solving the Schrödinger equation
      for the given potential.
    */
    gsl_m* get_matrix_no_onshell(qs::quantum_channel chn, bool rel_correction);

    /*
      This function populates the saved matrices (saved_matrices_) in the given channel.
      It is necessary to run this function BEFORE calling get_saved_matrix() in the given 
      channel. If it is not done, get_saved_matrix() will return an error and terminate.
    */
    void populate_saved_mtx(qs::quantum_channel chn, bool rel_correction);

    /*
      This functions output should be identical to get_matrix(), but should take less time
      since it computes the full matrix as a sum of the saved matrices with the 
      appropriate LECs.
    */
    gsl_m* get_saved_matrix(double q_on_shell, qs::quantum_channel chn, bool rel_correction);
    

    /*
     * Print info about LECs and params in the potential
     */
    void print_LECs_and_params_info();
    
    /*
     * ************************************************************************
     * This part contains functionality to export the potential matrix elements
     * in an harmonic oscillator basis.
     * ************************************************************************
     *
     */
    
    /*
     *
     * Save harmonic oscillator matrix elements.
     *
     */
    void save_ho_me_decomp(std::string save_dir, int Nmax, int hbar_omega,
            bool from_saved_mtx);
    void print_meta_data(std::string file_name, int Nmax, int hbar_omega,
            bool append);
    void save_ho_me_chns(std::string file_name, int Nmax, 
            int hbar_omega, bool from_saved_mtx, bool all_chns, 
            std::vector<qs::quantum_channel> chns = 
            std::vector<qs::quantum_channel>(), 
            bool print_zero_in_unused_chn = false);

    int  save_ho_me_from_saved(std::string file_name, std::vector<double> LECs);

    // <l|V|lp>_sjt. l - outgoing, lp - ingoing.
    void get_chn_block_from_qn(int L, int Lp, int S, int J, int T, 
            qs::quantum_channel* chn, int* block_indexT);

    double compute_HO_matrix_el(int no, int Lo, int ni, 
        int Li, int S, int J, int T, int Tz, double* p_grid, 
        double* w_grid, int num_grid_points, double mN, double Omega);

    void get_sub_matrix();

 };
#endif
