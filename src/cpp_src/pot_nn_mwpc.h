
#include "quantum_states.h"
#include "Constants.h"
#include "gsl_matrix.h"
#include <iostream>
#include <cmath>/*
   pot_nn_mwpc.h
   This file is part of the nn_mwpc project.

   The potential class creats nucleon-nucleon potential objects in a specific
   partial wave. The potential object is flexible in the sense that you can
   choose to create potentials that includes hand-picked contributions.
   By precomputing matrices the function call to compute the potential matrix
   for a given set of LECs should be quite efficient.

   The class ha functionality to be able to compute and return specific matrix
   elements.

   Oliver Thim 2021-09 --
   Department of Physics, Chalmers
*/

#ifndef POT_NN_MWPC_H
#define POT_NN_MWPC_H

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "gsl_matrix.h"
#include "Term.h"
#include "Constants.h"
#include "quantum_states.h"

//#define ENABLE_DEBUG


class Potential_mwpc
{
private:
   // The matrices are saved in a specific channel and the class
   // should be able to store a list of W's for different channels.
   // By calling populate_saved_matrices() the matrices are stored in memory.
   std::map<qs::quantum_channel, std::map<std::string, gsl_matrix*>, qs::comp> saved_matrices_;
   
   // Variables to store the Gauss-Legendre grid for the potential
   double* p_grid_;
   double* w_grid_;
   std::size_t mom_grid_size_;

   // Momentum cutoff for the potential. This is used to regulate 
   // the potential elements for solving the LS-equation in a more
   // stable way.
   double cutoff_Lambda_;

   // This is the maximum J of the channels the potential class can
   // calculate matrix elements in, since Lagandre polynomials are just
   // precomputed and stored to this order.
   unsigned int J_max_;

   // Pointer to the array of the stored polynomials
   double** stored_Legendre_polynomials_;

   unsigned int N_GLI_PWA_; // Number of points in Gauss-Legandre angular integration

   gsl_integration_fixed_workspace* int_ang_; // Need to be saved to not delete pointers
   double* z_mesh; // GL integration points
   double* w_z_mesh; // GL integration weights
   unsigned int len_z_mesh; // GL integration number of points
private:
   /*
      ---------------------------
      Some private helper-methods
      ---------------------------
   */

   /*
      This method clears the saved matrices from memory.
   */
   void clear_saved_matrices();

   /*
      This method computes the \mathcal{A}-integrals (see Erkelenz). The method takes an
      array of v_alpha-values (see Erkelenz) to speed up the computation.
   */
   double compute_A_integral(double qi, double qo, int J,int l, std::vector<double> v_alpha_arr);
   
   /*
      This method performs a partial-wave-decomposition (partial-wave-analysis) of the 
      potential terms v_alpha funciton. This is done by first computing the A-intagrals
      using compute_A_integral(), and then computing the matrix elemetnts of the form

      V_arr = [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
      where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc

      following the formulas in Erkelenz for a given spin structure. (NOTE some formulas are wrong
      is is pointed out in an appendix in one of Machleidts papers)
   */
   void pwa(double qi,double qo, bool coupled, int J,double A_m,double A_p,double A_0,double A_1,std::string spin_struct,bool isovector,double* V_arr);
   
   
   double pot_OPEP_mom(double qo,double qi, double z);
public:
   
   /*
      Some class variables that ideally should be private with getters and setters, but in the
      interest of computation time they are public for the time beeing. (The user of this 
      code is trusted!)
   */
   std::vector<std::string> LEC_names_; // List of ALL possible LEC names
   std::unordered_map<std::string, double> LECs_; // List of lecs  and their values
   std::vector<std::string> LECs_in_use_; // List of lec names of lecs in use in this potential
   std::vector<Term> terms_in_pot_; // Terms in the potential

   /* 
      Constructor
   */
   Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA = 96,double* p_grid = nullptr, double* w_grid = nullptr, std::size_t grid_size = 0,unsigned int J_max = 0, double cutoff_Lambda = 450.0);

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
   void calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr);
   
   // TODO: implement
   double calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz);
 
   /* 
      This function returns a mom_grid_size_ + 1 x mom_grid_size_ + 1 matrix IF the channel is uncoupled
      and 2*mom_grid_size_ + 2 x 2*mom_grid_size_ +2 if the channel is coupled.

      The last momentum point (in the respective blocks in the coupled case) is the on-shell
      q_on_shell passed to the function.
   */
   gsl_matrix* get_matrix(double q_on_shell,qs::quantum_channel chn, bool rel_correction);

   /*
      This function returns a mom_grid_size_ x mom_grid_size_  matrix IF the channel is uncoupled
      and 2*mom_grid_size_ x 2*mom_grid_size_  if the channel is coupled. You probaly dont want 
      relativistic factors here!!

      This function is contrary to get_matrix() useful when solving the Schrödinger equation
      for the given potential.
   */
   gsl_matrix* get_matrix_no_onshell(qs::quantum_channel chn, bool rel_correction);

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
   gsl_matrix* get_saved_matrix(double q_on_shell, qs::quantum_channel chn, bool rel_correction);
 
 };

#endif
