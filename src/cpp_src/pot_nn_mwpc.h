/*
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
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" // GL integration
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
   // By calling precompute() the matrices are stored in memory.
   // You can ask the class which channels are stored in memory
   std::map<qs::quantum_channel, std::map<std::string, gsl_matrix*>, qs::comp> saved_matrices_; // List of matrices that are saved by the class
   
   // Variables to store the gauss legendre grid for the potential
   double* p_grid_;
   double* w_grid_;
   std::size_t mom_grid_size_;
   double cutoff_Lambda_;

   unsigned int J_max_;
   double** stored_Legendre_polynomials_;
   // Variables to store
   unsigned int N_GLI_PWA_;
   gsl_integration_fixed_workspace* int_ang_;
   double* z_mesh;
   double* w_z_mesh;
   unsigned int len_z_mesh;
private:

   /*
      opep_get_el
      V(...) returns a list of length 6
      The elements are on the form [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
      where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc
   */
   void clear_saved_matrices();

   double compute_A_integral(double qi, double qo, int J,int l, std::vector<double> v_alpha_arr);
   void pwa(double qi,double qo, bool coupled, int J,double A_m,double A_p,double A_0,double A_1,std::string spin_struct,bool isovector,double* V_arr);
   double pot_OPEP_mom(double qo,double qi, double z);
public:
   
   std::vector<std::string> LEC_names_;
   std::unordered_map<std::string, double> LECs_; // List of lecs  and their values. Public for now...
   std::vector<std::string> LECs_in_use_; // List of lec names of lecs in use in this potential
   std::vector<Term> terms_in_pot_; // Terms in the potential

   // Constructor
   Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA = 96,double* p_grid = nullptr, double* w_grid = nullptr, std::size_t grid_size = 0,unsigned int J_max = 0);

   // Destructor 
   ~Potential_mwpc();

   void calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr);
   double calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz);
 
   /* 
      This function returns a mom_grid_size_ x mom_grid_size_ matrix IF the channel is uncoupled
      and 2*mom_grid_size_ x 2*mom_grid_size_ if the channel is coupled.
   */
   gsl_matrix* get_matrix(double q_on_shell,qs::quantum_channel chn, bool rel_correction);
   void populate_saved_mtx(qs::quantum_channel chn, bool rel_correction);
   gsl_matrix* get_saved_matrix(double q_on_shell, qs::quantum_channel chn, bool rel_correction);
 /*
      Returns a list of potential elements given the lecs. In some terms the LECs
      do not enter like \alpha_i W_i and therefor the matrix element depends on
      the vales of the lecs. For OPE and the LO contact terms this is not the
      case.
   */
   //std::vector<double> get_mtx_el(Term dia, double pp, double p, Channel chn);

   /*
      This function will popolate the W_list by computing the potential for
      one neq 0 lec at the time. This means that the full potential mastrix can
      be obtained by a lecs * W.
   */
   //int populate(bool include_rel_fac, bool include_cutoff);

   /*
      This function geives the whole potential matrix for a geven set of lecs.
      This function call will be fast if the W's are precomputed and saved in
      memory during runtime.
   */

   //double** get_mtx(LECs lecs);
};

#endif