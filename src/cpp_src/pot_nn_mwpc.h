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
#include <vector>
#include "gsl_sf_legendre.h" // Legendre polynomials
#include "gsl_integration.h" // GL integration
#include "Constants.h"
#include <iostream>
#include <cmath>

class Potential_mwpc
{
private:
   //std::vector<double**> W_list; // List of matrices that are saved by the class
   //LECs lecs; // List of lecs
   //std::vector<Terms> terms; // Terms in the potential
   //double* p_grid;
   //double* w_grid;
   //std::vector<double> phys_constants;

   static const int N_GLI_PWA = 96;
   gsl_integration_fixed_workspace* int_ang;
public:

   /*
      opep_get_el
      V(...) returns a list of length 6
      The elements are on the form [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
      where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc
   */
   void opep_get_el(double qi,double qo, bool coupled, int J,double* output);
public:
   // Constructor
   Potential_mwpc();

   // Destructor 
   ~Potential_mwpc();
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
