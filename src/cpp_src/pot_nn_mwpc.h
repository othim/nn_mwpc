/*
pot_nn_mwpc.h

This file is part of the nn_mwpc project.

Oliver Thim 2021-09 --
Department of Physics, Chalmers
*/
#include <vector>


class Potential_mwpc
{
private:
   std::vector<double**> W_list; // List of matrices that are saved by the class
   LECs lecs; // List of lecs
   std::vector<Terms> terms; // Terms in the potential
   double* p_grid;
   double* w_grid;
   std::vector<double> phys_constants;
private:

public:
   Potential_mwpc();

   // Returns a list of potential elements given the lecs. In some terms the LECs
   // do not enter like \alpha_i W_i and therefor the matrix element depends on
   // the vales of the lecs. For OPE and the LO contact terms this is not the
   // case.
   std::vector<double> get_mtx_el(Term dia, double pp, double p, Channel chn);

   /*
      This function will popolate the W_list by computing the potential for
      one neq 0 lec at the time. This means that the full potential mastrix can
      be obtained by a lecs * W.
   */
   int populate();

   /*
      This function geives the whole potential matrix for a geven set of lecs.
      This function call will be fast if the W's are precomputed and saved in
      memory during runtime.
   */
   double** get_mtx(LECs lecs);
};
