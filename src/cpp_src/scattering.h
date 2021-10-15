/*
    scattering.h
    This file is part of the nn_mwpc project.

    This header contains functions to compute NN scattering
    observalbes from phase shifts. The conventions used for 
    phase shifts, M/K-matrix elements and observables should
    be stated in connection to the function.

    Oliver Thim 2021-09 --
    Department of Physics, Chalmers
*/
#include "quantum_states.h"
#include "gsl_sf_legendre.h"
#include "wigxjpf.h"
#include <vector>
#include <complex>
#include <iostream>


/*
    This function computes the M-matrix elements from phase shifts.
    (The M/K-matrix is basically just a matrix of scattering amplitudes)

    The exact formula that is implemented is stated in the README file.
    Generally the conventions in the book by Taylor is used here.
    Matrix can generally be complex, and phase shifts can either be
    real or complex. This function is written for REAL phase shifts.

    The M-matrix is defined as a matrix of the scattering amplitudes 
    M(\theta,\phi) := { f(p' lambda' <- p,lambda) }_lambda lambda' where lambda lambda' 
    are general spin projections. The amplitude, f, is related to the
    T-matrix as f = -(2 pi)^2 \mu T(on shell). See Taylor.
*/
std::complex<double> get_M_matrix_p(std::vector<qs::quantum_channel> chns_vec,std::vector<Phase_shifts_chn> phase_shifts_vec, unsigned int s, int mo, int mi, double cos_theta, double q_on_shell,double rho_T,unsigned int l_max);

/*
    This function is the same as get_M_matrix_p but computes the M-matrix directly from 
    the T-matrix elements.

    The T-matrix elements are in the form [Tmm, Tmp, Tpp] if coupled and if uncoupled only Tmm
    is non-zero.
*/
std::complex<double> get_M_matrix_T(std::vector<qs::quantum_channel> chns_vec,std::vector<std::complex<double>*>T_on_shell_vec, double q_on_shell, unsigned int s, int mo, int mi, double cos_theta,unsigned int l_max);

double compute_Saclay_amplitudes();

double compute_observebles();