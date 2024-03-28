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
#include "physics_helpers.h"
#include "gsl_sf_legendre.h"
#include "gsl_blas.h"
#include "wigxjpf.h"
#include <vector>
#include <complex>
#include <iostream>
#include "gsl_matrix.h"

namespace sc 
{
    

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
std::complex<double> get_M_matrix_p(std::vector<qs::quantum_channel> chns_vec,
        std::vector<Phase_shifts_chn> phase_shifts_vec, int s, int mo, int mi, 
        double cos_theta, double q_on_shell,double rho_T,int l_max);

/*
    This function is the same as get_M_matrix_p but computes the M-matrix directly from 
    the T-matrix elements.

    The T-matrix elements are in the form [Tmm, Tmp, Tpp] if coupled and if uncoupled only Tmm
    is non-zero.
*/
std::complex<double> get_M_matrix_T(std::vector<qs::quantum_channel> chns_vec,
        std::vector<std::complex<double>*> T_on_shell_vec, double q_on_shell, 
        int s, int mo, int mi, double cos_theta,int l_max);


/*
    This function computes the Saclay amplitudes as defined in:
    Formalism of nucleon-nucleon elastic scattering experiments. 
    Journal de Physique, 1978, 39 (1), pp.1-32.
    This is done by first computing all relevant M-matrix elements and then 
    forming the amplitudes as certain combinations of them.
*/
std::vector<std::complex<double> > compute_Saclay_amplitudes(
        std::vector<qs::quantum_channel> chns_vec, 
        std::vector<Phase_shifts_chn> phase_shifts_vec, 
        double theta, double q_on_shell,double rho_T, int l_max,
        ph::constants_struct* program_const);

/*
 * Same function as above but the T-matrix elements is an argument
 * instead of the phase shifts.
 */
std::vector<std::complex<double> > compute_Saclay_amplitudes(
        std::vector<qs::quantum_channel> chns_vec, 
        std::vector<std::complex<double>*> T_on_shell_vec, 
        double theta, double q_on_shell,double rho_T, int l_max,
        ph::constants_struct* program_const);

/*
 * Helper function for computing the saclay amplitudes given
 * the different M-matrix elements
 */
std::vector<std::complex<double>> saclay_amplitudes_from_M_elements(
        std::complex<double> M_pp, std::complex<double> M_00,
        std::complex<double> M_pm, std::complex<double> M_s,
        std::complex<double> M_p0, std::complex<double> M_0p,double theta,
        double MeVm2_to_mbarn);

/*
    This function computes observables from Saclay amplitudes.
    The observables are specified as follows.

    Eg I_0000 <-> obs="I 0000"
       C_nn00 <-> obs="C nn00"
       C_llll <-> obs="C llll"
       etc.
    Total cross sections also:
        obs="SGT"
        obs="SGTT"
        obs="SGTL"
    This follows the notation in Table III in 
    Formalism of nucleon-nucleon elastic scattering experiments. 
    Journal de Physique, 1978, 39 (1), pp.1-32.
    
    Here:
        l \propto k_i + k_f
        m \propto k_f - k_i
        n \propto k_i x k_f

    Other notations for these vectors can be:

    q_on_shell in the on shell momentum in MeV
*/
double compute_observable(std::vector<std::complex<double> > sac_amp, 
        double q_on_shell, std::string obs, ph::constants_struct* program_const);

/*
 * This function computes the observables that are defined in the lab frame.
 * The formalisam follows
 *
    Formalism of nucleon-nucleon elastic scattering experiments. 
    Journal de Physique, 1978, 39 (1), pp.1-32.
 */
double compute_observable_lab(std::vector<std::complex<double> > sac_amp, 
        double q_on_shell, std::string obs, double angle, 
        ph::constants_struct* program_const);
/*
 * This function computes the total cross section SGT. The return unit is mb.
 */
double compute_total_cross_section(std::vector<qs::quantum_channel> chns_vec, 
    std::vector<std::complex<double>*> T_on_shell_vec,double q_on_shell,
    int l_max,bool optical_thm,
    ph::constants_struct* program_const);

/* 
 * This function gives the M-matrix for the given on shell energy
 * and cm scattering angle
 */ 
gsl_matrix_complex* get_M_matrix(std::vector<qs::quantum_channel> chns_vec,
    std::vector<Phase_shifts_chn> phase_shifts_vec, double q_on_shell, 
    double theta, double rho_T, int l_max);

/*
 * This function computes the observable trace as in 
 * Formalism of nucleon-nucleon elastic scattering experiments. 
 * Journal de Physique, 1978, 39 (1), pp.1-32.
 */ 
 double get_observables(gsl_matrix_complex* sigma_i_1, gsl_matrix_complex* sigma_i_2,
        gsl_matrix_complex* sigma_o_1, gsl_matrix_complex* sigma_o_2,
        gsl_matrix_complex* M_matrix);

/*
 * Funtion to compute the S-matrix from Stapp phase shifts
 *
 * The function returns the 11, 12 and 22 element of the S-matrix
 */
std::complex<double>* S_from_Stapp(double dm, double dp, double eps);

/*
 * Funtion to compute the S-matrix from BB phase shifts.
 *
 * The function returns the 11, 12 and 22 element of the S-matrix
 */
std::complex<double>* S_from_BB(double dm, double dp, double eps);


std::vector<std::complex<double>*> T_from_phase_shifts(std::vector<Phase_shifts_chn> phases, std::vector<qs::quantum_channel> chns_vec,double rho_T);

}
