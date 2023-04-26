#ifndef PHYSICS_HELPERS_H
#define PHYSICS_HELPERS_H

#include <iostream>
#include <cmath>
#include "gsl_matrix.h"
#include "gsl_integration.h"
#include "gsl_sf_legendre.h"
#include "gsl_complex.h"
#include "gsl_complex_math.h"
#include "gsl_matrix_complex_double.h"
#include "wigxjpf.h"
#include "quantum_states.h"
#include "Constants.h"
#include "gsl_eigen.h"
#include "gsl_blas.h"
#include <stdio.h>

namespace ph {
/*
 * This function initialize relevant data that needs to be in memory 
 * to be able to use the functions.
 */
void physics_helpers_init();

/*
 * This function should be called when the use of these functions
 * are done to clean up some saved data.
 */
void physics_helpers_free();


// Get Pauli matrices

/*
 * Returns the pauli matrix, x,y or z depending on the 
 * value of axis='x','y' or 'z'.
 */
gsl_matrix_complex* get_Pauli_matrix(char axis); 


/* Construct kronecker product of matrices, where the matrix m2 is placed
 * at each position in the matrix m1.
 */ 
gsl_matrix_complex* kronecker_product(gsl_matrix_complex* m1, gsl_matrix_complex* m2);


// Take trace
gsl_complex trace(gsl_matrix_complex* m);


// Get GC-coefficients unsing wigxjpf
double CG_coeff(int J_2, int M_2, int j1_2, int j2_2, int m1_2, int m2_2);

/*
 * This function fills the passed arrays with Gauss-Legendre quadratur points
 * for the intervall [0,infty) with the given scale.
 */
void gauss_legendre_inf_mesh(unsigned int Number_of_points, double scale,double** p,double** w);

/*
 * This function fills the passed arrays with Gauss-Legendre quadratur points
 * for the intervall [min,max]
 */
void gauss_legendre_finite_mesh(unsigned int Number_of_points, double min,
        double max,double** p,double** w);


typedef struct 
{
   gsl_vector_complex* eigenvalues;
   gsl_matrix_complex* eigenvectors;
} eigen_t;

/*
 * This function solves the SE for a given parital-wave channel and
 * produces eigenvales and eigenvectors. The method is exact diagonalization 
 * of the Hamiltonian in the given quadreture basis.
 * The potential is assumed to be in a partial wave basis with normalization 
 * <p'|p> = \delta(p'-p)/p^2, so the factor \pi/2 from Landau is removed.
 * Furthermore the 3-momentum normalization is assumed to be
 * <p',p> = \delta^3(p'-p) without a (2\pi)^3 factor, which means that a potential
 * that is written in the above partial wave basis might also need factor of
 * (2\pi)^(+/-)3 to convert to the correct basis.
 */
eigen_t solve_SE(double* p, double* w, unsigned int numer_of_grid_points,
        qs::quantum_channel chn, const gsl_matrix* V);

/*
 * Converts radians to degrees
 */
double rad_to_deg(double in);

/*
 * Functions to print gsl matrices
 */

void print_m(gsl_matrix* matrix);
void print_v(gsl_vector* vec);

void print_m(gsl_matrix_complex* matrix);
void print_v(gsl_vector_complex* vec);

void print_m_complex_to_file(std::string path,gsl_matrix_complex* matrix);

/* 
 *
 * Matrix operations
 *
 *
 */

/*
 * This function implements the onshell multiplication, meaning 
 * that the matrix multiplixation does not include the last 
 * column/row of the matrix. This is achieved by first multiplying
 * the matrices as usual and then subtracting the error that is 
 * induced. The time lost is negligable compared to just doing
 * an ordinary multiplication.
 *
 * res = m1*m2
 * !!NOTE!!
 * m1, m2, res needs to be distinct matrices. You cannot have eg.
 * m1 <- m1*m2
 *
 * This function is tested againts on_shell_mult_bf() which is a loop
 * brute force verion of the original sum that we want to compute.
 */
void on_shell_mult(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res);
/*
 * Ordinary matrix multiplization
 */
void mult(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res);
/*
 * This is the brute force version of the 'ph::on_shell_mult()' function
 */
void on_shell_mult_bf(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res);

/*
 * Function that copies the content on the real matrix into the 
 * complex matrix.
 */
void make_matrix_complex(gsl_matrix_complex* Mz,gsl_matrix* M);

void matrix_set(gsl_matrix* m,         size_t i, size_t j, double el);
void matrix_set(gsl_matrix_complex* m, const size_t i, const size_t j, 
        const gsl_complex el);
void matrix_set(gsl_matrix_complex* m, const size_t i, const size_t j, 
        const double el);

/*
 * This function assumes that el is the real part if the matrix element
 */
void matrix_set(gsl_matrix_complex* m, const size_t i, const size_t j, 
        const double el);

gsl_matrix*         matrix_alloc(const size_t i, const size_t j, gsl_matrix* m);
gsl_matrix_complex* matrix_alloc(const size_t i, const size_t j, gsl_matrix_complex* m);

void matrix_free(gsl_matrix* m);
void matrix_free(gsl_matrix_complex* m);

void matrix_set_zero(gsl_matrix* m);
void matrix_set_zero(gsl_matrix_complex* m);

void matrix_add(gsl_matrix* m1, gsl_matrix* m2);
void matrix_add(gsl_matrix_complex* m1, gsl_matrix_complex* m2);

void matrix_sub(gsl_matrix* m1, gsl_matrix* m2);
void matrix_sub(gsl_matrix_complex* m1, gsl_matrix_complex* m2);

void matrix_memcpy(gsl_matrix* m1, gsl_matrix* m2);
void matrix_memcpy(gsl_matrix_complex* m1, gsl_matrix_complex* m2);

void matrix_scale(gsl_matrix* m1, const double scale);
void matrix_scale(gsl_matrix_complex* m1, const double scale);

void matrix_from_vector(gsl_matrix_complex* M,gsl_vector_complex* vec);

double get_mN(int Tz);
}
#endif


