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
 * <p'|p> = (pi/2)*\delta(p'-p)/p^2 (as in Landau).
 */
eigen_t solve_SE(double* p, double* w, unsigned int numer_of_grid_points,qs::quantum_channel chn, const gsl_matrix* V);

/*
 * Converts radians to degrees
 */
double rad_to_deg(double in);

/*
 * Functions to print gsl matrices
 */

void print_m(gsl_matrix* matrix);
void print_v(gsl_vector* vec);

void print_m_complex(gsl_matrix_complex* matrix);
void print_v_complex(gsl_vector_complex* vec);
}
#endif


