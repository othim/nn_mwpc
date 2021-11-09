#ifndef PHYSICS_HELPERS_H
#define PHYSICS_HELPERS_H

#include <iostream>
#include <cmath>
#include "gsl_matrix.h"
#include "gsl_complex.h"
#include "gsl_complex_math.h"
#include "gsl_matrix_complex_double.h"
#include "wigxjpf.h"
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

#endif


