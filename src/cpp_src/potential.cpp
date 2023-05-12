#include "potential.h"


template potential<gsl_matrix>::potential();
template potential<gsl_matrix_complex>::potential();

template potential<gsl_matrix>::~potential();
template potential<gsl_matrix_complex>::~potential();


template <class gsl_m>
potential<gsl_m>::potential()
{}

template <class gsl_m>
potential<gsl_m>::~potential()
{}




