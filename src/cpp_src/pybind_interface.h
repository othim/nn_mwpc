/* pybind_interface.cpp
 * This file is a part of the nn_mwpc project.
 *
 * This file defines a class and some functions that will be acessable from
 * python. The class will serves as a simplified interface to the code.
 *
 * Oliver Thim 2021-12-
 * Department of Physics, Chalmers
 */
#ifndef PYBIND11_INTERFACE
#define PYBIND11_INTERFACE
// includes
#include "pybind11/pybind11.h"
#include <algorithm>
#include <iostream>
#include <ctime>

#include "pot_nn_mwpc.h"
#include "quantum_states.h"
#include "LS_Solver.h"
#include "gsl_sf_legendre.h" 
#include "gsl_integration.h" 
#include "scattering.h"
#include "physics_helpers.h"
#include "pot_ext.h"

/* This class will be acessed from python through the bindings in pybind11.
 * The pybind code will be written in sucha a way that C++ will always have
 * ownership over the objects that are created even if it is created in 
 * python.
 *
 * By compiling into the package nn-mwpc and importing the code in
 * python with 'import nn-mwpc' all the functionality will be reached
 * by eg 'myobj = nn-mwpc.nn_mwpc_interface(<args>)' and an object is 
 * created from the class. Then one can call 
 * 'phases = myobj.compute_phase_shifts(<args>)' in python and get phase 
 * shifts.
 *
 * You can't access all the functionality from the code without modifying the
 * code in this class. The main idea is to be able to easily acess key 
 * pre-determined features from python, that are explicitly implemented here.
 * 
 * The initialization of the object can take some time since pre-computations 
 * are made. The idea is that the functions that will be called multiple times
 * should be as fast as possible.
 */
class nn_mwpc_interface
{
private:
/* Here variables that needs to be initialized are stored
 */
    int J_;
public:
    nn_mwpc_interface(int J);
    //~nn_mwpc_interface();

    double* compute_phase_shifts();
    
    double compute_obs();

    //std::string print_settings();
};
/* Here is the code for the pybind11 interface
 */

namespace py = pybind11;


PYBIND11_MODULE(nn_mwpc, m) 
{
    py::class_<nn_mwpc_interface>(m,"nn_mwpc_interface")
        .def(py::init<int>())
        .def("compute_obs", &nn_mwpc_interface::compute_obs);
}
#endif
