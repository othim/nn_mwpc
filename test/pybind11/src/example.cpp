#include <pybind11/pybind11.h>
#include "header.h"


struct data
{
    double d;
    int i;
};

static data d;

void init(int k)
{
    d = {.d=2.0, .i=4};
}

int add(int i, int j) {
    return d.d + d.i;
}


int sub_from_h(int i, int j)
{
    return sub(i,j);
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    
    m.def("init",&init, " ");
    m.def("add", &add, "A function which adds two numbers");
    m.def("sub_from_h", &sub_from_h, "Funtion from other file");
}

// To compile (on MacOS):
// g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC $(python3-config --includes) -I/opt/anaconda3/envs/nn-mwpc-env/include example.cpp -o example.so   
// The undefined dynamic_lookup are not needed on linux
