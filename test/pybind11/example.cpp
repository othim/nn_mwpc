#include <pybind11/pybind11.h>

struct data
{
    double d;
    int i;
};

static data d;

void init(int k)
{
    d = {.d=3.0, .i=4};
}

int add(int i, int j) {
    return d.d + d.i;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    
    m.def("init",&init, " ");
    m.def("add", &add, "A function which adds two numbers");
}

// To compile (on MacOS):
// g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC $(python3-config --includes) -I/opt/anaconda3/envs/nn-mwpc-env/include example.cpp -o examp│toliver@elan.local /opt/anaconda3/envs/nn-mwpc-env 
// le.so   
