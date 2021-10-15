/*
    This file contains some data structures to easily handle
    and construct quantum states in NN systems.
*/
#if !defined(QUANTUM_STATES_H)
#define QUANTUM_STATES_H 1

// This is the possible quantum channels if [T,H]=0.
namespace qs 
{
    struct quantum_channel 
    {
        unsigned int J; unsigned int S; int tz; bool coupled;
    };

    struct comp
    {
        bool operator() (quantum_channel a, quantum_channel b) const 
        {
            return a.J < b.J;
        }
    };
};

struct Phase_shifts_chn {double delta_p; double delta_m; double epsilon; double delta_uncoupled;};

struct lo_li {int lo; int li;};

#endif