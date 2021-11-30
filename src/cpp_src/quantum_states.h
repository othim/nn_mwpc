/*
    This file contains some data structures to easily handle
    and construct quantum states in NN systems.

    Functions that are handy for constructing states are implemented in the 
    quantum_states.cpp file.
*/
#if !defined(QUANTUM_STATES_H)
#define QUANTUM_STATES_H 1

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <utility>
#include <cstdlib>

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
            return std::make_tuple(a.J,a.S,a.tz,a.coupled) < std::make_tuple(b.J,b.S,b.tz,b.coupled);
        }
    };

    struct quantum_NN_state
    {
        unsigned int J; unsigned int L; unsigned int S; unsigned int T; int Tz; int pi;
    };
};

struct Phase_shifts_chn {double delta_p; double delta_m; double epsilon; double delta_uncoupled;};

struct lo_li {unsigned int lo; unsigned int li;};

// For key in map of in_out_state
struct q_chn {unsigned int J; unsigned int S; unsigned int T; int Tz; int pi;};

struct comp_chn
{
    bool operator() (q_chn a, q_chn b) const 
    {
        return std::make_tuple(a.J,a.S, a.T,a.Tz,a.pi) < std::make_tuple(b.J,b.S, b.T,b.Tz,b.pi);
    }
};

// Possible bra ket state
struct in_out_state 
{ unsigned int J; unsigned int Li; unsigned int Lo; unsigned int S; 
    unsigned int T; int Tz; int pi;    
};

/*
*   ---------------------------------------------------------------
*   Functions that are used to manipulate these data structures and 
*   construct quantum states.
*   ---------------------------------------------------------------
*/
/*
    This function creates the possible discrete quantum states for the NN-system
    that is expressed in the basis
    |p> x |j,l,s,mj,t,tz>, where the last part of the state is the discrete part.
    The Puli principle gives that (-1)^{j+l+t} = -1
*/
std::vector<qs::quantum_NN_state> get_states_NN(unsigned int J_max, unsigned int J_min, int Tz_min,
    int Tz_max, bool print);

/*
*   This function groups the states into quantum channels 
*/
std::vector<qs::quantum_channel> get_channels(std::vector<qs::quantum_NN_state> states,bool print);


#endif
