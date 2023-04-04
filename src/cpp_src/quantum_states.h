/*
    This file contains some data structures to easily handle
    and construct quantum states in NN systems.

    Functions that are handy for constructing states are implemented in the 
    quantum_states.cpp file.
*/
#if !defined(QUANTUM_STATES_H)
#define QUANTUM_STATES_H 1

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <utility>
#include <cstdlib>
#include <string>
// This is the possible quantum channels if [T,H]=0.
namespace qs 
{
    struct quantum_channel 
    {
         int J; int S; int T; int Tz; bool coupled;
    };

    struct comp
    {
        bool operator() (quantum_channel a, quantum_channel b) const 
        {
            return std::make_tuple(a.J,a.S,a.Tz,a.coupled) < 
                std::make_tuple(b.J,b.S,b.Tz,b.coupled);
        }
    };

    struct quantum_NN_state
    {
        int J; int L; int S; int T; int Tz; int pi;

    };

};

struct Phase_shifts_chn {double delta_p; double delta_m; double epsilon; 
    double delta_uncoupled;};

struct lo_li {int lo; int li;};

// For key in map of in_out_state
struct q_chn {int J; int S; int T; int Tz; int pi;};

struct comp_chn
{
    bool operator() (q_chn a, q_chn b) const 
    {
        return std::make_tuple(a.J,a.S, a.T,a.Tz,a.pi) < 
            std::make_tuple(b.J,b.S, b.T,b.Tz,b.pi);
    }
};

// Possible bra ket state
struct in_out_state { int J; int Li; int Lo; int S; int T; int Tz; int pi; };

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
    The Puli principle gives that (-1)^{s+l+t} = -1
*/
std::vector<qs::quantum_NN_state> get_states_NN(int J_max, int J_min, 
        int Tz_min, int Tz_max, bool print);

/*
*   This function groups the states into quantum channels 
*/
std::vector<qs::quantum_channel> get_channels(
        std::vector<qs::quantum_NN_state> states,bool print);

std::string quantum_channel_to_string(qs::quantum_channel chn); 

#endif
