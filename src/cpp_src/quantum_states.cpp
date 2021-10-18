#include "quantum_states.h"


std::vector<qs::quantum_NN_state> get_states_NN(unsigned int J_max, unsigned int J_min, int Tz_min,
    int Tz_max, bool print)
{
    std::vector<qs::quantum_NN_state> states;

    for (int Tz = Tz_min; Tz < Tz_max+1; Tz ++)
    {
        for (unsigned int J = J_min; J < (J_max + 1); J++)
        {
            for (unsigned int S = 0; S < 2; S++)
            {
                for (unsigned int L = std::abs(int(J-S)); L < J+S+1; L++)
                {
                    for (unsigned int T = std::abs(Tz); T < 2; T++)
                    {
                        if ( (L+S+T) % 2 != 0) // Pauli principle
                        {
                            qs::quantum_NN_state qs = 
                            {.J = J,
                            .L=L,
                            .S=S,
                            .T=T,
                            .Tz=Tz,
                            .pi = (int)std::pow(-1,L),
                            };
                            if (print) {
                                std::cout << "State: " << "J=" << J << " L=" << L << " S=" << S << " T="
                                    << T << " Tz=" << Tz << " pi=" << qs.pi << std::endl;
                            }
                            states.push_back(qs);
                        }
                    }
                }
            }
        }
    }
    return states;
}

bool comp_in_out_states(const in_out_state &f, const in_out_state &s)
{
    if (f.S < s.S) {
        return true;
    } else {
        return false;
    } 
    /*
    else if (f.S == s.S && f.J < s.J)
    {
        return true;
    } else if (f.S == s.S && f.J == s.J && f.pi < s.pi)
    {
        return true;
    } else 
    {
        return false;
    }     */
}

std::vector<qs::quantum_channel> get_channels(std::vector<qs::quantum_NN_state> states, bool print)
{
    // Conserved quantum numbers are J,S,T,Tz,pi. 
    // For now these are hard coded in


    // Loop through the bra and ket states

    // If quantum numbers conserved save that to the list of channels
    // Increase a counter so that we know if its coupled or not
    // Get a list of channels, not necessarily in the same order as in the
    // Python code, but it does not matter 

    std::map<q_chn,std::vector<in_out_state>,comp_chn> states_in_out; // To not have duplicates
    //std::unordered_map<double, int> states_in_out;

    for (std::size_t i = 0; i < states.size(); i++)
    {
        qs::quantum_NN_state bra = states[i];

        for (std::size_t j = 0; j < states.size(); j++)
        {
            qs::quantum_NN_state ket = states[j];

            if (bra.J == ket.J && bra.S == ket.S && bra.T == ket.T && bra.Tz == ket.Tz &&
                bra.pi == ket.pi)
            {
                // Insert this state in the correct channel that is only created once
                in_out_state ios= 
                {
                    .J = bra.J,
                    .Li=ket.L,
                    .Lo=bra.L,
                    .S = bra.S,
                    .T = bra.T,
                    .Tz = bra.Tz,
                    .pi = bra.pi
                };

                q_chn chn = 
                {
                    .J = bra.J,
                    .S = bra.S,
                    .T = bra.T,
                    .Tz = bra.Tz,
                    .pi = bra.pi
                };
                
                if (states_in_out.count(chn) == 0) // Key is not present
                {
                    std::vector<in_out_state> state;
                    state.push_back(ios);
                    states_in_out.insert(std::make_pair (chn,state));
                } else
                {
                    states_in_out[chn].push_back(ios);
                }
            }
        }
    }

   // std::vector<in_out_state> states_vec(states_in_out.begin(),states_in_out.end());
    // If print
    if (print)
    {
        std::cout << std::endl;
        for (auto it = states_in_out.begin(); it != states_in_out.end();++it)
        {
            for (int i=0; i < (*it).second.size(); i++)
            {
                in_out_state sm = (*it).second[i];
                std::cout << "(J=" << sm.J << " Li=" << sm.Li << " Lo=" << sm.Lo 
                    << " S=" << sm.S << " T=" << sm.T << 
                        " Tz=" << sm.Tz << " pi=" << sm.pi << ")" << std::endl;
            }
            std::cout << std::endl;
        }
    }
    std::cout << "Done printing" << std::endl;
    // By this point all possible in and out states
    // are identified and they now needs to be grouped
    // by the conserved quantum numbers. Maybe use the library
    // #include <range/v3/all.hpp>

    // Sorting can be done with the container std::multimap
    // with the appropriate compare instruction.

    // Sorting: S, J, pi mm,mp,pm,pp in l
    
    std::vector<qs::quantum_channel> chns;
    /*std::vector<std::vector<qs::quantum_channel> > blocks;
    in_out_state lp = states[0];
    unsigned int counter = 0;

    std::vector<qs::quantum_channel> block;

    for (std::size_t i = 1; i < states.size(); i++) 
    {
        std::pair<in_out_state, int> cp = states[i];

        // If this is true the states are in the same channel.
        if (cp.first.J == lp.first.J && cp.first.S == lp.first.S && 
            cp.first.Tz == lp.first.Tz && cp.first.pi == lp.first.pi)
        {
            counter++;
            block.push_back(cp.first); // Add the current state to the block
        } else 
        {
            block.push_back(lp.first); // Add the current state to the block
            coup = false;
            if (counter>2) : coup = true;
            qs::quantum_channel new_chn = {.J = lp.first.J,.S = lp.first.S,lp.first.Tz,.coupled = coup};
            chns.push_back(new_chn);
            
            if (print)
            {
                // Print info on the created channel
                std::cout << "channel" << (chn.size()-1) << ": J=" << new_chn.J << " S=" << new_chn.S 
                    << " Tz=" << new_chn.tz << " Pi=" << lp.first.pi << std::endl;
                
                // Print the channels
                for (int i = 0; i < chns.size(); i++)
                {
                    std::cout << "   Lo=" << chns[i].Li << " Li=" << chns[i].Lo << std::endl;
                }                
            }
            blocks.push_back(block);
            block.clear();
            counter = 0;
            lp = cp;
        }
          if (states_in_out.size()==0)
                {
                    states_in_out.push_back(ios);
                } else 
                {
                    bool inserted = false;
                    for (auto it = states_in_out.begin(); it != states_in_out.end(); ++it)
                    {
                        if (comp_in_out_states(ios,*it)) // iso < states[j] insert before
                        {
                            std::cout << comp_in_out_states(ios,*it) << " " << ios.S << "<" << (*it).S << std::endl; 
                            states_in_out.insert(it,ios);
                            inserted = true;
                        }
                    }
                    if (inserted==false) 
                    {
                        states_in_out.push_back(ios);
                    }
    } */
    return chns;
}