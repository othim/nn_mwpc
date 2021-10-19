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

std::vector<qs::quantum_channel> get_channels(std::vector<qs::quantum_NN_state> states, bool print)
{
    // Conserved quantum numbers are J,S,T,Tz,pi. 
 
    std::map<q_chn,std::vector<in_out_state>,comp_chn> channels; // To not have duplicates

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
                    .pi = bra.pi,
                };
                
                if (channels.count(chn) == 0) // Key is not present
                {
                    std::vector<in_out_state> state;
                    state.push_back(ios);
                    channels.insert(std::make_pair (chn,state));
                } else // If key is present, add this state to the coupled channel
                {
                    channels[chn].push_back(ios);
                }
            }
        }
    }

    // If print
    // Just want the vector of quantum_channels
    std::vector<qs::quantum_channel> chns;

    if (print)
    {
        std::cout << std::endl << "Conserved quantum numbers: J,S,T,Tz,pi" << std::endl;
        int j = 0;
        for (auto it = channels.begin(); it != channels.end();++it)
        {
            q_chn key = (*it).first;
            bool coup = false;
            if ((*it).second.size() >1) {
                coup = true;
            }
            
            // Chreate a channel and append it to the list of channels
            qs::quantum_channel chn = {.J = key.J, .S = key.S, .tz = key.Tz, coupled = coup};
            chns.push_back(chn);

            std::cout << "Channel " << j << ": J=" << key.J << " S=" << key.S << " T=" << key.T << 
                " Tz=" << key.Tz << " pi=" << key.pi << " coup=" << coup << std::endl;

            for (int i=0; i < (*it).second.size(); i++)
            {
                in_out_state sm = (*it).second[i];
                std::cout << "   (J=" << sm.J << " Li=" << sm.Li << " Lo=" << sm.Lo 
                    << " S=" << sm.S << " T=" << sm.T << 
                        " Tz=" << sm.Tz << " pi=" << sm.pi << ")" << std::endl;
            }
            std::cout << std::endl;
            j++;
        }
         std::cout << "Done printing" << std::endl;
    }
    // Return the list of qs::quantum_channel
    return chns;
}