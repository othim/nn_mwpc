#include "quantum_states.h"


std::vector<qs::quantum_NN_state> get_states_NN(int J_max, int J_min, 
        int Tz_min, int Tz_max, bool print)
{
    std::vector<qs::quantum_NN_state> states;

    for (int Tz = Tz_min; Tz < Tz_max+1; Tz ++)
    {
        for (int J = J_min; J < (J_max + 1); J++)
        {
            for (int S = 0; S < 2; S++)
            {
                for (int L = std::abs(int(J-S)); L < J+S+1; L++)
                {
                    for (int T = std::abs(Tz); T < 2; T++)
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
                                std::cout << "State: " << "J=" << J 
                                    << " L=" << L << " S=" << S << " T="
                                    << T << " Tz=" << Tz << " pi=" << qs.pi 
                                    << std::endl;
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


std::vector<qs::quantum_NN_HO_state> get_states_NN_HO(int N_max,
        int Tz_min, int Tz_max, bool print)
{
    std::vector<qs::quantum_NN_HO_state> states;
    for (int Tz = Tz_min; Tz < Tz_max+1; Tz ++)
    {
        for (int n=0; n<N_max/2+1; n++)
        {
            for (int L=0; L<N_max-2*n+1; L++)
            {
                for (int S=0; S<2; S++)
                {
                    for (int T = std::abs(Tz); T < 2; T++)
                    {
                        if ( (L+S+T) % 2 != 0)
                        {
                            for (int J = std::abs(int(L-S)); J < L+S+1; J++)
                            {
                                qs::quantum_NN_HO_state qs = 
                                {
                                .n=n,
                                .L=L,
                                .S=S,
                                .J = J,
                                .T=T,
                                .Tz=Tz,
                                .pi = (int)std::pow(-1,L)
                                };
                                if (print) {
                                    std::cout << "State: " << "n=" << n << " L=" << L 
                                        << " S=" << S << " J=" << J << " T="
                                        << T << " Tz=" << Tz << " pi=" << qs.pi 
                                        << std::endl;
                                }
                                states.push_back(qs);
                            }
                        }
                    }
                }
            }
        }
    }
    return states;
}



std::vector<qs::quantum_channel> get_channels(
        std::vector<qs::quantum_NN_state> states, bool print)
{
    // Conserved quantum numbers are J,S,T,Tz,pi. 
    // To not have duplicates
    std::map<q_chn,std::vector<in_out_state>,comp_chn> channels; 
    for (std::size_t i = 0; i < states.size(); i++)
    {
        qs::quantum_NN_state bra = states[i];

        for (std::size_t j = 0; j < states.size(); j++)
        {
            qs::quantum_NN_state ket = states[j];

            if (bra.J == ket.J && bra.S == ket.S && bra.T == ket.T 
                    && bra.Tz == ket.Tz && bra.pi == ket.pi)
            {
                // Insert this state in the correct channel 
                // that is only created once
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
                } else // If key is present, add this state
                {
                    channels[chn].push_back(ios);
                }
            }
        }
    }

    // If print
    // Just want the vector of quantum_channels
    std::vector<qs::quantum_channel> chns;
    
    if (print) {
        std::cout << std::endl << "Conserved quantum numbers: J,S,T,Tz,pi" 
            << std::endl;
    }
    int j = 0;
    for (auto it = channels.begin(); it != channels.end();++it)
    {
        q_chn key = (*it).first;
        bool coup = false;
        if ((*it).second.size() >1) {
            coup = true;
        }
        
        // Chreate a channel and append it to the list of channels
        qs::quantum_channel chn = {.J = key.J, .S = key.S, .T = key.T, 
            .Tz = key.Tz, .coupled = coup};
        chns.push_back(chn);
        if (print) {

            std::cout << "Channel " << j << ": J=" << key.J << " S=" << key.S 
                << " T=" << key.T << " Tz=" << key.Tz << " pi=" << key.pi << 
                " coup=" << coup << std::endl;
        }
        for (int i=0; i < (int)(*it).second.size(); i++)
        {
            in_out_state sm = (*it).second[i];
            if (print) {
            std::cout << "   (J=" << sm.J << " Li=" << sm.Li << " Lo=" << sm.Lo 
                    << " S=" << sm.S << " T=" << sm.T << 
                    " Tz=" << sm.Tz << " pi=" << sm.pi << ")" << std::endl;
            }
        }
        if (print) {
            std::cout << std::endl;
        }
        j++;
    }
    if (print)
    {
        std::cout << "Done printing" << std::endl;
    }
    // Return the list of qs::quantum_channel
    return chns;
}


std::string quantum_channel_to_string(qs::quantum_channel chn)
{
    char LS[9] = {'S','P','D','F','G','H','I','K','L'};
    if (chn.J > 8)
    {
        std::cout << "J to high in quantum_chennal_to_string" << std::endl;
        
        std::string s;
        return s;
    }
    int s = (int)(2*chn.S + 1);
    std::string s_s = std::to_string(s);
    std::string s_j = std::to_string(chn.J);
    
    // For L
    std::string s_LS;
    if (chn.coupled)
    {
        s_LS = std::string(1,LS[chn.J-1]) + "-" + std::string(1,LS[chn.J+1]);       
    } else if (!chn.coupled && chn.T == 1 && chn.J == 0 && chn.S == 1) // 3P0
    {
        s_LS = "P";
    } else
    {
        s_LS = LS[chn.J];  
    }

    return s_s + s_LS + s_j;
}

std::string quantum_channel_to_string_2(qs::quantum_channel chn)
{
    std::stringstream s; 
    s << "J=" << chn.J << "_S=" << chn.S 
        << "_T=" << chn.T << "_Tz=" << chn.Tz << 
        "_coup=" << chn.coupled;
    return s.str();
}

bool chn_in_chns(qs::quantum_channel chn, std::vector<qs::quantum_channel> chns)
{
    for (int i = 0; i < chns.size(); i++)
    {
        if (chns[i].J       == chn.J &&
            chns[i].S       == chn.S &&
            chns[i].T       == chn.T &&
            chns[i].Tz      == chn.Tz &&
            chns[i].coupled == chn.coupled) 
        {
            return true;
        }
    }
    return false;
}
