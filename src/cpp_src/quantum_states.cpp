#include "quantum_states.h"


std::vector<qs::quantum_states> get_states_NN(unsigned int J_max, unsigned int J_min, int Tz_min,
    int Tz_max, bool print)
{
    std::vector<qs::quantum_states> states;

    for (int Tz = Tz_min; Tz < Tzmax+1; Tz ++)
    {
        for (int J = J_min; J < J_max +1, J++)
        {
            for (int S = 0; S < 2; S++)
            {
                for (int L = std::abs(J-S); L < J+S+1; L++)
                {
                    for (int T = std::abs(Tz); T < 2; T++)
                    {
                        if (L+S+T %2 != 0) // Pauli principle
                        {
                            qs::quantum_state qs = 
                            {.J = J,
                            .L=L,
                            .S=S,
                            .T=T,
                            .Tz=Tz,
                            .pi = (int)std::pow(-1,L);
                            }
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

std::vector<qs::quntum_channel> get_channels(std::vector<qs::quantum_states> states)
{
    // Conserved quantum numbers are J,S,T,Tz,pi. 
    // For now these are hard coded in
    std::vector<qs::quantum_channel>
    for (std::size_t i = 0; i < states.size(); i++)
    {
        qs::quantum_NN_state bra = states[i];
        for (std::size_t j = 0; j < states.size(); j++)
        {
            qs::quantum_NN_state ket = states[j];

            if (bra.J == ket.J && bra.S == ket.S && bra.T == ket.T && bra.Tz == ket.Tz
                bra.pi = ket.pi)
            {

            }
        }
    }


}