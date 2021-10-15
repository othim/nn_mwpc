#include "scattering.h"


/*
    This function computes the on-shell T-matrix elements from a vector of
    phase shifts. Phase shifts are in the Stapp convention
*/
std::vector<std::complex<double>*> T_from_phase_shifts(std::vector<Phase_shifts_chn> phases, std::vector<qs::quantum_channel> chns_vec,double rho_T)
{
    // T = i/(2\rho_T) * (S-1)
    // This is from inverting the on-shell relation between T and S
    const std::complex<double> imag_u(0.0,1.0);
    std::complex<double> fac = imag_u/(2.0*rho_T);

    
    std::vector<std::complex<double>*> T_vec;

    for (std::size_t i = 0; i < phases.size(); i++)
    {
        std::complex<double> T[3] = { (std::complex<double>)0.0,(std::complex<double>)0.0,
            (std::complex<double>)0.0 };

        if (chns_vec[i].coupled == false)
        {
            T[0] = fac*(std::exp(imag_u*2.0*phases[i].delta_uncoupled) - 1.0);
        } else
        {
            double two_eps = 2.0*phases[i].epsilon;
            double dm = phases[i].delta_m;
            double dp = phases[i].delta_p;
            
            T[0] = fac*(std::cos(two_eps)*std::exp(imag_u * 2.0* dm ) - 1.0); // Tmm
            T[1] = fac*(imag_u*std::sin(two_eps)*std::exp(imag_u*(dm+dp))); // Tmp
            T[2] = fac*(std::cos(two_eps)*std::exp(imag_u * 2.0* dm ) - 1.0); // Tpp
        }
        T_vec.push_back(&T[0]);
    }
    return T_vec;
}



std::complex<double> get_M_matrix_p(std::vector<qs::quantum_channel> chns_vec,
    std::vector<Phase_shifts_chn> phase_shifts_vec, unsigned int s, int mo, int mi, 
        double cos_theta, double q_on_shell,double rho_T,unsigned int l_max)
{
    #ifdef DEGUB_SCATTERING
        std::cout << "get_M_matrix_p()" << std::endl;
    #endif

    // Convert phase shifts back to T-matrix
    // Just needs to done once...
    std::vector<std::complex<double>*> T_vec = T_from_phase_shifts(phase_shifts_vec,chns_vec,rho_T);
    
    // Call get_M_Matrix_T
    //return std::complex<double>(0,0);
    return get_M_matrix_T(chns_vec,T_vec,q_on_shell,s,mo,mi,cos_theta,l_max);
}

/*
    Computes the (lo,li) pairs that are possible from
    the quantum numbers J and S.
*/
std::vector<lo_li> get_ls(unsigned int J, unsigned int S)
{
    std::vector<lo_li> elements;
    if (S == 1 && J != 0)
    {
        elements.push_back((lo_li) {.lo = J-1, .li = J-1});
        elements.push_back((lo_li) {.lo = J-1, .li = J+1});
        elements.push_back((lo_li) {.lo = J+1, .li = J-1});
        elements.push_back((lo_li) {.lo = J+1, .li = J+1});
    } else if (S == 0) // Ok if J==0 here
    {
        elements.push_back((lo_li) {.lo = J, .li = J});
    } else if (S == 1 && J == 0)
    {
        elements.push_back((lo_li) {.lo = J+1, .li = J+1});
    } else 
    {
        std::cerr << "Error in get_ls()" << std::endl;
    }

    return elements;
}

std::complex<double> get_M_matrix_T(std::vector<qs::quantum_channel> chns_vec,
    std::vector<std::complex<double>*> T_on_shell_vec, double q_on_shell,unsigned int s, int mo, int mi, double cos_theta,unsigned int l_max)
{ 
    const std::complex<double> imag_u(0.0,1.0);  
    // Move to some constructor-ish
    wig_table_init(2*100, 9);
    wig_temp_init(2*100);

    // Precompute some spherical harmonics for the given cos_theta
    double* sph_arr = (double*) malloc(gsl_sf_legendre_array_n(l_max)*sizeof(double)); // Not including coplex e^{im\phi} phase. But we set \phi = 0 anyway!

    gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM,l_max,cos_theta,sph_arr);

    std::complex<double> result = 0;
    
    
    // Sum over channels, j-sum
    for (std::size_t i = 0; i < chns_vec.size(); i++) //TODO
    {
        //std::cout << "i=" << i << std::endl;
        qs::quantum_channel current_chn = chns_vec[i];
        unsigned int J = current_chn.J;
        // Sum over the correct s-values
        if (current_chn.S == s)
        {
            // Compute (lo,li) pairs for these quantum numbers
            std::vector<lo_li> Ls = get_ls(J,s); 
            //std::cout << "Ls-len: " << Ls.size() << std::endl;
            // Loop over allowed pairs (lo,li)
            for (std::size_t j = 0; j < Ls.size(); j++)
            {   
                unsigned int li = Ls[j].li;
                unsigned int lo = Ls[j].lo;
                std::cout << "li:" << li << " lo:" << lo << std::endl; 
                // Check if mi and mo are compatible with channel
                if (!(abs(mi-mo) > lo || abs(mi) > current_chn.S || abs(mi) > J))
                {
                    // Compute Y_lm and wiegner 3j symbols
                    double y_lm = sph_arr[gsl_sf_legendre_array_index(lo, mi-mo)];
                    double wig1 = wig3jj(2*  lo , 2*  s , 2*  J ,
                                         2*  mi-mo , 2*  mo , 2*  mi );

                    double wig2 = wig3jj(2*  li , 2*  s , 2*  J ,
                                         2*  0 , 2*  mi , 2*  -mi );
                    
                    // Take correct T-matrix element 
                    std::complex<double> T_el;

                    if (!(Ls.size() == 1))
                    {
                        T_el = T_on_shell_vec[i][0]; // Uncoupled case
                    } else
                    {
                        if (lo == li)
                        {
                            if (current_chn.J == li-1)
                            {
                                T_el = T_on_shell_vec[i][2]; // Tpp
                            } else
                            {
                                T_el = T_on_shell_vec[i][0]; // Tmm
                            }
                        } else
                        {
                            T_el = T_on_shell_vec[i][1]; // Tmp
                        }
                    }
                    result += std::pow(imag_u,(li-lo)) * (std::complex<double>) (2*J+1)*sqrt(2*li+1)*y_lm*wig1*wig2*T_el;
                    //std::cout << "Result" << result << std::endl;
                } // end if
            } // end loop over (lo,li)
        } // end S=s
    } // end for chn
    free(sph_arr);
    // Add factor in front 
    result *= (imag_u) * (std::complex<double>) (-sqrt(M_PI)/q_on_shell);
    // Return
    return result;
}

