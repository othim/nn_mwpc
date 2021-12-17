#include "scattering.h"


/*
    This function computes the on-shell T-matrix elements from a vector of
    phase shifts. Phase shifts are in the Stapp convention.
*/
std::vector<std::complex<double>*> T_from_phase_shifts(std::vector<Phase_shifts_chn> phases, std::vector<qs::quantum_channel> chns_vec,double rho_T)
{
    // T = i/(2\rho_T) * (S-1)
    // This is from inverting the on-shell relation between T and S
    const std::complex<double> imag_u(0.0,1.0);
    std::complex<double> fac = 1.0;//imag_u/(2.0*rho_T);

    //std::cout << fac << std::endl;
    
    std::vector<std::complex<double>*> T_vec;

    for (std::size_t i = 0; i < phases.size(); i++)
    {
        std::complex<double>* T = new std::complex<double>[3]; // TODO needs to be deallocated
        T[0] = (std::complex<double>)0.0;
        T[1] = (std::complex<double>)0.0;
        T[2] = (std::complex<double>)0.0;

        if (chns_vec[i].coupled == false)
        {
            T[0] = fac*(std::exp(imag_u*2.0*phases[i].delta_uncoupled) - 1.0);
        } else
        {
            double two_eps = 2.0*phases[i].epsilon;
            double dm = phases[i].delta_m;
            double dp = phases[i].delta_p;
            
            T[0] = (std::complex<double>)(fac*(std::cos(two_eps)*std::exp(imag_u * 2.0* dm ) - 1.0)); // Tmm
            T[1] = (std::complex<double>)(fac*(imag_u*std::sin(two_eps)*std::exp(imag_u*(dm+dp)))); // Tmp
            T[2] = (std::complex<double>)(fac*(std::cos(two_eps)*std::exp(imag_u * 2.0* dp ) - 1.0)); // Tpp
        
        }
        T_vec.push_back(&T[0]);
    }
    
    
    /*std::cout << "-------" << std::endl;
    for (int i = 0; i < T_vec.size(); i++)
    {
       std::cout << T_vec[i][0] << " " << T_vec[i][1] << " " <<  T_vec[i][2] << std::endl;
    }
    std::cout << "-------" << std::endl;*/
    return T_vec;
}



std::complex<double> get_M_matrix_p(std::vector<qs::quantum_channel> chns_vec,
    std::vector<Phase_shifts_chn> phase_shifts_vec, int s, int mo, int mi, 
        double cos_theta, double q_on_shell,double rho_T,int l_max)
{
    #ifdef DEGUB_SCATTERING
        std::cout << "get_M_matrix_p()" << std::endl;
    #endif

    // Convert phase shifts back to T-matrix
    // Just needs to done once...
    std::vector<std::complex<double>*> T_vec = T_from_phase_shifts(phase_shifts_vec,chns_vec,rho_T);
    
    // Call get_M_Matrix_T
    //return std::complex<double>(0,0);
    return sc::get_M_matrix_T(chns_vec,T_vec,q_on_shell,s,mo,mi,cos_theta,l_max);
}

/*
    Computes the (lo,li) pairs that are possible from
    the quantum numbers J and S.
*/
std::vector<lo_li> get_ls(int J, int S, bool coupled)
{
    //std::cout << J << " " << S << " " << coupled << std::endl;
    std::vector<lo_li> elements;
    if (S == 1 && J != 0)
    {
        if (coupled)
        {
            elements.push_back((lo_li) {.lo = J-1, .li = J-1});
            elements.push_back((lo_li) {.lo = J-1, .li = J+1});
            elements.push_back((lo_li) {.lo = J+1, .li = J-1});
            elements.push_back((lo_li) {.lo = J+1, .li = J+1});
        } else 
        {
            elements.push_back((lo_li) {.lo = J, .li = J});
        }
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
    std::vector<std::complex<double>*> T_on_shell_vec, double q_on_shell,int s, int mo, int mi, double cos_theta,int l_max)
{ 
    const std::complex<double> imag_u(0.0,1.0);  
    // Move to some constructor-ish
    wig_table_init(2*100, 9);
    wig_temp_init(2*100);

    // Precompute some spherical harmonics for the given cos_theta
    double* sph_arr = (double*) malloc(gsl_sf_legendre_array_n(l_max)*sizeof(double)); // Not including coplex e^{im\phi} phase. But we set \phi = 0 anyway!
   /* 
    if (cos_theta < -1 || cos_theta > 1)
    {
        std::cout << "Error" << std::endl;
        std::cout << cos_theta << std::endl;
        double a;
        std::cin >> a;
    }*/
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,l_max,cos_theta,1,sph_arr);

    std::complex<double> result = 0;
    
    
    // Sum over channels, j-sum
    for (std::size_t i = 0; i < chns_vec.size(); i++) //TODO
    {
        //std::cout << "i=" << i << std::endl;
        qs::quantum_channel current_chn = chns_vec[i];
        int J = current_chn.J;
        bool coupled = current_chn.coupled;
        // Sum over the correct s-values
        if (current_chn.S == s)
        {
            // Compute (lo,li) pairs for these quantum numbers
            std::vector<lo_li> Ls = get_ls(J,s,coupled); 
            //std::cout << "Ls-len: " << Ls.size() << std::endl;
            // Loop over allowed pairs (lo,li)
            for (std::size_t j = 0; j < Ls.size(); j++)
            {   
                int li = Ls[j].li;
                int lo = Ls[j].lo;
                //std::cout << "li:" << li << " lo:" << lo << std::endl; 

                // Check if mi and mo are compatible with channel
                if (!(abs(mi-mo) > lo || abs(mi) > s || abs(mo) > s || abs(mi) > J))
                {
                    // Compute Y_lm and wiegner 3j symbols
                    /*if (std::abs(mi - mo)> 1)
                    {
                        std::cout << "Error" << std::endl;
                        
                        std::cout << s << " " << mo << " " << mi << std::endl;
                        double a;
                        std::cin >> a;
                    }*/

                    // This is necessary since gsl only calculates the 
                    // elements of positive m
                    double y_lm;
                    if (-(mi-mo)<0) {
                        //y_lm = gsl_sf_legendre_sphPlm(lo,(mi-mo), cos_theta);
                        y_lm = std::pow(-1,(mi-mo))*sph_arr[gsl_sf_legendre_array_index(lo, (mi-mo))];
                    } else {
                        //y_lm = std::pow(-1, -(mi-mo))*gsl_sf_legendre_sphPlm(lo,-(mi-mo),cos_theta);
                        
                        y_lm = sph_arr[gsl_sf_legendre_array_index(lo, -(mi-mo))];
                    }

                    double wig1 = wig3jj(2*  lo , 2*  s , 2*  J ,
                                         2*  (mi-mo) , 2*  mo , -2*  mi );

                    double wig2 = wig3jj(2*  li , 2*  s , 2*  J ,
                                         2*  0 , 2*  mi , 2*  -mi );
                    
                    // Take correct T-matrix element 
                    std::complex<double> T_el;

                    if ((Ls.size() == 1))
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
                    /*std::cout << "T_el: " << T_el << std::endl;
                    std::cout << "y_lm: " << y_lm << std::endl;
                    std::cout << "wig1: " <<  wig1 << std::endl;
                    std::cout << "wig2: " << wig2 << std::endl;*/
                    
                    // OLD
                    std::complex<double> add = std::pow(imag_u,(li-lo)) * (std::complex<double>) (2.0*J+1)*sqrt(2*li+1)*y_lm*wig1*wig2*T_el;
                    // ---
                    
                    //std::complex<double> add = std::pow(imag_u,(li-lo)) * (std::complex<double>) (2.0*J+1)*sqrt(1/(4*M_PI))*y_lm*wig1*wig2*T_el;
                    
                    result += add;
                } // end if
            } // end loop over (lo,li)
        } // end S=s
    } // end for chn
    free(sph_arr);
    wig_temp_free();
    wig_table_free();


    // Add factor in front 
    result *= (imag_u) * (std::complex<double>) (-sqrt(M_PI)/q_on_shell);
    // Return


    return result;
}


std::vector<std::complex<double> > compute_Saclay_amplitudes(std::vector<qs::quantum_channel> chns_vec,
    std::vector<Phase_shifts_chn> phase_shifts_vec, double theta, double q_on_shell, double rho_T, int l_max)
{
    // Compute M-matrix elements
    std::complex<double> M_pp = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)1,(int)1,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_00 = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)0,(int)0,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_pm = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)1,(int)-1,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_s =
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)0,(int)0,(int)0,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_p0 = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)1,(int)0,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_0p = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)0,(int)1,std::cos(theta),q_on_shell,rho_T,l_max);
    
    // Compute Saclay amplitudes as in eq 2.14 in 
    // Formalism of nucleon-nucleon elastic scattering experiments. 
    // Journal de Physique, 1978, 39 (1), pp.1-32.

    std::complex<double> a = (1.0/2.0) * (M_pp + M_00 - M_pm);
    std::complex<double> b = (1.0/2.0) * (M_pp + M_s + M_pm);
    std::complex<double> c = (1.0/2.0) * (M_pp - M_s + M_pm);

    std::complex<double> d = (1.0/(2*std::cos(theta)))*(-M_pp + M_00 + M_pm);
    //std::complex<double> d = -(1.0/(std::sqrt(2)*std::sin(theta)))*(M_p0 + M_0p);
    const std::complex<double> imag_u(0.0,1.0);  
    std::complex<double> e = (imag_u/std::sqrt(2))*(M_p0 - M_0p);  
    
    double fac = std::sqrt(constants::MeVm2_to_mbarn);
    
    //std::cout << "M_00: " << M_00*fac << " M_pm: " << M_pm*fac << " M_s: " << M_s*fac << std::endl;
    //std::cout << "M_pp: " << M_pp*fac << " M_10: " << M_p0*fac << " M_01:" << M_0p*fac << std::endl;
    
    std::vector<std::complex<double> > sac_amp;
    
    // This factor is to convert from MeV^{-1} to mbarn^{0.5} 
    sac_amp.push_back(a*fac);
    sac_amp.push_back(b*fac);
    sac_amp.push_back(c*fac);
    sac_amp.push_back(d*fac);
    sac_amp.push_back(e*fac);
    return sac_amp;
}


double compute_observable(std::vector<std::complex<double> > sac_amp,std::string obs)
{
    // The vector of sac_amp contains a,b,c,d,e in order
    std::complex<double> a = sac_amp[0];
    std::complex<double> b = sac_amp[1];
    std::complex<double> c = sac_amp[2];
    std::complex<double> d = sac_amp[3];
    std::complex<double> e = sac_amp[4];
    
    double DSG = (1.0/2.0)*(std::abs(a)*std::abs(a) + std::abs(b)*std::abs(b) + std::abs(c)*std::abs(c)
            + std::abs(d)*std::abs(d) + std::abs(e)*std::abs(e));
    
    if (obs == "I 0000" || obs == "C nnnn")
    {
        return DSG;
    } else if (obs == "C nn00" || obs == "A 00nn")
    {
        return (1.0/(2.0*DSG))*(std::abs(a)*std::abs(a) - std::abs(b)*std::abs(b) - std::abs(c)*std::abs(c)
            + std::abs(d)*std::abs(d) + std::abs(e)*std::abs(e));
    } else if (obs == "D n0n0" || obs == "D 0n0n")
    {
        return (1.0/(2.0*DSG))*(std::abs(a)*std::abs(a) + std::abs(b)*std::abs(b) - std::abs(c)*std::abs(c)
            - std::abs(d)*std::abs(d) + std::abs(e)*std::abs(e));
    } else if (obs == "K 0nn0" || obs == "K n00n")
    {
        return (1.0/(2.0*DSG))*(std::abs(a)*std::abs(a) - std::abs(b)*std::abs(b) + std::abs(c)*std::abs(c)
            - std::abs(d)*std::abs(d) + std::abs(e)*std::abs(e));
    } else if (obs == "C llll" || obs == "C mmmm")
    {
        return (1.0/(2.0*DSG))*(std::abs(a)*std::abs(a) + std::abs(b)*std::abs(b) + std::abs(c)*std::abs(c)
            + std::abs(d)*std::abs(d) - std::abs(e)*std::abs(e));
    } else if (obs == "P n000" || obs == "P 0n00" || obs == "A 00n0" || obs == "A 000n" 
        || obs == "C nnn0" || obs == "C nn0n" || obs == "M n0nn" || obs == "N 0nnn")
    {
        return std::real(std::conj(a)*e)/DSG;
    } else if (obs == "C lllm" || obs == "C llml" || obs == "C lmmm" || obs == "C mlmm")
    {
        return std::imag(std::conj(a)*e)/DSG;
    } else if (obs == "C lmll" || obs == "C mlll" || obs == "C mmlm" || obs == "C mmml")
    {
        return -std::imag(std::conj(a)*e)/DSG;
    } else if (obs == "C lnl0" || obs == "C mnm0" || obs == "C nl0l" || obs == "C nm0m" ||
        obs == "M m0mn" || obs == "N 0lnl" || obs == "N 0mnm")
    {
        return std::real(std::conj(b)*e)/DSG;
    } else if (obs == "D l0m0" || obs == "D 0l0m" || obs == "C nlnm" || obs == "C lnmn")
    {
        return std::imag(std::conj(b)*e)/DSG;
    } else if (obs == "D m0l0" || obs == "D 0m0l" || obs == "C mnln" || obs == "C nmnl")
    {
        return -std::imag(std::conj(b)*e)/DSG;
    } else if (obs == "C nll0" || obs == "C nmm0" || obs == "C ln0l" || obs == "C mn0m" 
        || obs =="M l0nl" || obs == "M m0nm" || obs == "N 0lln" || obs =="N 0mmn")
    {
        return std::real(std::conj(c)*e)/DSG;
    } else if (obs == "K 0lm0" || obs == "K l00m" || obs == "C nlmn" || obs == "C lnnm")
    {
            return std::imag(std::conj(c)*e)/DSG;
    } else if (obs == "K m00l" || obs == "K 0ml0" || obs == "C mnnl" || obs == "C nmln")
    {
        return -std::imag(std::conj(c)*e)/DSG;
    } else if (obs == "C lln0" || obs == "C ll0n" || obs =="M n0ll" || obs == "N 0nll")
    {
            return -std::real(std::conj(d)*e)/DSG; // OBS
    } else if (obs == "C mmn0" || obs == "C mm0n" || obs == "M n0mm" || obs == "N 0nmm")
    {
            return std::real(std::conj(d)*e)/DSG; // OBS
    } else if (obs == "C lm00" || obs == "C ml00" || obs == "C mlnn" || obs == "C lmnn")
    {
        return std::imag(std::conj(d)*e)/DSG;
    } else if (obs == "A 00lm" || obs == "A 00ml" || obs == "C nnml" || obs == "C nnlm")
    {
            return -std::imag(std::conj(d)*e)/DSG;
    } else if (obs == "D m0m0" || obs == "D 0m0m" || obs == "C nlnl" || obs == "C lnln")
    {
        return std::real(std::conj(a)*b + std::conj(c)*d)/DSG;
    } else if (obs =="C mnl0" || obs == "C nm0l")
    {
        return std::imag(std::conj(a)*b + std::conj(c)*d)/DSG;
    } else if (obs == "M l0mn" || obs == "N 0lnm")
    {
        return -std::imag(std::conj(a)*b + std::conj(c)*d)/DSG;
    } else if (obs == "D l0l0" || obs == "D 0l0l" || obs == "C nmnm" || obs == "C mnmn")
    {
        return std::real(std::conj(a)*b - std::conj(c)*d)/DSG;
    } else if (obs == "C lnm0" || obs == "C nl0m")
    {
        return -std::imag(std::conj(a)*b - std::conj(c)*d)/DSG;
    } else if (obs == "M m0ln" || obs == "N 0mnl")
    {
        return std::imag(std::conj(a)*b - std::conj(c)*d)/DSG;
    } else if (obs == "K 0mm0" || obs == "K m00m" || obs == "C nlln" || obs == "C lnnl")
    {
        return std::real(std::conj(a)*c + std::conj(b)*d)/DSG;

    } else if (obs == "C nlm0" || obs == "C ln0m")
    {
        return -std::imag(std::conj(a)*c - std::conj(b)*d)/DSG;
    } else if (obs == "M m0nl" || obs == "N 0mln")
    {
        return (std::imag(std::conj(a)*c - std::conj(b)*d))/DSG;
    } else if (obs == "C mm00" || obs == "A 00mm")
    {
        //std::cout << "here" << std::endl;
        return (std::real(std::conj(a)*d + std::conj(b)*c))/DSG;
    } else if (obs == "C nnll" || obs == "C llnn")
    {
        return -std::real(std::conj(a)*d + std::conj(b)*c)/DSG;
    } else if (obs == "C lmn0" || obs == "C ml0n")
    {
        return -std::imag(std::conj(a)*d + std::conj(b)*c)/DSG;
    } else if (obs =="M n0lm" || obs == "N 0nml") 
    {
        return std::imag(std::conj(a)*d + std::conj(b)*c)/DSG;
    } else if (obs == "C ll00" || obs == "A 00ll")
    {
        return -std::real(std::conj(a)*d - std::conj(b)*c)/DSG;
    } else if (obs == "C mmnn" || obs == "C nnmm")
    {
        return std::real(std::conj(a)*d - std::conj(b)*c)/DSG;
    } else if (obs == "C mln0" || obs == "C lm0n")
    {
        return -std::imag(std::conj(a)*d - std::conj(b)*c)/DSG;
    } else if (obs == "M n0ml" || obs == "N 0nlm")
    {
        return std::imag(std::conj(a)*d - std::conj(b)*c)/DSG;
    } else {
        std::cout << "Unknown observable" << std::endl;
        return 0;
    }    
}

double compute_total_cross_section(std::vector<qs::quantum_channel> chns_vec, 
    std::vector<Phase_shifts_chn> phase_shifts_vec,double q_on_shell,double rho_T,int l_max)
{
    double theta = 0.0;
    // Compute M-matrix elements
    std::complex<double> M_pp = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)1,(int)1,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_00 = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)0,(int)0,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_pm = 
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)1,(int)1,(int)-1,std::cos(theta),q_on_shell,rho_T,l_max);
    std::complex<double> M_s =
        get_M_matrix_p(chns_vec,phase_shifts_vec,(int)0,(int)0,(int)0,std::cos(theta),q_on_shell,rho_T,l_max);
 
    std::complex<double> a = (1.0/2.0) * (M_pp + M_00 - M_pm);
    std::complex<double> b = (1.0/2.0) * (M_pp + M_s + M_pm);

    double sigma = ((2*M_PI)/q_on_shell)*std::imag(a+b);

    return sigma*constants::MeVm2_to_mbarn; // Convert to milli barn
}

gsl_matrix_complex* get_M_matrix(std::vector<qs::quantum_channel> chns_vec,
    std::vector<Phase_shifts_chn> phase_shifts_vec, double q_on_shell, double theta, double rho_T, int l_max)
{
    gsl_matrix_complex* M = gsl_matrix_complex_alloc(4,4);
    
    // The M-matrix in the basis |+->_1 x |+->_2 is as follows
    // 1:- - + +
    // 2:- + - +  2  1
    //   * * * *  -  -
    //   * * * *  +  -
    //   * * * *  -  +
    //   * * * *  +  +
    // These matrix elements M_(-,+) can be calculated from
    // the function get_M_matrix_T with the insetion of apropriate 
    // Clebsh-Gordan coefficients.
    
    // Compute the basis change matrix

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // Calculate indices
            int m1r = (int)(2*((int)i / (int)2 - 1)); // -1 or +1
            int m1c = (int)(2*((int)j / (int)2 - 1)); // - || -
                
            int m2r = (int)(2*((int)i % (int)2 - 1)); // - || -
            int m2c = (int)(2*((int)j % (int)2 - 1)); // - || -
                     
            // Calculate the matrix elements 
            gsl_complex el = gsl_complex_rect(0.0,0.0);

            for (int S = 0; S < 2; S++) // S =0,1
            {
                // Note that this is with a factor 2.
                int Mi = m1r + m2r; // Row is out states
                int Mo = m1c + m2c; // Column is in states
                
                if (Mi > S*2 || Mo > S*2)
                {
                    continue;
                }
                std::complex<double> el_tmp = get_M_matrix_p(chns_vec, phase_shifts_vec,
                        S, (int)((double)Mo/2.0), (int)((double)Mi/2.0), std::cos(theta), q_on_shell, rho_T, l_max);
                std::cout << "el" << el_tmp << std::endl;
                gsl_complex M_el = gsl_complex_rect(std::real(el_tmp),std::imag(el_tmp));
                gsl_complex tmp = gsl_complex_mul(gsl_complex_rect(
                        ph::CG_coeff(2*S, Mi, 1, 1, m1r,m2r)*
                        ph::CG_coeff(2*S, Mo, 1, 1, m1c,m2c),0.0),M_el);
                std::cout << "tmp" << GSL_REAL(tmp) << "," << GSL_IMAG(tmp) << std::endl;
                el = gsl_complex_add(el,tmp);
            }            
            std::cout << "mel" << GSL_REAL(el) << "," << GSL_IMAG(el) << std::endl;
            gsl_matrix_complex_set(M,i,j,el);     
        }
    }
    return M;
}

double get_observables(gsl_matrix_complex* sigma_i_1, gsl_matrix_complex* sigma_i_2,
        gsl_matrix_complex* sigma_o_1, gsl_matrix_complex* sigma_o_2,
        gsl_matrix_complex* M_matrix)
{
    // Set up tensor products of the sigma-matrices 
    gsl_matrix_complex* in_tensor_prod  = ph::kronecker_product(sigma_i_1, sigma_i_2);
    gsl_matrix_complex* out_tensor_prod = ph::kronecker_product(sigma_o_1, sigma_o_2);
    //std::cout << "1" << std::endl;
    // Multiply the matrices
    gsl_complex alpha = gsl_complex_rect(1.0,0.0);
    gsl_complex beta  = gsl_complex_rect(0.0,0.0);

    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M_matrix->size1, M_matrix->size2);
    // sigma_in M sigma_out M^dagger -> out_tensor_prod
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, out_tensor_prod,M_matrix, beta,tmp );
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   alpha, M_matrix,tmp, beta, out_tensor_prod);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   alpha, in_tensor_prod,out_tensor_prod, beta, tmp);
    

    // Take the trace
    gsl_complex t = ph::trace(tmp);
    std::cout << GSL_REAL(t) << " " << GSL_IMAG(t) << std::endl;
    gsl_matrix_complex_free(tmp);
    // Return 
    return GSL_REAL(t)/4.0; // The normalization factor
}
