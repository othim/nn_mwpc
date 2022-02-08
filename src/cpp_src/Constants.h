#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1


namespace constants
{
    // Constants (January 2022)
    // ----------------------------------
    const double gA  = 1.29; // Same as Andreas
    const double fpi = 92.4; // Same as Andreas
    const double mpi = 138.039; // Average of +,-,0 
    const double Mp  = 938.2720880259; // NIST
    const double Mn  = 939.5654203856; // NIST
    const double inv_fm_to_MeV = 197.3269804; // NIST
    // ----------------------------------

    // Nijmegen constants
    //const double Mp  = 938.27231;
    //const double Mn  = 939.56563;
    //const double mpi = 139.5675;
    //const double inv_fm_to_MeV = 197.326971941683;
    
    const double MeVm2_to_mbarn = inv_fm_to_MeV*inv_fm_to_MeV*10.0;
};

#endif
