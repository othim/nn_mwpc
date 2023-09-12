#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1


namespace constants
{
    #if !defined(ANDREAS_CONST) && !defined(NIJM_CONST)
    // Constants (January 2022)
    // ----------------------------------
    /*
    const double gA  = 1.275; // PDG (2022-03
    const double fpi = 92.1; // PDG, ch 71 (2022-03)
    const double mpi = 138.039; // Average of +,-,0 (PDG 2022-03) 
    const double Mp  = 938.2720880259; // NIST
    const double Mn  = 939.5654203856; // NIST
    const double inv_fm_to_MeV = 197.3269804; // NIST
    */
    // ----------------------------------
    
    // for benchmark
    const double gA  = 1.29;
    const double fpi = 92.4;
    const double mpi = 138.0;
    const double Mp = 939;
    const double Mn = 939;
    const double inv_fm_to_MeV = 197.3269804; // NIST
    #endif

    #ifdef ANDREAS_CONST
    // These constants are for testing the phase shifts
    const double gA  = 1.29;
    const double fpi = 92.4;
    const double mpi = 138.039;
    const double Mp  = 938.2720;
    const double Mn  = 939.5653;
    const double inv_fm_to_MeV = 197.3269804; // NIST
    #endif
    #ifdef NIJM_CONST
    // Nijmegen constants
    const double gA  = 1.29;
    const double fpi = 92.4;
    const double Mp  = 938.27231;
    const double Mn  = 939.56563;
    const double mpi = 139.5675;
    const double inv_fm_to_MeV = 197.326971941683;
    #endif

    const double MeVm2_to_mbarn = inv_fm_to_MeV*inv_fm_to_MeV*10.0;
};

#endif
