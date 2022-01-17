#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1


namespace constants
{
    const double gA  = 1.29;//1.289;
    const double fpi = 92.4; // 92.2
    //const double mpi = 138.04;//138.039;
    //const double Mp  = 938.27208;
    //const double Mn  = 939.56542;
    //const double inv_fm_to_MeV = 197.3269789;

    // Nijmegen constants
    const double Mp  = 938.27231;
    const double Mn  = 939.56563;
    const double mpi = 139.5675;
    const double inv_fm_to_MeV = 197.326971941683;
    
    const double MeVm2_to_mbarn = inv_fm_to_MeV*inv_fm_to_MeV*10.0;
};

#endif
