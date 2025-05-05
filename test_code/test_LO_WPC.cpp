
/*
 *
 *
 *
 * This file contains code that computes phase shifts and the PB observable 
 * in WPC.
 *
 *
 * Oliver Thim 2024
 *
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <vector>
#include "pybind_interface.h"
#include "pybind_interface_dwb.h"


/*
 * Function declarations
 */
void check_observable(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string, bool print);

void check_chn( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string, bool print);

void check_binding( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, int chn_number, bool print);

void new_log();
void print_to_log(std::string data);

void print_color(std::string data, std::string color);

struct arg_struct
{
    std::vector<std::string> args;
    std::vector<std::string> arg_value;
};
arg_struct parse_arguments(int argc, char** argv);


/*
 * Constants
 */

double TOL = 0.03;
bool print = true;
// Location of test data
std::string DATA_DIR = "../data/"; 

// Where test logfiles are written
std::string OUT_DIR  = "./out/";
std::string LOG_DIR  = "./test_logfiles/";

int main(int argc, char** argv)
{
    arg_struct ar = parse_arguments(argc,argv);
    
    // Set log dir
    LOG_DIR += ar.arg_value[1];

    // Create new log-file
    new_log();

    bool print = false;
    if (ar.arg_value[2] == "true")
    {
        print = true;
    }
    /*
     * ************************************************************************
     *                             Constants
     * ************************************************************************
     */

    unsigned int number_of_p_points = 100;   // Number of momentum-grid points
    int          J_max              = 20;
 
    double       Lambda             = 500.0;
    bool         finite_grid        = false;
    int          cut_pow            = 6;
    bool         sharp_cutoff       = false;
    bool         pre_comp_pot       = true;
    bool         rel_corr           = true;
    bool         inc_weights_in_pot = false;
    bool         cut_on_shell       = true;


    // Define the program constants (All in MeV)
    // ****************************
    ph::constants_struct* program_const = new ph::constants_struct;


    program_const->fpi            = 92.4;
    program_const->mpi            = 138.039;        // Average of +,-,0 (PDG 2022-03) 
    program_const->Mp             = 938.2720; 
    program_const->Mn             = 939.5653; 
    program_const->inv_fm_to_MeV  = 197.326971941683;    
    program_const->MeVm2_to_mbarn = program_const->inv_fm_to_MeV
        *program_const->inv_fm_to_MeV*10.0;

    // ****************************

    // Do pre-computations for the physics helper functions
    ph::physics_helpers_init();
    

    check_binding(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot, 
            rel_corr, number_of_p_points, finite_grid, inc_weights_in_pot, 
            cut_on_shell, program_const, 3, print);
    // Check angular observables
    check_observable(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot, 
            rel_corr, number_of_p_points, finite_grid, inc_weights_in_pot, 
            cut_on_shell, program_const, "PB", print);
    

    std::vector<std::string> chn_strings = {"1S0","1P1","3P0","3P1","3S1",
        "3D1","E1","3P2","3F2","E2"};
    
    for (auto& s : chn_strings) 
    {
        check_chn(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot, 
                rel_corr, number_of_p_points, finite_grid, inc_weights_in_pot, 
                cut_on_shell, program_const, s, false);
    }
    // Free allocated memory
    ph::physics_helpers_free();
    delete program_const;
    return 0;
}

void check_observable(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string, bool print)
{
    std::cout << "Testing observables with the WPC_LO potential" << std::endl;
    std::cout << "------------------------------------------------\n\n"; 

    nn_mwpc_interface interface = nn_mwpc_interface("WPC_LO",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,
         program_const->fpi,program_const->mpi,program_const->Mp,
         program_const->Mn,program_const->inv_fm_to_MeV);   

    double Tl = 30.0;
    double C1S0	= -0.1/100.0; 
    double C3S1	= -0.13/100.0;
    //double C3S1	= -0.072/100.0;
    
    // Open file
    std::string data = "../data/PB_30_MeV_Andreas_original.txt";
    std::ifstream infile(data);
    if (infile.is_open())
    {
        std::cout << "File" + data + " loaded: OK" << std::endl;
    } else
    {
        std::cout << "File" + data + ": Failed" << std::endl;
    }
    // Read the first line, the heading
    std::string heading;
    std::getline(infile, heading);
    //std::cout << heading << std::endl;
    
    double theta_obs[180];
    double data_obs[180];
    double errors[180];

    int k=0;
    double t_cm, t_lab,PB, th;
    while(infile >> t_cm >> t_lab >> th >> PB >> th >> th >> th >> th) 
    {
        data_obs[k]   = PB;
        theta_obs[k]  = t_cm;
        //std::cout << t_cm << "   " << PB << std::endl;
        k++;
    }
    
    //
    // Open outfile
    //
    std::ofstream myfile;
    std::string filename;

    std::string path = OUT_DIR + "out_" + obs_string + "_" + std::to_string((int)Tl) + ".txt"; 
    myfile.open(path);
        
    std::string obs_string2;
    if (obs_string == "DSG") {
        obs_string2 = "I 0000";
    } else if (obs_string == "PB") {
        obs_string2 = "P n000";
    } else if (obs_string == "DT") {
        obs_string2 = "K 0nn0";
    } else if (obs_string == "PT") {
        obs_string2 = "P 0n00";
    } else if (obs_string == "CKK") {
        obs_string2 = "A 00mm";
    } else if (obs_string == "AYY") {
        obs_string2 = "C nn00";
    } else if (obs_string == "AZZ") {
        obs_string2 = "A 00kk";
    } else {
        obs_string2 = obs_string; // SGT, SGTL, SGTT
    }
    std::vector<double> LECs = {C1S0,C3S1};
    interface.solve_LS(Tl,LECs);
    std::cout << "Angle \t obs \t correct \t abs. rel. error" << std::endl;    
    for (int i=0; i<180; i++)
    {
        double theta_deg = theta_obs[i];
        double obs       = interface.compute_observable(obs_string2,theta_deg);
        //std::cout << theta_deg << std::endl;
        if (data_obs[i] != 0) {
            errors[i] = std::abs((data_obs[i] - obs)/data_obs[i]);
        } else {
            errors[i] = 0;
        }
        myfile << theta_deg << "\t" << obs << "\t" << data_obs[i] << "\t" << 
            errors[i] << std::endl;
        if (print)
        {
            std::cout << theta_deg << "\t" << obs << "\t" << data_obs[i] << "\t" << 
                errors[i] << std::endl;
        }
    }
    myfile.close();
    double mean_error = 0;
    for (int i=0; i<180; i++)
    {  
        mean_error += errors[i];
    }
    mean_error /= 180.0;
    std::cout << "Mean absolute relative error: " << mean_error/180.0  
        << std::endl;
    std::cout << "Maximum error: " << *(std::max_element(errors, errors + 180)) 
        << std::endl;
}

void check_chn( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string chn_string, bool print)
{
    print = false;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Testing  phase shifts with the WPC_LO potential." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Channel: " << chn_string << std::endl;
    print_to_log("Channel: " + chn_string + "\n");
    std::vector<std::string> files;
    
    // Make constants the same
    // -----------------------
    // Set the constants to the values Andreas use
    
    double C1S0	= -0.1/100.0; 
    //double C1S0	= 0;//-0.1/100.0; 
    double C3S1	= -0.13/100.0;
    
    nn_mwpc_interface interface = nn_mwpc_interface("WPC_LO",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,
         program_const->fpi,program_const->mpi,program_const->Mp,
         program_const->Mn,program_const->inv_fm_to_MeV);   
    // -----------------------
    
    // Open file
    std::string data = "../data/phase_shifts_Andreas_original_part.txt";
    std::ifstream infile(data);
    if (infile.is_open())
    {
        std::cout << "File" + data + " loaded: OK" << std::endl;
    } else
    {
        std::cout << "File" + data + ": Failed" << std::endl;
    }

    // Read and save the data to arrays
    double D_E[350];
    double D_1S0[350];
    double D_3S1[350];
    double D_3D1[350];
    double D_E1[350];
    double D_1P1[350];
    double D_3P1[350];
    double D_3P0[350];
    double D_3P2[350];
    double D_3F2[350];
    double D_E2[350];
    
    const int nrow = 350;
    
    double* computed_tmp = new double[nrow];

    double E, a1S0, a3S1, aE1, a1P1, a3P1, a3P0, a3P2, aE2, a3F2, a3D1;
    int k = 0;
    //# Tlab,np 1S0,np 3S1, np E1, np 1P1, np 3P1, np 3P0, np 3P2, np E2,    
    //  np 3D1, np 3F2 

    // Read the first line, the heading
    std::string heading;
    std::getline(infile, heading);
    //std::cout << heading << std::endl;
    while(infile >> E >> a1S0 >> a3S1 >> aE1 >> a1P1 >> a3P1 >>
            a3P0 >> a3P2 >> aE2 >> a3D1 >> a3F2) 
    {
        D_E[k]   = E;
        D_1S0[k] = a1S0;
        
        D_3S1[k] = a3S1;
        D_3D1[k] = a3D1;
        D_E1[k]  = aE1;
        
        D_1P1[k] = a1P1;
        D_3P1[k] = a3P1;
        D_3P0[k] = a3P0;
        
        D_3P2[k] = a3P2;
        D_3F2[k] = a3F2;
        D_E2[k]  = aE2;

        /*std::cout << E << "   " << a1S0 << "   " << a3S1 << "   " << 
            aE1 << "   " << a1P1 << "   " << a3P1 << "   " << 
            a3P0 << "   " << a3P2 << "   " << "   " << aE2 << "   " 
            << a3D1 << "   " << a3F2 << std::endl;
        */
        k++;
        //std::cout << E << "   " << phase << std::endl;
    }
    
   
    
    std::string chn_string_full = chn_string;
    if (chn_string=="3S1" || chn_string=="3D1" ||chn_string=="E1")
    {
        chn_string_full = "3S-D1";
    } else if (chn_string=="3P2" || chn_string=="3F2" || chn_string=="E2")
    {
        chn_string_full = "3P-F2";
    }
    // Take just the relevant channel
    int chn_number = 0;
    for (int i = 0; i < (int)interface.get_chn_len(); i++)
    {
        if (interface.get_chn_LS_term(i) == chn_string_full) {
            chn_number = i;
            break;
        }
    }
    
    double max_err  = 0;
    double mean_err = 0;
    if (print)
    {
        std::cout << std::setw(5) << "T_lab" << "   " << std::setw(10) << 
            "data" << "   " << std::setw(10) << "C_phase" 
            << "   " << std::setw(8) << "err" << std::endl;
    }
    for (int i = 0; i < 350; i++)
    {
        double Tl = (double)(i+1);
        std::vector<double> LECs = {C1S0,C3S1};
        std::vector<double> phases_vec = interface.compute_phase_shift(chn_number,Tl,LECs);
        Phase_shifts_chn  phases;
        phases.delta_p = phases_vec[0];
        phases.delta_m = phases_vec[1];
        phases.epsilon = phases_vec[2];
        phases.delta_uncoupled = phases_vec[3];


        double* D_phase = nullptr;
        double C_phase  = 0.0;
        if (chn_string == "1S0")
        {
            D_phase = D_1S0;
            C_phase = phases.delta_uncoupled*180.0/M_PI;
            //std::cout << Tl << "\t" << C_phase << std::endl;
            
        } else if (chn_string == "1P1")
        {
            D_phase = D_1P1;
            C_phase = phases.delta_uncoupled*180.0/M_PI;
        } else if (chn_string == "3P0")
        {
            D_phase = D_3P0;
            C_phase = phases.delta_uncoupled*180.0/M_PI;
        } else if (chn_string == "3P1")
        {
            D_phase = D_3P1;
            C_phase = phases.delta_uncoupled*180.0/M_PI;
        } else if (chn_string == "3S1")
        {
            D_phase = D_3S1;
            C_phase = phases.delta_m*180.0/M_PI;
        } else if (chn_string == "3D1")
        {
            D_phase = D_3D1;
            C_phase = phases.delta_p*180.0/M_PI;
        } else if (chn_string == "E1")
        {
            D_phase = D_E1;
            C_phase = phases.epsilon*180.0/M_PI;
        } else if (chn_string == "3P2")
        {
            D_phase = D_3P2;
            C_phase = phases.delta_m*180.0/M_PI;
        } else if (chn_string == "3F2")
        {
            D_phase = D_3F2;
            C_phase = phases.delta_p*180.0/M_PI;
        } else if (chn_string == "E2")
        {
            D_phase = D_E2;
            C_phase = phases.epsilon*180.0/M_PI;
        }
        computed_tmp[i] = C_phase;
        double err = std::abs(D_phase[i] - C_phase);
        if (print)
        {
        std::cout << std::setw(5) << D_E[i] << "   " << std::setprecision(8) 
            << std::setw(10) << 
            D_phase[i] << "   " << std::setw(10) << C_phase 
            << "   " << std::setw(8) << err << std::endl;
        
        }
        max_err = std::max(max_err, err);
        mean_err += err; 

    }
    
    std::cout << "Mean error (deg): " << mean_err/350.0 << std::endl <<
        "Max error (deg): " << max_err << std::endl;

    double tol = 2e-4;
    if (mean_err/350.0<tol && max_err<tol)
    {
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "******** Test: OK (abs.tol=" << tol << ") ********" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        print_to_log("TEST: SUCCESS\n\n");
        std::cout << std::endl << std::endl << std::endl;
    } else 
    {   
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "******** Test: FAILED (abs.tol=" << tol << ") ********" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << std::endl << std::endl << std::endl;
        print_to_log("TEST: FAILED\n\n");

    }
}

void check_binding( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, int chn_number, bool print)
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Testing binding energy with the WPC_LO potential." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    print_to_log("Testing binding energy 3S-D1\n");
    
    // Make constants the same
    // -----------------------
    // Set the constants to the values Andreas use
    
    double C1S0	= -0.1/100.0; 
    double C3S1	= -0.13/100.0;
    //double C3S1	= -0.072/100.0;
    
    nn_mwpc_interface interface = nn_mwpc_interface("WPC_LO",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,
         program_const->fpi,program_const->mpi,program_const->Mp,
         program_const->Mn,program_const->inv_fm_to_MeV);   
    // -----------------------
        
    std::vector<double> LECs = {C1S0,C3S1};
    std::vector<double> eigenvalues = interface.compute_binding_energy(chn_number,LECs);
    std::cout << "Eigenvalues in 3S-D1 LO WPC:" << std::endl;
    for (auto e : eigenvalues) 
    {
        if (e<0.0) {
            std::cout << e << std::endl;
            print_to_log("Computed: " + std::to_string(e)+", old: (-19.8288)\n\n");
        }
    }
    std::cout << "(-19.8288)" << std::endl;

}

arg_struct parse_arguments(int argc, char** argv)
{

    arg_struct ar;
    ar.args = {"-t","-f","-p","-h"};
    ar.arg_value = {"phase","test_LO_WPC.out","false","false"};

    bool help = false;
    for (int i=0; i<ar.args.size(); i++)
    {
        for (int j=0; j<argc; j++)
        {
            std::string s(argv[j]);
            if (s == ar.args[i])
            {
                if (s == "-h")
                {
                    help = true;
                } else {
                    ar.arg_value[i] = argv[j+1];
                }
            }
        }
    }
    
    if (help)
    {
        std::cout << "Arguments:" << std::endl;
        for (int i=0; i<ar.args.size(); i++)
        {
            std::cout << ar.args[i] << " : " << ar.arg_value[i] << std::endl;
        }
        std::exit(0);
    }
    std::cout << std::endl;
    return ar;
}

void new_log()
{
    std::ofstream myfile;
    myfile.open(LOG_DIR);
    myfile << "New log:" << LOG_DIR << std::endl;
    std::time_t result = std::time(nullptr);
    myfile << std::asctime(std::localtime(&result)) << std::endl;
    myfile.close();
}

void print_to_log(std::string data)
{
    std::ofstream myfile;
    myfile.open(LOG_DIR, std::ios_base::app);
    myfile << data;
    myfile.close();
}

void print_color(std::string data, std::string color)
{
    if (color == "green")
    {
        std::cout << "\033[0;32m"+data+"\033[0m\n";
    } else if (color == "red")
    {
        std::cout << "\033[0;31m"+data+"\033[0m\n";
    }
}
