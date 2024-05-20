/*
 *
 *
 *
 * This file contans code that test the implementation of phase shifts,
 * scattering observables and M-matrix elementes by using the nijmegen1
 * potential.
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
        unsigned int J_max_in_pot,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string, bool print);

void check_observable_noang(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string,bool print);

void check_M_el(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, bool print);

void check_binding( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, int chn_number, bool print);

void new_log();
void print_to_log(std::string data);

void print_color(std::string data, std::string color);

void load_data(std::string path,double* theta_obs,double* data_obs);
void load_data_M(std::string path,double* theta,double* data_real,
        double* data_imag);

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

    double       scale              = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points     = 76;    // Number of points in angular integration
    unsigned int number_of_p_points = 100;   // Number of momentum-grid points
    unsigned int J_max_in_pot       = 50;    // Maximum J that is stored for L-polynomials
    int          J_max              = 14;
    int          J_min              = 0;
    int          Tz_min             = 0;
    int          Tz_max             = 0;
 
    double       Lambda             = 5000.0;
    bool         finite_grid        = false;
    double       finite_grid_max    = 10000.0;
    int          cut_pow            = 6;
    bool         sharp_cutoff       = false;
    bool         pre_comp_pot       = false;
    bool         rel_corr           = false;
    bool         inc_weights_in_pot = false;
    bool         cut_on_shell       = true;
    double sharp_cutoff_add = 0;


    // Define the program constants (All in MeV)
    // ****************************
    ph::constants_struct* program_const = new ph::constants_struct;


    program_const->fpi            = 92.4;
    program_const->mpi            = 139.5675;        // Average of +,-,0 (PDG 2022-03) 
    program_const->Mp             = 938.27231; 
    program_const->Mn             = 939.56563; 
    program_const->inv_fm_to_MeV  = 197.326971941683;    
    program_const->MeVm2_to_mbarn = program_const->inv_fm_to_MeV
        *program_const->inv_fm_to_MeV*10.0;

    // ****************************

    // Do pre-computations for the physics helper functions
    ph::physics_helpers_init();
    

    check_binding(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot,
            rel_corr,number_of_p_points, finite_grid, inc_weights_in_pot, 
            cut_on_shell, program_const,3,print);
    // Check angular observables
    
    std::vector<std::string> obs_strings = {"DSG","PB","DT","AYY","AZZ","PT",
    "CKK"};

    for (auto obs_string : obs_strings)
    {
        check_observable(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot,
                rel_corr,number_of_p_points, finite_grid, inc_weights_in_pot, 
                cut_on_shell, program_const,obs_string,print);
    }
    
    // Check total cross sections
    
    std::vector<std::string> obs_string2 = {"SGT","SGTL","SGTT"};
    for (auto obs_string : obs_string2)
    {
        check_observable_noang(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot,
                rel_corr,number_of_p_points, finite_grid, inc_weights_in_pot, 
                cut_on_shell, program_const,obs_string,print);
    }
    

    // Check M-matrix elements
    check_M_el(J_max, Lambda, cut_pow, sharp_cutoff, pre_comp_pot,
            rel_corr,number_of_p_points, finite_grid, inc_weights_in_pot, 
            cut_on_shell, program_const,print);
    
    
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
    std::cout << "Testing observables with the nijmegen1 potential" << std::endl;
    std::cout << "------------------------------------------------\n\n"; 

    nn_mwpc_interface interface = nn_mwpc_interface("nijmegen1",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,program_const);   

    double energies1[3] = {10.0,50.0, 200.0};
    double energies2[3] = {125.0, 200.0,350.0};
    double* energies; 
    if (obs_string == "DSG")
    {
        energies = energies2;
    } else 
    {
        energies = energies1;
    }
    for (int j=0; j<3; j++)
    {
        double Tl = energies[j];
        std::cout << "Tlab=" << Tl << std::endl;
        //
        // Load data    
        //
        double theta_obs[180];
        double data_obs[180];
        double errors[180];

        std::string path;
        path = DATA_DIR + "np_" + obs_string 
            + "_" + std::to_string((int)Tl) + "_nijm1.txt";

        load_data(path,&theta_obs[0],&data_obs[0]);
        //
        // Open outfile
        //
        std::ofstream myfile;
        std::string filename;

        path = OUT_DIR + "out_" + obs_string + "_" + std::to_string((int)Tl) + ".txt"; 
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

        interface.solve_LS_ext_pot(Tl);
        std::cout << "Angle \t obs \t correct \t abs. rel. error" << std::endl;    
        for (int i=0; i<180; i++)
        {
            double theta_deg = theta_obs[i];
            double obs   = interface.compute_observable(obs_string2,theta_deg);
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
}

void check_M_el(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, bool print)
{
    std::cout << "Testing observables with the nijmegen1 potential" << std::endl;
    std::cout << "------------------------------------------------\n\n"; 

    nn_mwpc_interface interface = nn_mwpc_interface("nijmegen1",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,program_const);   

    double energies[3] = {10.0,50.0, 200.0};
    for (int j=0; j<3; j++)
    {
        double Tl = energies[j];
        std::cout << "Tlab=" << Tl << std::endl;
        //
        // Load data    
        //
        double theta[180];
        double data_real[180];
        double data_imag[180];
        double err_real[180];
        double err_imag[180];
        
        interface.solve_LS_ext_pot(Tl);

        std::string path;
        for (int S=0; S<2; S++)
        {
        for (int Mo=0; Mo<S+1; Mo++)
        {
        for (int Mi=0; Mi<S+1; Mi++)
        {
            if (S==1)
            {
                path = DATA_DIR + "np_M_" +  
                    std::to_string(Mo) + std::to_string(Mi)
                    + "_" + std::to_string((int)Tl) + "_nijm1.txt";
            } else 
            {
                path = DATA_DIR + "np_M_ss_" + 
                    std::to_string((int)Tl) + "_nijm1.txt";
            }
            load_data_M(path,&theta[0],&data_real[0],&data_imag[0]);
            if (print)
            {
                std::cout << 
                    "Angle \t th_real \t th_imag \t abs. error (real) \t (imag)" 
                    << std::endl;   
            }
            for (int i=0; i<179; i++)
            {
                std::complex<double>
                    M_el = interface.compute_M_element(Tl,
                    theta[i],S,Mo,Mi);
                
                M_el *= std::sqrt(program_const->MeVm2_to_mbarn);
                if (data_real[i] != 0 && data_imag[i] != 0) {
                    err_real[i] = std::abs(data_real[i] - std::real(M_el));
                    err_imag[i] = std::abs(data_imag[i] - std::imag(M_el));
                } else {
                    err_real[i] = 0;
                    err_imag[i] = 0;
                }
                if (print)
                {
                    std::cout << theta[i] << "\t" << std::real(M_el) << "\t" 
                        << std::imag(M_el) << "\t" << 
                        err_real[i] << "\t" << err_imag[i] << std::endl;
                }

            }
            std::cout << "Maximum error: " << 
                *(std::max_element(err_real,err_real + 179)) << "," << 
                *(std::max_element(err_imag,err_imag + 179)) 
                << "\n\n";
        }
        }
        }
        std::cout << "\n\n";
    }
}

void check_observable_noang(
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, std::string obs_string,bool print)
{
    std::cout << "Testing observables with the nijmegen1 potential" << std::endl;
    std::cout << "------------------------------------------------\n\n"; 

    nn_mwpc_interface interface = nn_mwpc_interface("nijmegen1",J_max,
         Lambda,cut_pow,sharp_cutoff,pre_comp_pot,rel_corr,number_of_p_points,
         finite_grid, inc_weights_in_pot, cut_on_shell,program_const);   

    //
    // Load data    
    //
    double energies[3];
    double data_obs[3];
    double errors[3];

    std::string path;
    path = DATA_DIR + "np_" + obs_string + "_nijm1.txt";

    load_data(path,&energies[0],&data_obs[0]);
    //
    // Open outfile
    //
    std::ofstream myfile;
    std::string filename;

    path = OUT_DIR + "out_" + obs_string + ".txt"; 
    myfile.open(path);
    for (int j=0; j<3; j++)
    {
        double Tl = energies[j];
        std::cout << "Tlab=" << Tl << std::endl;
        
            
        interface.solve_LS_ext_pot(Tl);
        std::cout << "Tlab \t obs \t correct \t abs. rel. error" << std::endl;    
            
        double obs   = interface.compute_observable(obs_string,0.0);
        if (data_obs[j] != 0) {
            errors[j] = std::abs((data_obs[j] - obs)/data_obs[j]);
        } else {
            errors[j] = 0;
        }
        myfile << Tl << "\t" << obs << "\t" << data_obs[j] << "\t" << 
            errors[j] << std::endl;
        std::cout << Tl << "\t" << obs << "\t" << data_obs[j] << "\t" << 
            errors[j] << std::endl;
    }
    
    myfile.close();
    double mean_error = 0;
    for (int i=0; i<3; i++)
    {  
        mean_error += errors[i];
    }
    mean_error /= 3.0;
    std::cout << "Mean absolute relative error: " << mean_error/3.0  
        << std::endl;
    std::cout << "Maximum error: " << *(std::max_element(errors, errors + 3)) 
        << std::endl;
}

void check_binding( 
        unsigned int J_max,double Lambda, int cut_pow,bool sharp_cutoff,
        bool pre_comp_pot, bool rel_corr,unsigned int number_of_p_points, bool 
        finite_grid, bool inc_weights_in_pot, bool cut_on_shell, 
        ph::constants_struct* program_const, int chn_number, bool print)
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Testing binding energy with the nijmegen1 potential." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    
    double* p_grid;
    double* w_grid;
    bool FINITE_GRID = false;
    double E_nijmegen = -2.224575;
    Lambda = 10000.0;
    double scale = 100.0;
    
    std::vector<qs::quantum_NN_state> states = 
        get_states_NN(J_max, 0, 0, 0, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    qs::quantum_channel chn = chns[3]; // 3S1-3D1
    
    ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    
    Potential_ext<gsl_matrix> nijmegen = Potential_ext<gsl_matrix>(
            p_grid, number_of_p_points, Lambda, &nijm_correct_arg);
    
    gsl_matrix* pot_V_mtx = nijmegen.get_matrix_no_onshell(chn,false);
    ph::eigen_t diag_res = ph::solve_SE(p_grid, w_grid, 
            number_of_p_points, chn, pot_V_mtx,program_const->Mn,program_const->Mp);
    
    std::cout << "The eigenvalues in " << quantum_channel_to_string(chn) << " are: "
        << std::endl;
    std::cout << "Nijmegen1 binding energy: -2.224575 MeV" << std::endl;
    
    for (int i = 0; i < (int)diag_res.eigenvalues->size; i++)
    {
        double E = GSL_REAL(gsl_vector_complex_get(diag_res.eigenvalues,i));
        if (E < 0.0)
        {
            std::cout << "This code: " << E << " MeV" << std::endl;
            std::cout << "Relative error: " << (E-E_nijmegen)/E_nijmegen << std::endl; 
        }
    }
    gsl_matrix_complex_free(diag_res.eigenvectors);
    gsl_vector_complex_free(diag_res.eigenvalues); 

    delete[] pot_V_mtx;
}

void load_data(std::string path,double* theta_obs,double* data_obs)
{
    std::ifstream infile(path);
    if (infile.is_open())
    {
        std::cout << "File" + path + " loaded: OK" << std::endl;
    } else
    {
        std::cout << "File" + path + ": Failed" << std::endl;
    }

    // Read and save the data to arrays
    double theta, obs;
    int k = 0;
    while(infile >> theta >> obs)
    {
        theta_obs[k] = theta;
        data_obs[k]  = obs;
        k++;
    }
}

void load_data_M(std::string path,double* theta,double* data_real,
        double* data_imag)
{
    std::ifstream infile(path);
    if (infile.is_open())
    {
        std::cout << "File" + path + " loaded: OK" << std::endl;
    } else
    {
        std::cout << "File" + path + ": Failed" << std::endl;
    }

    // Read and save the data to arrays
    double th, real, imag;
    int k = 0;
    std::string str;
    std::getline(infile,str);
    std::getline(infile,str);
    std::getline(infile,str);
    while(infile >> th >> real >> imag)
    {
        theta[k] = th;
        data_real[k]  = real;
        data_imag[k]  = imag;
        k++;
    }
}

arg_struct parse_arguments(int argc, char** argv)
{

    arg_struct ar;
    ar.args = {"-t","-f","-p","-h"};
    ar.arg_value = {"phase","./test_logfiles/out_default.txt","false","false"};

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
    myfile << "New log" << std::endl;
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
