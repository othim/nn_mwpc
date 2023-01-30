/*
 * born_approx.cpp
 * This file is part of the nn_mwpc project.
 *
 * This file contatins the implementation of the function in the header
 * born_approx.h.
 *
 * Oliver Thim 2022-11 --
 * Department of Physics, Chalmers
 */
#include "born_approx.h"

void dress_in_weights(gsl_matrix_complex* M,double* p,double* w,
        int mom_grid_size);

void matrix_from_vector(gsl_matrix_complex* M,gsl_vector_complex* vec);

void pow_matrix_mult(gsl_matrix_complex* m,int pow,gsl_matrix_complex* res);

gsl_matrix_complex* dwba::pw_T_BA(int start_order,int stop_order, gsl_matrix_complex* V, 
        gsl_matrix_complex* G0)
{
    bool print = false;
    if (print)
    {
        std::cout << "V:" << std::endl;
        ph::print_m(V);
    
        std::cout << "G0:" << std::endl;
        ph::print_m(G0);
    }
    gsl_matrix_complex* G0V = gsl_matrix_complex_alloc(V->size1,V->size2);
    ph::mult(G0,V,G0V);
    
    if (print)
    {
        std::cout << "G0V:" << std::endl;
        ph::print_m(G0V);
    }

    gsl_matrix_complex* tmp1 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* tmp2 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* res = gsl_matrix_complex_alloc(V->size1,V->size2);

    gsl_matrix_complex_set_zero(res);

    for (int i=start_order; i<stop_order+1;i++)
    {
        // Compute (G0*V)^i
        if (i==0)
        {
            gsl_matrix_complex_add(res,V);
        } else
        {
            pow_matrix_mult(G0V,i,tmp1);
            
            //std::cout << "(G0V)^i:" << std::endl;
            //ph::print_m(tmp1);
            
            ph::mult(V,tmp1,tmp2);
            
            //std::cout << "CGV:" << std::endl;
            //ph::print_m(tmp2);
            // Add to the result
            gsl_matrix_complex_add(res,tmp2);
            
        }
            
        //std::cout << "res:" << std::endl;
        //ph::print_m(res);
    }
    // Deallocate 
    gsl_matrix_complex_free(G0V);
    gsl_matrix_complex_free(tmp1);
    gsl_matrix_complex_free(tmp2);

    return res;
}

gsl_matrix_complex* dwba::pw_T_DWBA(int order, 
        gsl_matrix_complex* T_I, gsl_matrix_complex* V_I,
        gsl_matrix_complex* V_II, gsl_matrix_complex* G0)
{
    // If the order is trivial
    if (order==0) 
    {
        gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
        gsl_matrix_complex_memcpy(tmp,T_I);
        return tmp;
    }
    // Get the Möller wave operators
    gsl_matrix_complex* omega_p = dwba::pw_moller_plus(T_I, G0);
    gsl_matrix_complex* omega_m_dagger = dwba::pw_moller_minus_dagger(T_I, G0);
    
    
    // Allocate the G1 matrix
    gsl_matrix_complex* G1 = gsl_matrix_complex_alloc(T_I->size1,
            T_I->size2);
    
    // G1 = G0 + G0*T_I*G0 = (1 + G0*T_I)*G0
    ph::mult(omega_p,G0,G1);
    //ph::mult(G0,omega_m_dagger,G1);

    // Perform the sum with the help of the Born approx code
    // VGV_sum = V_II + V_II*G1*V_II + ...
    gsl_matrix_complex* VGV_sum = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);


    // Call the BA code with the V_II potential ans the full propagator
    // to compute V_II + V_II*G1*V_II + ...
    int start = 0;
    int stop = order-1;
    VGV_sum = pw_T_BA(start,stop,V_II,G1);

    // Multiply the VGVGV...V sum with the Möller operators from left and right
    // This is what I call F(...) in the notes
    
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    F(omega_p,omega_m_dagger,VGV_sum,tmp);

    // Add the result to the leading order T_I matrix
    // T = T_I + tmp
    gsl_matrix_complex_add(tmp,T_I); // tmp <- tmp + T_I

    // Remove all temporary matrices
    gsl_matrix_complex_free(G1);
    gsl_matrix_complex_free(omega_p);
    gsl_matrix_complex_free(omega_m_dagger);
    gsl_matrix_complex_free(VGV_sum);

    return tmp;
}

gsl_matrix_complex* dwba::pw_moller_plus(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0)
{
    // Omega_p = 1 + G0*T_I
    
    // Create an identity matrix
    gsl_matrix_complex* id = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    gsl_matrix_complex_set_identity(id);
    
    //gsl_matrix_complex_set(id,id->size1-1,id->size2-1,gsl_complex_rect(0.0,0.0));
    
    // Multiply G0*T
    gsl_matrix_complex* omega_p = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    ph::mult(G0,T_I,omega_p);
    
    // Add them
    gsl_matrix_complex_add(omega_p,id);
    
    gsl_matrix_complex_free(id);
    return omega_p;
}


gsl_matrix_complex* dwba::pw_moller_minus_dagger(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0)
{
    // Omega^dagger_m = 1 + T_I*G0
    
    // Create an identity matrix
    gsl_matrix_complex* id = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    gsl_matrix_complex_set_identity(id);
    
    //gsl_matrix_complex_set(id,id->size1-1,id->size2-1,gsl_complex_rect(0.0,0.0));

    // Multiply G0*T
    gsl_matrix_complex* omega_p = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    ph::mult(T_I,G0,omega_p);
    
    // Add them
    gsl_matrix_complex_add(omega_p,id);
    
    gsl_matrix_complex_free(id);

    return omega_p;
}

/*
 *
 * Methods outside the namespace that is conisidered helper methods
 * that should not generally be acessed from outside of this file.
 *
 */

void pow_matrix_on_shell_mult(gsl_matrix_complex* M,int pow,gsl_matrix_complex* res)
{
    if (pow==0)
    {
        gsl_matrix_complex_set_identity(res);
        return;
    } else if (pow==1)
    {
        gsl_matrix_complex_memcpy(res,M);
        return;
    } 
    else 
    {
        gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);
        gsl_matrix_complex_set_identity(tmp);
        
        gsl_matrix_complex_memcpy(res,M);
        for (int i = 0; i < pow-1; i++)
        {
            ph::on_shell_mult(res,M,tmp);
            gsl_matrix_complex_memcpy(res,tmp);

        }
        gsl_matrix_complex_free(tmp);
        return;
    }
}

void pow_matrix_mult(gsl_matrix_complex* M,int pow,gsl_matrix_complex* res)
{
    if (pow==0)
    {
        gsl_matrix_complex_set_identity(res);
        return;
    } else if (pow==1)
    {
        gsl_matrix_complex_memcpy(res,M);
        return;
    } 
    else 
    {
        gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);
        gsl_matrix_complex_set_identity(tmp);
        
        gsl_matrix_complex_memcpy(res,M);
        for (int i = 0; i < pow-1; i++)
        {
            ph::mult(res,M,tmp);
            gsl_matrix_complex_memcpy(res,tmp);

        }
        gsl_matrix_complex_free(tmp);
        return;
    }
}

void F(gsl_matrix_complex* omega_p,gsl_matrix_complex* omega_m_dagger,
        gsl_matrix_complex* M,gsl_matrix_complex* res)
{
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);

    // left
    ph::mult(omega_m_dagger,M,tmp);
    // right
    ph::mult(tmp,omega_p,res);

    gsl_matrix_complex_free(tmp);   
}


void dwba::make_tests(std::string chn_string)
{

    std::cout << "Testing the Born and DW Born approximation" << std::endl;
    
    

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 76; // Number of points in angular integration
    unsigned int number_of_p_points = 1000; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    bool REL_CORR = false;
    bool CUT_ON_SHELL = true;
    
    int J_max = 2;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = false;
    
    int cut_pow = 10000000;
    double C1S0	= -0.01/100.0; // contact term C1S0 for lambda = 450 [MeV]
    
    double Lambda = 450.0;
    bool FINITE_GRID = true;
    
    double Tl = 1.0; // MeV
    // ---------------------------------
    // ---------------------------------
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    
    // ---------------------------------
    // ---------------------------------
    
    /*
     * Construct the momentum grid
     */
    double* p_grid;
    double* w_grid;
    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,
                Lambda + 300.0,&p_grid,&w_grid);
    } else 
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }

    /*
     * Construct the potential
     */
    
    std::vector<std::string> terms_1S0;
    terms_1S0.push_back("C1S0");
    
    Potential_mwpc Pot_1S0 = Potential_mwpc(terms_1S0,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,true,CUT_ON_SHELL);
    
    Potential_mwpc Pot_1S0_nogrid = Potential_mwpc(terms_1S0,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,false,CUT_ON_SHELL);
    
    
    std::vector<std::string> terms;
    terms.push_back("Yamaguchi_1S0");
    Potential_mwpc Pot_Yam = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,true,CUT_ON_SHELL);
    
    Potential_mwpc Pot_Yam_nogrid = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,false,CUT_ON_SHELL);
    
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot_1S0.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
        Pot_1S0_nogrid.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
        Pot_Yam.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
        Pot_Yam_nogrid.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
    }

    
    Pot_1S0.LECs_["C1S0"] = C1S0;
    Pot_1S0_nogrid.LECs_["C1S0"] = C1S0;

    Pot_Yam.LECs_["Yamaguchi_1S0"] = 100.0;
    Pot_Yam_nogrid.LECs_["Yamaguchi_1S0"] = 100.0;
    
    // Solve for the T-matrix
    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    //double rho_T
    qs::quantum_channel chn = chns[0];

    // Get the channel
    for (int i = 0; i < (int)chns.size(); i++)
    {
        chn = chns[i];
        if (quantum_channel_to_string(chn) == chn_string) {
            break;
        }
    }

    // Print 
    std::cout << "Computing in chn: " << quantum_channel_to_string(chn) << std::endl;
    // Solve the distorted wave problem
    
        
    // Get the on-shell momenta and reduced mass
    LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
    std::cout << std::setprecision(16) << "mu = " << mu << " q_on_shell = " << q_on_shell << std::endl;
    gsl_matrix* V_1S0  = Pot_1S0.get_saved_matrix(q_on_shell, chn, REL_CORR);
    gsl_matrix* V_1S0_nogrid  = Pot_1S0_nogrid.get_saved_matrix(q_on_shell, chn, REL_CORR);
    //gsl_matrix* V_Yam_nogrid  = Pot_Yam_nogrid.get_saved_matrix(q_on_shell, chn, REL_CORR);
    
    
    std::cout << "Solving exact 1S0" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::complex<double>* tt = solver.solve_in_chn_T_Telem(Tl,chn,V_1S0_nogrid);
    std::cout << "T (exact) = " << tt[3] << std::endl;
    std::cout << "|T|^2 (exact) = " << std::pow(std::abs(tt[3]),2) << std::endl;
    std::cout << "-----------------" << std::endl;
    
    std::cout << "Making potential complex" << std::endl;
    gsl_matrix_complex* V_1S0_z  = gsl_matrix_complex_alloc(V_1S0->size1,V_1S0->size2);
    ph::make_matrix_complex(V_1S0_z,V_1S0);
    
    
    //V_1S0  = Pot_Yam.get_saved_matrix(q_on_shell, chn, REL_CORR);
    //ph::make_matrix_complex(V_1S0_z,V_1S0);
    
    // Get the propagator matrix
    std::cout << "Computing the propagator" << std::endl;
    gsl_vector_complex* prop_vec = solver.setup_G0_vector_complex(q_on_shell,
            chn.coupled,mu);
    
    std::cout << "-----------------" << std::endl;
    std::cout << "-----------------" << std::endl;
    gsl_matrix_complex* G0 = gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    matrix_from_vector(G0,prop_vec);
 
    
    std::ofstream myfile;
    //std::string DATA_DIR = "~/Documents/phd/projects/dwb/data/";
    std::string DATA_DIR = "../../../projects/dwb/data/";
    std::string filename = DATA_DIR +"born_approx.txt"; 
    myfile.open(filename);
    myfile << "Np =" << number_of_p_points << std::endl;
    myfile << "Order, |T|^2" << std::endl;
    for (int ord = 0; ord < 11; ord++)
    {
        gsl_matrix_complex* T_ba = dwba::pw_T_BA((int)0,ord,V_1S0_z,G0);
        gsl_complex onT = gsl_matrix_complex_get(T_ba,number_of_p_points,number_of_p_points);
        std::cout << "order=" << ord << ", in BA" << std::endl; 
        std::cout << "T= " << GSL_REAL(onT) << "," << GSL_IMAG(onT) << std::endl;
        double T2 = std::pow(GSL_REAL(onT),2)+ std::pow(GSL_IMAG(onT),2);
        std::cout <<"|T|^2 = " <<  T2 << std::endl << std::endl;
        
        // Write to file in format order, |T|^2
        myfile << ord << "   " << T2 << std::endl;
    }
    myfile.close();

    gsl_matrix_complex_free(V_1S0_z);
    gsl_matrix_complex_free(G0);
    gsl_vector_complex_free(prop_vec);
    gsl_matrix_free(V_1S0);
}


void dwba::make_tests_DWBA(std::string chn_string)
{

    std::cout << "Testing the Born and DW Born approximation" << std::endl;
    

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 76; // Number of points in angular integration
    unsigned int number_of_p_points =100; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    bool REL_CORR = false;
    bool CUT_ON_SHELL = true;
    
    int J_max = 2;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = false;
    
    int cut_pow = 100000000;
    //double C1S0	= -0.01/100.0; // contact term C1S0 for lambda = 450 [MeV]
    
    double Lambda = 1000000.0;
    bool FINITE_GRID = false;
    
    double Tl = 1.0; // MeV
    
    double lam   = 100.0;
    double lam_t = -10000.0;
    double beta  = 40.0;
    double gamma = 30.0;
    //double lam_t = 0.0;
    std::string DATA_DIR = "../../../projects/dwb/data/";
    // ---------------------------------
    // ---------------------------------
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    std::vector<qs::quantum_NN_state> states = get_states_NN(J_max, J_min, Tz_min, Tz_max, print);
     
    // Construct the quantum scattering channels from the states
    std::cout << "Contruction scattering channels..." << std::endl;
    std::vector<qs::quantum_channel> chns = get_channels(states, print);   
    
    // ---------------------------------
    // ---------------------------------
    
    /*
     * Construct the momentum grid
     */
    double* p_grid;
    double* w_grid;
    if (FINITE_GRID)
    {
        ph::gauss_legendre_finite_mesh(number_of_p_points,0,
                Lambda + 300.0,&p_grid,&w_grid);
    } else 
    {
        ph::gauss_legendre_inf_mesh(number_of_p_points,scale,&p_grid,&w_grid);
    }

    /*
     * Construct the potential
     */
    
    
    std::vector<std::string> terms;
    terms.push_back("Yamaguchi_1S0");
    Potential_mwpc Pot_Yam = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,true,CUT_ON_SHELL);
    
    Potential_mwpc Pot_Yam_nogrid = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,false,CUT_ON_SHELL);
    
    Pot_Yam_nogrid.params_["Yamaguchi_beta"] = beta;
    Pot_Yam.params_["Yamaguchi_beta"] = gamma;

    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot_Yam.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
        Pot_Yam_nogrid.populate_saved_mtx(chn,REL_CORR); // Realtivistic factor on
    }
    
    // Solve for the T-matrix
    LS_Solver solver = LS_Solver(number_of_p_points,p_grid,w_grid,FINITE_GRID);
   
    double q_on_shell;
    double mu;
    //double rho_T
    qs::quantum_channel chn = chns[0];

    // Get the channel
    for (int i = 0; i < (int)chns.size(); i++)
    {
        chn = chns[i];
        if (quantum_channel_to_string(chn) == chn_string) {
            break;
        }
    }

    // Print 
    std::cout << "Computing in chn: " << quantum_channel_to_string(chn) << std::endl;
    // Solve the distorted wave problem
    
        
    // Get the on-shell momenta and reduced mass
    LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);
    std::cout << std::setprecision(16) << "mu = " << mu << " q_on_shell = " << q_on_shell << std::endl;
    
    
    Pot_Yam_nogrid.LECs_["Yamaguchi_1S0"] = lam_t;
    Pot_Yam_nogrid.params_["Yamaguchi_beta"] = beta;
    gsl_matrix* V_Yam_nogrid  = Pot_Yam_nogrid.get_matrix(q_on_shell, chn, REL_CORR);
    
    Pot_Yam_nogrid.LECs_["Yamaguchi_1S0"] = lam;
    Pot_Yam_nogrid.params_["Yamaguchi_beta"] = gamma;
    gsl_matrix* tmp  = Pot_Yam_nogrid.get_matrix(q_on_shell, chn, REL_CORR);
    
    gsl_matrix_add(V_Yam_nogrid,tmp);

    std::cout << "Solving exact 1S0" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::complex<double>* tt = solver.solve_in_chn_T_Telem(Tl,chn,V_Yam_nogrid);
    //gsl_matrix_free(V_Yam_nogrid);
    std::cout << "T (exact) = " << tt[3] << std::endl;
    std::cout << "|T|^2 (exact) = " << std::pow(std::abs(tt[3]),2) << std::endl;
    std::cout << "-----------------" << std::endl;
    
    std::cout << "Making potential complex" << std::endl;
    gsl_matrix_complex* V_1S0_I  = gsl_matrix_complex_alloc(V_Yam_nogrid->size1,V_Yam_nogrid->size2);
    gsl_matrix_complex* V_1S0_II  = gsl_matrix_complex_alloc(V_Yam_nogrid->size1,V_Yam_nogrid->size2);
    
    Pot_Yam.params_["Yamaguchi_beta"] = beta;
    Pot_Yam.LECs_["Yamaguchi_1S0"] = lam_t;
    gsl_matrix* V_Yam_I  = Pot_Yam.get_matrix(q_on_shell, chn, REL_CORR);
    ph::make_matrix_complex(V_1S0_I,V_Yam_I);
    ph::print_m_complex_to_file(DATA_DIR+"VI.txt",V_1S0_I);
    
    Pot_Yam.params_["Yamaguchi_beta"] = gamma;
    Pot_Yam.LECs_["Yamaguchi_1S0"] = lam;
    gsl_matrix* V_Yam_II  = Pot_Yam.get_matrix(q_on_shell, chn, REL_CORR);
    ph::make_matrix_complex(V_1S0_II,V_Yam_II);
    ph::print_m_complex_to_file(DATA_DIR+"VII.txt",V_1S0_II);
    
    // Get the propagator matrix
    std::cout << "Computing the propagator" << std::endl;
    gsl_vector_complex* prop_vec = solver.setup_G0_vector_complex(q_on_shell,
            chn.coupled,mu);
    
    std::cout << "-----------------" << std::endl;
    std::cout << "-----------------" << std::endl;
    gsl_matrix_complex* G0 = gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    matrix_from_vector(G0,prop_vec);
 
    // Store some matrices to file
    ph::print_m_complex_to_file(DATA_DIR+"G0.txt",G0);
    
    std::ofstream myfile;
    //std::string DATA_DIR = "~/Documents/phd/projects/dwb/data/";
    std::string filename = DATA_DIR +"DWBA"+std::to_string((int)number_of_p_points)+ ".txt"; 
    myfile.open(filename);
    myfile << "Np =" << number_of_p_points << std::endl;
    myfile << "Order, |T|^2" << std::endl;
    for (int ord = 0; ord < 11; ord++)
    {
        Pot_Yam_nogrid.params_["Yamaguchi_beta"] = beta;
        Pot_Yam_nogrid.LECs_["Yamaguchi_1S0"] = lam_t;
        gsl_matrix* V_Yam_nogrid  = Pot_Yam_nogrid.get_matrix(q_on_shell, chn, REL_CORR);
        gsl_matrix_complex* T_I = solver.solve_in_chn_T_fullT(Tl,chn,V_Yam_nogrid);
        dress_in_weights(T_I,p_grid,w_grid,(int)number_of_p_points);
        ph::print_m_complex_to_file(DATA_DIR+"TI.txt",T_I);
        gsl_complex onT_I = gsl_matrix_complex_get(T_I,number_of_p_points,number_of_p_points);
        std::cout << "T_I= " << GSL_REAL(onT_I) << "," << GSL_IMAG(onT_I) << std::endl;

        gsl_matrix_complex* T_DWBA = dwba::pw_T_DWBA(ord,T_I,V_1S0_I,V_1S0_II,G0);
        gsl_complex onT = gsl_matrix_complex_get(T_DWBA,number_of_p_points,number_of_p_points);
        
        gsl_matrix_complex_free(T_I);
        // Something is wrong with this memory...
        gsl_matrix_complex_free(T_DWBA); 
        gsl_matrix_free(V_Yam_nogrid);

        std::cout << "order=" << ord << ", in DWBA" << std::endl; 
        std::cout << "T= " << GSL_REAL(onT) << "," << GSL_IMAG(onT) << std::endl;
        double T2 = std::pow(GSL_REAL(onT),2)+ std::pow(GSL_IMAG(onT),2);
        std::cout <<"|T|^2 = " <<  T2 << std::endl << std::endl;
        
        // Write to file in format order, |T|^2
        myfile << ord << "   " << T2 << std::endl;
    }
    myfile.close();

    gsl_matrix_complex_free(V_1S0_I);
    gsl_matrix_complex_free(V_1S0_II);
    gsl_matrix_complex_free(G0);
    gsl_vector_complex_free(prop_vec);
    gsl_matrix_free(V_Yam_I);
    gsl_matrix_free(V_Yam_II);
}

/*
 * This function adds the weights and momenta to the matrix M
 */
void dress_in_weights(gsl_matrix_complex* M,double* p,double* w,
        int mom_grid_size)
{
    for (int i = 0; i<mom_grid_size + 1; i++) 
    {
        for (int j = 0; j<mom_grid_size + 1; j++) 
        {
            double weights_momenta = 1.0;
            // If ingoing momenta is of shell
            if (j<mom_grid_size)
            {
            weights_momenta *= std::sqrt(w[j])*p[j];
            }

            // If outgoing momenta is of shell
            if (i<mom_grid_size)
            {
            weights_momenta *= std::sqrt(w[i])*p[i];
            }
            
            if ((int)M->size1 == mom_grid_size+1)
            {   
                gsl_complex fac = gsl_complex_rect(weights_momenta,0.0);
                
                gsl_complex el  = gsl_matrix_complex_get(M,i,j);
                gsl_matrix_complex_set(M,i,j,gsl_complex_mul(el,fac));
            } else 
            {
                gsl_complex fac = gsl_complex_rect(weights_momenta,0.0);
                // Do all four elemets
                gsl_complex el_00  = gsl_matrix_complex_get(M,i,j);
                gsl_matrix_complex_set(M,i,j,gsl_complex_mul(el_00,fac));
                
                gsl_complex el_01  = gsl_matrix_complex_get(M,i,j+mom_grid_size+1);
                gsl_matrix_complex_set(M,i,j+mom_grid_size+1,gsl_complex_mul(el_01,fac));
                
                gsl_complex el_10  = gsl_matrix_complex_get(M,i+mom_grid_size+1,j);
                gsl_matrix_complex_set(M,i+mom_grid_size+1,j,gsl_complex_mul(el_10,fac));
                
                gsl_complex el_11  = gsl_matrix_complex_get(M,i+mom_grid_size+1,j+mom_grid_size+1);
                gsl_matrix_complex_set(M,i+mom_grid_size+1,j+mom_grid_size+1,gsl_complex_mul(el_11,fac));
            }
        }
    }
}


void matrix_from_vector(gsl_matrix_complex* M,gsl_vector_complex* vec)
{
    gsl_matrix_complex_set_zero(M);
    // Set the diagonal values to the vecotr values
    for (int i = 0; i < (int)M->size1; i++)
    {
        gsl_matrix_complex_set(M,i,i,gsl_vector_complex_get(vec,i));
    }
}


/*
 *
    // Get the on-shell momenta and reduced mass
    LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);

    gsl_matrix* V_Pot  = Pot.get_saved_matrix(q_on_shell, chn, REL_CORR);
    gsl_matrix* V_grid = Pot_grid.get_saved_matrix(q_on_shell, chn, REL_CORR);
    gsl_matrix* V_3P0  = Pot_3P0.get_saved_matrix(q_on_shell, chn, REL_CORR);
    gsl_matrix* V_1S0  = Pot_1S0.get_saved_matrix(q_on_shell, chn, REL_CORR);
    
    std::complex<double>* T_elem = solver.solve_in_chn_T_Telem(Tl, chn, V_Pot);
    double fac = 1.0;
    std::cout << "Full potential" << std::endl;
    std::cout << "T: " <<  fac*T_elem[0] << "   " << fac*T_elem[1] << 
            "   " << fac*T_elem[2] << "   " << fac*T_elem[3] << std::endl;
    
    
    // Do it in distorted wave parturbation theory
    // -------------------------------------------

    // First solve the problem where C3P0 is zero
    Pot.LECs_["C3P0"] = 0;
    gsl_matrix_free(V_Pot);
    V_Pot  = Pot.get_saved_matrix(q_on_shell, chn, REL_CORR);
    std::complex<double>* T_elem_0 = solver.solve_in_chn_T_Telem(Tl, chn, V_Pot);
    gsl_matrix_complex* T_full     = solver.solve_in_chn_T_fullT(Tl,chn,V_Pot);

    std::cout << "C3P0=0" << std::endl;
    std::cout << "T: " <<  fac*T_elem_0[0] << "   " << fac*T_elem_0[1] << 
            "   " << fac*T_elem_0[2] << "   " << fac*T_elem_0[3] << std::endl;
    
    
    std::cout << "Dressing T in the weights and momenta" << std::endl;
    
    dress_in_weights(T_full,p_grid,w_grid,(int)number_of_p_points);
    
    // Make potential complex (this is already dresses in the weights)
    std::cout << "Making potential complex" << std::endl;
    gsl_matrix_complex* V_pot_z  = gsl_matrix_complex_alloc(V_Pot->size1,V_Pot->size2);
    gsl_matrix_complex* V_grid_z = gsl_matrix_complex_alloc(V_Pot->size1,V_Pot->size2);
    gsl_matrix_complex* V_3P0_z  = gsl_matrix_complex_alloc(V_Pot->size1,V_Pot->size2);
    gsl_matrix_complex* V_1S0_z  = gsl_matrix_complex_alloc(V_Pot->size1,V_Pot->size2);

    // Get the propagator matrix
    std::cout << "Computing the propagator" << std::endl;
    gsl_vector_complex* prop_vec = solver.setup_D_vector_complex(q_on_shell,
            chn.coupled,mu);
    
    gsl_matrix_complex* G0 = gsl_matrix_complex_alloc(prop_vec->size,prop_vec->size);
    matrix_from_vector(G0,prop_vec);

    int order = 1;
    gsl_matrix_complex* T_dwba = pw_T_DWBA(order,T_full,V_grid_z,V_3P0_z,G0);

    gsl_complex onT = gsl_matrix_complex_get(T_dwba,number_of_p_points,number_of_p_points);

    std::cout << "order=" << order << ", in DWBA" << std::endl; 
    std::cout << GSL_REAL(onT) << "," << GSL_IMAG(onT) << std::endl;
    std::cout << "\n\n\n";
    
    for (int ord = 0; ord < 1; ord++)
    {
        gsl_matrix_complex* T_ba = dwba::pw_T_BA((int)0,ord,V_3P0_z,G0);
        onT = gsl_matrix_complex_get(T_ba,number_of_p_points,number_of_p_points);
        std::cout << "order=" << ord << ", in BA" << std::endl; 
        std::cout << GSL_REAL(onT) << "," << GSL_IMAG(onT) << std::endl;
        std::cout <<"|T|^2 = " <<  std::pow(GSL_REAL(onT),2)+ std::pow(GSL_IMAG(onT),2) << std::endl << std::endl;
    }
    std::cout << "C3P0 = " << C3P0 << std::endl;
    std::cout << "(2*pi)^3*C3P0 = " << std::pow((2*M_PI),3)*C3P0 << std::endl;

    std::cout << "|C3P0|^2 = " << std::pow(C3P0,2) << std::endl;
    std::cout << "(2*pi)^6*|C3P0|^2 = " << std::pow((2*M_PI),6)*std::pow(C3P0,2) << std::endl;


    gsl_matrix_complex_free(T_dwba);
    gsl_matrix_complex_free(V_pot_z);
    gsl_matrix_complex_free(V_grid_z);
    gsl_matrix_complex_free(V_3P0_z);
    gsl_matrix_complex_free(V_1S0_z);
    gsl_matrix_complex_free(G0);
    gsl_vector_complex_free(prop_vec);
    gsl_matrix_complex_free(T_full);
    gsl_matrix_free(V_Pot);
    gsl_matrix_free(V_grid);
    gsl_matrix_free(V_3P0);
}*/
