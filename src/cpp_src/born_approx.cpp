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

gsl_matrix_complex* dwba::pw_T_BA(int start_order,int stop_order, gsl_matrix_complex* V, 
        gsl_matrix_complex* G0)
{
    gsl_matrix_complex* G0V = gsl_matrix_complex_alloc(V->size1,V->size2);

    ph::on_shell_mult(G0,V,G0V);
    
    gsl_matrix_complex* tmp1 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* tmp2 = gsl_matrix_complex_alloc(V->size1,V->size2);
    gsl_matrix_complex* res = gsl_matrix_complex_alloc(V->size1,V->size2);

    gsl_matrix_complex_set_zero(res);

    for (int i=start_order; i<stop_order+1;i++)
    {
        // Compute (G0*V)^i
        pow_matrix_on_shell_mult(G0V,i,tmp1);
        ph::on_shell_mult(V,tmp1,tmp2);

        // Add to the result
        gsl_matrix_complex_add(res,tmp2);
            
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
        return T_I;
    }
    // Get the Möller wave operators
    gsl_matrix_complex* omega_p = dwba::pw_moller_plus(T_I, G0);
    gsl_matrix_complex* omega_m_dagger = dwba::pw_moller_minus_dagger(T_I, G0);
    
    
    // Allocate the modified matrix
    gsl_matrix_complex* omega_p_G0 = gsl_matrix_complex_alloc(T_I->size1,
            T_I->size2);
        
    ph::on_shell_mult(omega_p,G0,omega_p_G0);

    // Perform the sum with the help of the Born approx code
    int start = 1;
    
    gsl_matrix_complex* VGV = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    pw_T_BA(start,order,V_II,omega_p_G0);
    
    // Multiply the VGVGV...V sum with the Möller operators from left and right
    // This is what I call F(...) in the notes
    
    gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    F(omega_p,omega_m_dagger,VGV,tmp);

    // Add the result to the leading order T_I matrix
    // T = T_I + tmp
    gsl_matrix_complex_add(tmp,T_I); // tmp <- tmp + T_I

    // Remove all temporary matrices
    gsl_matrix_complex_free(omega_p_G0);
    gsl_matrix_complex_free(omega_p);
    gsl_matrix_complex_free(omega_m_dagger);
    gsl_matrix_complex_free(VGV);

    return tmp;
}

gsl_matrix_complex* pw_moller_plus(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0)
{
    // Omega_p = 1 + G0*T
    //
    // Create an identity matrix
    gsl_matrix_complex* id = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    gsl_matrix_complex_set_identity(id);
    
    // Multiply G0*T
    gsl_matrix_complex* omega_p = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    ph::on_shell_mult(G0,T_I,omega_p);
    
    // Add them
    gsl_matrix_complex_add(omega_p,id);
    
    gsl_matrix_complex_free(id);
    return omega_p;
}


gsl_matrix_complex* pw_moller_minus_dagger(gsl_matrix_complex* T_I, 
        gsl_matrix_complex* G0)
{
    // Omega^dagger_m = 1 + T*G0
    //
    // Create an identity matrix
    gsl_matrix_complex* id = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    gsl_matrix_complex_set_identity(id);
    
    // Multiply G0*T
    gsl_matrix_complex* omega_p = gsl_matrix_complex_alloc(T_I->size1,T_I->size2);
    ph::on_shell_mult(T_I,G0,omega_p);
    
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
    } else 
    {
        gsl_matrix_complex* tmp = gsl_matrix_complex_alloc(M->size1,M->size2);
        gsl_matrix_complex_set_identity(tmp);
        
        gsl_matrix_complex_memcpy(res,M);
        for (int i = 0; i < pow; i++)
        {
            ph::on_shell_mult(res,M,tmp);
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
    ph::on_shell_mult(omega_m_dagger,M,tmp);
    // right
    ph::on_shell_mult(tmp,omega_p,res);

    gsl_matrix_complex_free(tmp);   
}


void dwba::make_tests(std::string chn_string)
{

    std::cout << "Testing the Born and DW Born approximation" << std::endl;
    
    

    // ------ CONSTANTS TO CHANGE ------
    // ---------------------------------
    double scale = 100.0; // Scale of momenutm grid MeV
    unsigned int ang_int_points = 76; // Number of points in angular integration
    unsigned int number_of_p_points = 120; // Number of momentum-grid points
    unsigned int J_max_in_pot = 50; // Maximum J that is stored for L-polynomials
    
    // Do precomputations
    ph::physics_helpers_init();
    // ---------------   
    
    // Construct the quantum states
    std::cout << "Constructing quantum states..." << std::endl;
    int J_max = 8;
    int J_min = 0;
    int Tz_min = 0;
    int Tz_max = 0;
    bool print = false;
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
    double Lambda = 450.0;
    bool FINITE_GRID = true;
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
    int cut_pow = 6;
    double C1S0	= -0.1/100.0; // contact term C1S0 for lambda = 450 [MeV]
    double C3S1	= -0.13/100.0; // contact term C3S1 for lambda = 450 [MeV]
    double C3P0 = 5e-8;
    double C3P2 = 0.0;    
    // Choose terms in LO MWPC potential
    std::vector<std::string> terms;
    terms.push_back("OPEP"); // To just test elements use just OPEP
    terms.push_back("C1S0");
    terms.push_back("C3S1");
    terms.push_back("C3P0");
    terms.push_back("C3P2");

    // Without the grid included
    Potential_mwpc Pot = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,false,true);
    // With grid included
    Potential_mwpc Pot_grid = Potential_mwpc(terms,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,true,true);

    // Potential with only the 3P0 LEC
    std::vector<std::string> terms_3P0;
    terms_3P0.push_back("C3P0");
    
    Potential_mwpc Pot_3P0 = Potential_mwpc(terms_3P0,ang_int_points,p_grid,w_grid,
            number_of_p_points,J_max_in_pot,Lambda, cut_pow, false,true,true);
    std::cout << "Saving potential matrices" << std::endl;
    for (auto chn : chns)
    {
        Pot.populate_saved_mtx(chn,true); // Realtivistic factor on
        Pot_grid.populate_saved_mtx(chn,true); // Realtivistic factor on
        Pot_3P0.populate_saved_mtx(chn,true); // Realtivistic factor on
    }

    // Set correct LECs
    Pot.LECs_["gA2"]  = 1.29*1.29;
    Pot.LECs_["C1S0"] = C1S0;
    Pot.LECs_["C3S1"] = C3S1;
    Pot.LECs_["C3P0"] = C3P0;
    Pot.LECs_["C3P2"] = C3P2;

    Pot_grid.LECs_["gA2"]  = 1.29*1.29;
    Pot_grid.LECs_["C1S0"] = C1S0;
    Pot_grid.LECs_["C3S1"] = C3S1;
    Pot_grid.LECs_["C3P0"] = 0;
    Pot_grid.LECs_["C3P2"] = C3P2;

    Pot_3P0.LECs_["C3P0"] = C3P0;
    
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
    
    double Tl = 30.0; // MeV
        
    // Get the on-shell momenta and reduced mass
    LS_Solver::get_mu_q_on_shell(Tl, chn, &mu, &q_on_shell);

    gsl_matrix* V_Pot  = Pot.get_saved_matrix(q_on_shell, chn, true);
    gsl_matrix* V_grid = Pot_grid.get_saved_matrix(q_on_shell, chn, true);
    gsl_matrix* V_3P0  = Pot_3P0.get_saved_matrix(q_on_shell, chn, true);
    
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
    V_Pot  = Pot.get_saved_matrix(q_on_shell, chn, true);
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

    gsl_matrix_complex_free(T_dwba);
    gsl_matrix_complex_free(V_pot_z);
    gsl_matrix_complex_free(V_grid_z);
    gsl_matrix_complex_free(V_3P0_z);
    gsl_matrix_complex_free(G0);
    gsl_vector_complex_free(prop_vec);
    gsl_matrix_complex_free(T_full);
    gsl_matrix_free(V_Pot);
    gsl_matrix_free(V_grid);
    gsl_matrix_free(V_3P0);
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
            
            if (M->size1 == mom_grid_size+1)
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
    // Set the diagonal values to the vecotr values
    for (int i = 0; i < M->size1; i++)
    {
        gsl_matrix_complex_set(M,i,i,gsl_vector_complex_get(vec,i));
    }

}

