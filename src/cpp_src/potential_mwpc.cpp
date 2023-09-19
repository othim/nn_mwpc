#include "potential_mwpc.h"

/* These declarations are neccessary to be able to have the declarations
 * in a separate .cpp file.
 */
template void Potential_mwpc<gsl_matrix>::populate_saved_mtx(qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix* Potential_mwpc<gsl_matrix>::get_saved_matrix(double q_on_shell,qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix* Potential_mwpc<gsl_matrix>::get_matrix(double q_on_shell,qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix* Potential_mwpc<gsl_matrix>::get_matrix_no_onshell(
        qs::quantum_channel chn, bool rel_correction);

template Potential_mwpc<gsl_matrix>::~Potential_mwpc();

template Potential_mwpc<gsl_matrix>::Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,double* p_grid, 
   double* w_grid, std::size_t mom_grid_size, unsigned int J_max, double cutoff_Lambda, int cut_pow,
   bool sharp_cutoff, bool inc_grid_weights_in_pot, bool cut_on_shell,
   std::string loop_reg, double lam_SFR);

template void Potential_mwpc<gsl_matrix_complex>::populate_saved_mtx(qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix_complex* Potential_mwpc<gsl_matrix_complex>::get_saved_matrix(double q_on_shell,qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix_complex* Potential_mwpc<gsl_matrix_complex>::get_matrix(double q_on_shell,qs::quantum_channel chn, 
        bool rel_correction);
    
template gsl_matrix_complex* Potential_mwpc<gsl_matrix_complex>::get_matrix_no_onshell(
        qs::quantum_channel chn, bool rel_correction);

template Potential_mwpc<gsl_matrix_complex>::~Potential_mwpc();

template Potential_mwpc<gsl_matrix_complex>::Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,double* p_grid, 
   double* w_grid, std::size_t mom_grid_size, unsigned int J_max, double cutoff_Lambda, int cut_pow,
   bool sharp_cutoff, bool inc_grid_weights_in_pot, bool cut_on_shell,
   std::string loop_reg, double lam_SFR);

template void Potential_mwpc<gsl_matrix_complex>::print_LECs_and_params_info();
template void Potential_mwpc<gsl_matrix>::print_LECs_and_params_info();

template double Potential_mwpc<gsl_matrix>::get_rel_cut(double p_in, 
        double p_out, double mu, bool rel_correction);

template double Potential_mwpc<gsl_matrix_complex>::get_rel_cut(double p_in, 
        double p_out, double mu, bool rel_correction);

template <class gsl_m>
Potential_mwpc<gsl_m>::Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,double* p_grid, 
   double* w_grid, std::size_t mom_grid_size, unsigned int J_max, double cutoff_Lambda, int cut_pow,
   bool sharp_cutoff, bool inc_grid_weights_in_pot, bool cut_on_shell,
   std::string loop_reg, double lam_SFR)
{
   // Init constants
   N_GLI_PWA_               = N_GLI_PWA;
   p_grid_                  = p_grid;
   w_grid_                  = w_grid;
   mom_grid_size_           = mom_grid_size;
   J_max_                   = J_max;
   cutoff_Lambda_           = cutoff_Lambda;
   cut_pow_                 = cut_pow;   
   sharp_cutoff_            = sharp_cutoff;
   inc_grid_weights_in_pot_ = inc_grid_weights_in_pot;
   cut_on_shell_            = cut_on_shell;
   loop_reg_                = loop_reg;
   lam_SFR_                 = lam_SFR;

   // Construct terms and append them to terms_in_pot
   for (std::size_t i = 0; i < terms.size(); i++)
   {
      terms_in_pot_.push_back(Term(terms[i]));
      #ifdef ENABLE_DEBUG
      std::cout << "Added " << terms[i] << " to terms_in_pot_" << std::endl;
      #endif
   }
   
   // Make list of all LECs and params and then fetch the ones in use
   // -------------------------------------------------------------------------
   LEC_names_.push_back("C1S0");
   LEC_names_.push_back("C3S1");
   LEC_names_.push_back("D1S0");
   LEC_names_.push_back("D3S1");
   LEC_names_.push_back("D1P1");
   LEC_names_.push_back("D3P0");
   LEC_names_.push_back("D3P1");
   LEC_names_.push_back("D3P2");
   LEC_names_.push_back("D_SD");
   LEC_names_.push_back("D_DS");
   LEC_names_.push_back("E1S0");
   LEC_names_.push_back("E3P0");
   LEC_names_.push_back("E3P2");
   LEC_names_.push_back("E_PF");
   LEC_names_.push_back("E_FP");
   LEC_names_.push_back("F1S0");
   LEC_names_.push_back("c1");
   LEC_names_.push_back("c3");
   LEC_names_.push_back("c4");


   LEC_names_.push_back("Yamaguchi_1S0");
   LEC_names_.push_back("Yamaguchi_lambda_00");
   LEC_names_.push_back("Yamaguchi_lambda_01");
   LEC_names_.push_back("Yamaguchi_lambda_10");
   LEC_names_.push_back("Yamaguchi_lambda_11");
   
   param_names_.push_back("gA");
   param_names_.push_back("d18");
   param_names_.push_back("Yamaguchi_beta");

   // Initialize all LECs to zero
   for (std::size_t i = 0; i < LEC_names_.size(); i++)
   {
      LECs_.insert( std::make_pair (LEC_names_[i],0.0) );
   }
   
   // Initialize all params to zero
   for (std::size_t i = 0; i < param_names_.size(); i++)
   {
      params_.insert( std::make_pair (param_names_[i],0.0) );
   }

   // Update list of LECs and params that are in this instancs of the potential 
   // Use a set to not have duplicates
   std::set<std::string> tmp_set_LECs;
   std::set<std::string> tmp_set_params;
   for (std::size_t j = 0; j < terms_in_pot_.size(); j++)
   {
      std::vector<std::string> term_string_LECs = terms_in_pot_[j].get_lecs_in_term();
      for (std::size_t i = 0; i < term_string_LECs.size(); i++)
      {
         tmp_set_LECs.insert(term_string_LECs[i]);
      }

      std::vector<std::string> term_string_params = terms_in_pot_[j].get_params_in_term();
      for (std::size_t i = 0; i < term_string_params.size(); i++)
      {
         tmp_set_params.insert(term_string_params[i]);
      }
   }
   // Set the the vector equal to the set
   LECs_in_use_.assign(tmp_set_LECs.begin(),tmp_set_LECs.end());
   params_in_use_.assign(tmp_set_params.begin(),tmp_set_params.end());

   #ifdef ENABLE_DEBUG
      std::cout << "LECs in potential: ";
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         std::cout << LECs_in_use_[i] << ", ";
      }
      std::cout << std::endl;
   #endif

   // -------------------------------------------------------------------------
   
   // Make grid for angular integration
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   int_ang_ = gsl_integration_fixed_alloc(T, N_GLI_PWA_, -1.0, 1.0, 0, 0);

   z_mesh   = gsl_integration_fixed_nodes(int_ang_);
   w_z_mesh = gsl_integration_fixed_weights(int_ang_);
  
   len_z_mesh = N_GLI_PWA_;

   // Store some Legendre Polynomials for J = 0,...,J_max

   stored_Legendre_polynomials_ = (double**) malloc( (J_max+1) * sizeof(double*));

   for (int i = 0; i < (int)J_max_+1; i++)
   {
      stored_Legendre_polynomials_[i] = (double*) malloc(len_z_mesh*sizeof(double));
      for (int j = 0; j < (int)len_z_mesh; j++)
      {
         stored_Legendre_polynomials_[i][j] = gsl_sf_legendre_Pl(i, z_mesh[j]);
      }
   }

   #ifdef ENABLE_DEBUG
      std::cout << "V object created" << std::endl;
   #endif
}

// Destructor
template <class gsl_m>
Potential_mwpc<gsl_m>::~Potential_mwpc()
{
   // Free all memory allocations
   gsl_integration_fixed_free(int_ang_);
  
   for (int i = 0; i < (int)J_max_; i++)
   {
      free(stored_Legendre_polynomials_[i]);
   }
   free(stored_Legendre_polynomials_);

   // Delete all the stored matrices    
    if (!saved_matrices_.empty())
    {
        for (auto it = saved_matrices_.begin(); it != saved_matrices_.end(); it++)
        {
            for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
            {
                ph::matrix_free(it->second[LECs_in_use_[i]]);
            }
        }
    }
    #ifdef ENABLE_DEBUG
        std::cout << "V object deleted" << std::endl;
    #endif
}


/*
   Computes the total isospin factor from \tau_1 \cdot \tau_2 from
   the constraint J+L+T = odd 
*/
template <class gsl_m>
int Potential_mwpc<gsl_m>::isoFac(int L,int S)
{
   int T = ((1-L-S) & 1); // L+S+T=odd
   int res = -3*(1-T) + T;
   //std::cout << "isoFac: " << res << std::endl;
   return res; // Return factor from \tau_i \cdot \tau_2 in terms of T
}


template <class gsl_m>
double Potential_mwpc<gsl_m>::calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz)
{
   // TODO implement
   return 0;
}

template <class gsl_m>
double Potential_mwpc<gsl_m>::compute_A_integral(double qi, double qo, int J,int l, 
        std::vector<double> v_alpha_arr)
{
   #ifdef ENABLE_DEBUG
   if (!(J<J_max_)) 
   {
      std::cerr << "Error in compute_A_integral(): J>=J_max, returning 0" << std::endl;
      return 0;
   }
   #endif
   double integral = 0;
   for (int i = 0; i < (int)len_z_mesh; i++)
   {
      integral += w_z_mesh[i] * v_alpha_arr[i] * 
         gsl_pow_int(z_mesh[i],l) * stored_Legendre_polynomials_[J][i];
      
      //integral += w_z_mesh[i] * pot_OPEP_mom(qi,qo,z_mesh[i]) * 
      //   gsl_pow_int(z_mesh[i],l) * gsl_sf_legendre_Pl(J, z_mesh[i]);
   }
   return integral*M_PI;
}

/*
   This function returns in the V_arr pointer, an array of lengt 6 which include the potential 
   elemetns V_arr = [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp] where 
   S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc.

   The potential elements are in the same normalization as in the input terms that are entered in 
   the class Term. This is a relativistic normalization (see README.md). This means that to, for example,
   construct a motential matrix in a momentum basis <p'|p> = \delta^3(p'-p) a relativistic factor
   needs to be added. The  (2\pi)^{-3} factor ARE included.
*/

template <class gsl_m>
void Potential_mwpc<gsl_m>::calc_element_V_arr(double qi,double qo, 
        qs::quantum_channel chn, double* V_arr)
{
    bool coupled = chn.coupled;
    int J = chn.J;
   // Define some variables and initialize V_arr to all zeros.
   double V_uncoupled_S0 = 0;
   double V_uncoupled_S1 = 0;
   double V_coupled_mm   = 0;
   double V_coupled_pm   = 0;
   double V_coupled_mp   = 0;
   double V_coupled_pp   = 0;
   for (std::size_t i = 0; i < terms_in_pot_.size(); i++)
   {  
      //std::cout << terms_in_pot_[i].get_term_name() << std::endl;
      if (!terms_in_pot_[i].well_def_pw())
      {
         #ifdef ENABLE_DEBUG
            std::clock_t start,end;
            start = std::clock();
         #endif
         // Compute v_alpha array. Just make this function call ONCE!
         //std::cout << "Term: " << terms_in_pot_[i].get_term_name() << std::endl;
         std::vector<double> v_alpha_arr = 
             terms_in_pot_[i].my_v_alpha(qi,qo,z_mesh,len_z_mesh,LECs_,params_,
                     chn,loop_reg_,lam_SFR_);

         // Call pwa with the correct spin structure from this term
         // This will fill up the array V_arr with the correct potential elements
         
         double tmp_arr[6];
         //std::cout << terms_in_pot_[i].get_spin_structure() << " " << 
         //    terms_in_pot_[i].get_isovector() << std::endl;
         pwa(qi,qo,coupled,J,terms_in_pot_[i].get_spin_structure(),
                 terms_in_pot_[i].get_isovector(),v_alpha_arr,&tmp_arr[0]);
         
         V_uncoupled_S0 += tmp_arr[0];
         V_uncoupled_S1 += tmp_arr[1];
         V_coupled_pp   += tmp_arr[2];
         V_coupled_mm   += tmp_arr[3];
         V_coupled_pm   += tmp_arr[4];
         V_coupled_mp   += tmp_arr[5];
         //std::cout << V_uncoupled_S0 << " " << V_uncoupled_S1 << std::endl;
      
         #ifdef ENABLE_DEBUG
            end = std::clock();
            std::cout << "Time taken to call pwa is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
            std::cout << " mu sec " << std::endl;
         #endif

      } else 
      {
         // Get JLS from the Term to know where to apend the lec value 
         LS_term LS_term = terms_in_pot_[i].get_LS_term();
         //std::cout << "LS-term: " << LS_term.J  << " " << LS_term.Li << " " << LS_term.Lo << " " <<
         //    LS_term.S << std::endl;
         //std::cout << "Initial: " << V_coupled_mm << " " << V_coupled_pp << " " << V_coupled_mp <<
         //   V_coupled_pm << " " << V_uncoupled_S0 << " " << V_uncoupled_S1 << std::endl; 
         // Find which ouput to append (if any)
         if (J==LS_term.J)
         {
            if (coupled)
            {
               if (LS_term.Li == J-1 && LS_term.Lo == J-1) // --
               {
                  V_coupled_mm += 
                      terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
               } else if (LS_term.Li == J+1 && LS_term.Lo == J+1) // ++
               {
                  V_coupled_pp += 
                      terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // -+
               {
                  V_coupled_mp += 
                      terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
               } else if (LS_term.Li == J+1 && LS_term.Lo == J-1) // +-
               {
                  V_coupled_pm += 
                      terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
               }
            } else
            {
               if (LS_term.Li == J) // Li=Lo=J for uncoupled
               {
                  if (LS_term.S == 0) // S0
                  {
                     V_uncoupled_S0 += 
                         terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
                     
                  } else if (LS_term.S == 1) // S1
                  {
                     V_uncoupled_S1 += 
                         terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
                  }
               } else if (LS_term.Li == 1 && LS_term.Lo == 1 && LS_term.J == 0 
                       && LS_term.S == 1) // 3P0 exception
               {
                   V_coupled_pp += 
                       terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_,params_);
               }
            } 

            //std::cout << "Final: " << V_coupled_mm << " " << V_coupled_pp << " " << V_coupled_mp <<
            //    V_coupled_pm << " " << V_uncoupled_S0 << " " << V_uncoupled_S1 << std::endl; 
         }
      } 
   }
   double fac = gsl_pow_int(2*M_PI,-3);
    
   /*
    * The minus sign on the last two potetnial terms is from the fact that 
    * we use the Machleidt convention. This comes from the fact that 
    * relation betweent the M-matrix and the partial wave amplitudes.
    * If you choose to have them positive, the off diagonal T-matrix
    * elements change sign, and hence the mixing angle. This means that 
    * to keep the observables unchanges the factor i^(l'-l) must be excluded
    * from the M-T_matrix relation.
    */
   V_arr[0] = V_uncoupled_S0*fac;
   V_arr[1] = V_uncoupled_S1*fac;
   V_arr[2] = V_coupled_pp*fac;
   V_arr[3] = V_coupled_mm*fac;
   V_arr[4] = -V_coupled_pm*fac;
   V_arr[5] = -V_coupled_mp*fac;
}
/*
   This funtion performs the partial wave decomposition, which depends on the tensor structure of the 
   respective term in the potential.
*/
template <class gsl_m>
void Potential_mwpc<gsl_m>::pwa(double qi,double qo, bool coupled, int J_int, 
        std::string spin_struct,
        bool isovector, std::vector<double>& v_alpha_arr,double* V_arr)
{
    //std::cout << "spin structure: " << spin_struct << std::endl;
    /*
     * Integral names:
     * ---------------
     * A_0 -> A^(J),(0)
     * A_1 -> A^(J),(1)
     * A_2 -> A^(J),(2)
     * A_M -> A^(J-1),(0)
     * A_P -> A^(J+1),(0)
     * A_M2 -> A^(J-2),(0)
     * A_P2 -> A^(J+2),(0)
     * A_M_1 -> A^(J-1),(1)
     * A_P_1 -> A^(J+1),(1)
     * A_M_2 -> A^(J-1),(2)
     * A_P_2 -> A^(J+1),(2)
     */

    /*
    double A_0 = compute_A_integral(qi,qo,J,0,v_alpha_arr);
    double A_1 = compute_A_integral(qi,qo,J,1,v_alpha_arr);
    double A_2 = compute_A_integral(qi,qo,J,2,v_alpha_arr);
    
    double A_M = compute_A_integral(qi,qo,J-1,0,v_alpha_arr);
    double A_P = compute_A_integral(qi,qo,J+1,0,v_alpha_arr);
    
    double A_M2 = compute_A_integral(qi,qo,J-2,0,v_alpha_arr);
    double A_P2 = compute_A_integral(qi,qo,J+2,0,v_alpha_arr);
    
    double A_M_1 = compute_A_integral(qi,qo,J-1,1,v_alpha_arr);
    double A_P_1 = compute_A_integral(qi,qo,J+1,1,v_alpha_arr);
    
    double A_M_2 = compute_A_integral(qi,qo,J-1,2,v_alpha_arr);
    double A_P_2 = compute_A_integral(qi,qo,J+1,2,v_alpha_arr);
    */

    // Case J to a double to avoid integrer errors
    double J = J_int;

    // Define some variables
    double V_uncoupled_S0 = 0;
    double V_uncoupled_S1 = 0;
    double V_coupled_mm   = 0;
    double V_coupled_pm   = 0;
    double V_coupled_mp   = 0;
    double V_coupled_pp   = 0;


    if (spin_struct == "C")
    {
        if (!coupled)
        {
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);

            V_uncoupled_S0 = 2.0*A_0;
            V_uncoupled_S1 = 2.0*A_0;
        } 
        
        if (coupled || J_int == 0)
        {
            double A_P = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            V_coupled_pp = 2.0*A_P;
            if (J_int > 0)
            {
                double A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);
                V_coupled_mm = 2.0*A_M;
                V_coupled_mp = 0.0;
                V_coupled_pm = 0.0;
            }
        }
    } 
    else if (spin_struct == "S")
    {
        if (!coupled)
        {
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);
            V_uncoupled_S0 = -6.0*A_0;
            V_uncoupled_S1 = 2.0*A_0;
        } 
        
        if (coupled || J_int == 0)
        {
            double A_P = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            V_coupled_pp = 2.0*A_P;
            if (J_int>0)
            {
                double A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);
                V_coupled_mm = 2.0*A_M;
                V_coupled_mp = 0.0;
                V_coupled_pm = 0.0;
            }
        }
    }
    else if (spin_struct == "LS")
    {
        if (!coupled)
        {
            double A_P = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            double A_M = 0;
            if (J_int>0)
            {
                A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);
            }
            V_uncoupled_S0 = 0;
            V_uncoupled_S1 = 2.0*qo*qi*(1.0/(2*J+1))*(A_P - A_M);
        } 
        
        if (coupled || J_int == 0)
        {
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);
            double A_P2 = compute_A_integral(qi,qo,J_int+2,0,v_alpha_arr);
            
            V_coupled_pp = 2.0*qo*qi*((J+2)/(2*J+3))*(A_P2 - A_0);
            if (J_int>0)
            {
                double A_M2 = 0;
                if (J_int>1)
                {
                    A_M2 = compute_A_integral(qi,qo,J_int-2,0,v_alpha_arr);
                }  
                V_coupled_mm = 2.0*qo*qi*((J-1)/(2*J-1))*(A_M2 - A_0);
                V_coupled_mp = 0.0;
                V_coupled_pm = 0.0;
            }
        }
    } 
    else if (spin_struct == "T")
    {
        // These potential elements are computed only if it in an uncoupled
        // channel
        if (!coupled)
        {
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);
            double A_1 = compute_A_integral(qi,qo,J_int,1,v_alpha_arr);
            double A_P = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            double A_M = 0;
            if (J_int>0)
            {
                A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);
            }

            V_uncoupled_S0 = 2.0 * (-(qo*qo+qi*qi)*A_0 + 2.0*qo*qi*A_1);
            V_uncoupled_S1 = 2.0 * ((qo*qo+qi*qi)*A_0 - 2.0*qo*qi*(1.0/(2.0*J+1.0))
                    *((double)J*A_P + (double)(J+1.0)*A_M));
        }
        
        // These matrix elements are computed if the channel is uncoupled OR
        // J=0, since for J=0 the 3P0 channel is uncoupled and come from the 
        // V_coupled_mm matrix element.
        if (coupled || J_int == 0)
        {
            double A_P = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);
            
            V_coupled_pp = (2.0/(2.0*J+1.0)) * (-(qo*qo+qi*qi)*A_P + 2.0*qo*qi*A_0);
            if (J_int>0)
            {
                double A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);

                V_coupled_mm = (2.0/(2.0*J+1.0)) * ((qo*qo+qi*qi)*A_M - 2*qo*qi*A_0);
                V_coupled_mp = (4.0*sqrt(J*(J+1.0))/(2.0*J+1.0)) * 
                    (qi*qi*A_P + qo*qo*A_M - 2*qo*qi*A_0);
                V_coupled_pm = (4.0*sqrt(J*(J+1))/(2.0*J+1.0)) * 
                    (qi*qi*A_M + qo*qo*A_P - 2*qo*qi*A_0);
            }
        }
    }
    else if (spin_struct == "sigmaL")
    {
        if (!coupled)
        {
            double A_0 = compute_A_integral(qi,qo,J_int,0,v_alpha_arr);
            double A_2 = compute_A_integral(qi,qo,J_int,2,v_alpha_arr);
            double A_M_1 = 0;
            if (J_int>0)
            {
                A_M_1 = compute_A_integral(qi,qo,J_int-1,1,v_alpha_arr);
            }
            double A_P_1 = compute_A_integral(qi,qo,J_int+1,1,v_alpha_arr);
            
            V_uncoupled_S0 = 2.0*qo*qo*qi*qi*(A_2-A_0);
            V_uncoupled_S1 = 2.0*qo*qo*qi*qi*(-A_0 + ((J-1.0)/(2.0*J+1.0))*
                    A_P_1 + ((J+2.0)/(2.0*J+1.0))*A_M_1);
        } 
        
        if (coupled || J_int == 0)
        {
            double A_1   = compute_A_integral(qi,qo,J_int,1,v_alpha_arr);
            double A_P   = compute_A_integral(qi,qo,J_int+1,0,v_alpha_arr);
            double A_P_2 = compute_A_integral(qi,qo,J_int+1,2,v_alpha_arr);

            V_coupled_pp = 2.0*qo*qo*qi*qi*( ((2.0*J+3)/(2.0*J+1))*A_P -
                        ((2.0)/(2.0*J+1))*A_1 - A_P_2);
            if (J_int>0)
            {
                double A_M = compute_A_integral(qi,qo,J_int-1,0,v_alpha_arr);
                double A_M_2 = compute_A_integral(qi,qo,J_int-1,2,v_alpha_arr);
                
                V_coupled_mm = 2.0*qo*qo*qi*qi*( ((2.0*J-1)/(2.0*J+1))*A_M +
                        ((2.0)/(2.0*J+1))*A_1 - A_M_2);
                V_coupled_mp = 4.0*qo*qo*qi*qi*(std::sqrt(J*(J+1.0))/
                        ((2.0*J+1.0)*(2.0*J+1.0)))*(A_P - A_M);
                V_coupled_pm = V_coupled_mp;
            }
        }
    }
    else if (spin_struct == "sigmak")
    {
        // TODO
    } 
    else
    {
        std::cout << "Unknown <spin_struct> passed to function <pwa> in Potential_mwpc object." 
            << std::endl;
    }
    
    // Add isospin factor if the term is an isovector
    if (isovector)
    {
        //std::cout << "isovector" << std::endl;
        V_uncoupled_S0 *= isoFac(J_int,0);
        V_uncoupled_S1 *= isoFac(J_int,1);
        V_coupled_mm   *= isoFac(J_int-1,1);
        V_coupled_pm   *= isoFac(J_int+1,1);
        V_coupled_mp   *= isoFac(J_int-1,1);
        V_coupled_pp   *= isoFac(J_int+1,1);
    }

    // Populate the array
    V_arr[0] = V_uncoupled_S0;
    V_arr[1] = V_uncoupled_S1;
    V_arr[2] = V_coupled_pp;
    V_arr[3] = V_coupled_mm;
    V_arr[4] = V_coupled_pm;
    V_arr[5] = V_coupled_mp;
}

/*
   This function is the key function to compute the potential matrix unsing the potential matrix elements
   conputed from calc_element_V_arr(). calc_element_V_arr() returns potential elements in a relativistic 
   momentum basis (see REDME.md) in which the LS equation has a factor of $(2\pi)^{-3} as well as
   a relativistic correction factor. To get the potential in the correclty normalized non-relativistic
   three-momentum basis <p'|p> = \delta(p'- p) these factors are included in the potential.

   Input:
      channel: (coupled, J,S)
      regulator

   Output:
      Potential matrix [V] including relativistic corrections (minimal relativity). The basis is a partial wave basis 
      normalized as (see README.md) springing from the three-momentum basis normalized as <p'|p> = \delta(p'- p).
      The normalization here affects the factors in the partial wave projected LS-equation as well as 
      the relation of R to phase shifts. (see README.md)
*/

template <class gsl_m>
gsl_m* Potential_mwpc<gsl_m>::get_matrix(double q_on_shell,qs::quantum_channel chn,bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "get_matrix()" << std::endl;
   #endif
   double mu = ph::get_mN(chn.Tz)/2.0; // Default

   // Allocate gsl matrices in the case of coupled and uncoupled channels.
   // The matrix becomes twise as large in the coupled case
   // In the construction of the saved matrix the desired on-shell momentum 
   // that the matrix element will be evaluated on are unknown. These matrix elements needs
   // to be computed at runtime
   gsl_m* matrix_data = nullptr;
   if (chn.coupled) {
      matrix_data = ph::matrix_alloc((2*mom_grid_size_ + 2),(2*mom_grid_size_ + 2),
              matrix_data);
   } else {
      matrix_data = ph::matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1,
              matrix_data);
   }

   // Set all elements to zero to be able to add the contribution of the on-shell 
   // part later. This ensures that the elements in these saved matrices are zero.
   ph::matrix_set_zero(matrix_data);

   // For each grid point, including qi, and qo
   
   // i: row index, j: column index
   // These loops populate the matrices everywhere except where
   // the on-shell part will go later.
   for (std::size_t i = 0; i < mom_grid_size_+1; i++)
   {
      for (std::size_t j = 0; j < mom_grid_size_+1; j++)
      {
         double V_arr[6]; // Array for data
         
         // Outgoing momentum is row index
         double p_in  = 0;
         double p_out = 0;

         if (j < mom_grid_size_) {
            p_in  = p_grid_[j];
         } else {
            p_in = q_on_shell;
         }
         if (i < mom_grid_size_) {
            p_out = p_grid_[i];
         } else {
            p_out = q_on_shell;
         }

         // Get the factor from the relativistic corrections
         // the cutoff and the grid
         double tot_fac = get_total_rel_cut_weight_factor(p_in,j,p_out,i,mu,rel_correction);

         //std::cout << " LECS: " << LECs_["gA2"] << " " << LECs_["C1S0"] << " " << LECs_["C3S1"] << std::endl;
         calc_element_V_arr(p_in,p_out,chn,&V_arr[0]);
         //std::cout << "Rel fac: " << rel_fac << std::endl;
         /*for (int i= 0; i < 6; i++)
         {
            std::cout << V_arr[i] << " ";
         }*/
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
                ph::matrix_set(matrix_data,j,i,V_arr[0]*tot_fac);
               //std::cout << "Pot el S0: " << V_arr[0]*rel_fac << std::endl;
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               if (chn.J != 0)
               {
                   ph::matrix_set(matrix_data,j,i,V_arr[1]*tot_fac);
               } else // For J=0,S=1,L=1 case
               {
                   ph::matrix_set(matrix_data,j,i,V_arr[2]*tot_fac); // Take pp element to get L=1
               }
               
            }
         } else 
         {
             //std::cout << "tot_fac=" << tot_fac << std::endl;
            // The matrix is constructed as [[mm,mp],[pm,pp]]
             ph::matrix_set(matrix_data,j,i,V_arr[3]*tot_fac); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
             ph::matrix_set(matrix_data,j,i+(mom_grid_size_+1),V_arr[5]*tot_fac); //mp
             ph::matrix_set(matrix_data,j+(mom_grid_size_+1),i,V_arr[4]*tot_fac); //pm
             ph::matrix_set(matrix_data,j+(mom_grid_size_+1),i+(mom_grid_size_+1),V_arr[2]*tot_fac); //pp

         }
      }
   }
   //std::cout << "Potential matrix from inside the potential:" << std::endl;
   //ph::print_m(matrix_data);
   // Return the matrix
   return matrix_data;
}

template <class gsl_m>
void Potential_mwpc<gsl_m>::clear_saved_matrices()
{
   #ifdef ENABLE_DEBUG
      std::cerr << "clear_saved_matrices()" << std::endl;
   #endif
   // Go through the map and remove the created matrices
   for (typename std::map<qs::quantum_channel, std::map<std::string, gsl_m*>, qs::comp>::iterator 
           it1=saved_matrices_.begin(); it1!=saved_matrices_.begin(); ++it1)
   {
      for (typename std::map<std::string, gsl_m*>::iterator it=it1->second.begin(); 
              it!=it1->second.end(); ++it)
      {
          ph::matrix_free(it->second);
      }
   }
   #ifdef ENABLE_DEBUG
      std::cerr << "clear_saved_matrices() - end" << std::endl;
   #endif
   
}

/*
   This function popolates the saved matrices in the given channel.
*/
template <class gsl_m>
void Potential_mwpc<gsl_m>::populate_saved_mtx(qs::quantum_channel chn, bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "populate_saved_mtx()" << std::endl;
   #endif
    
   clear_saved_matrices(); // Clears the allocated pointers

   // Get the constant
   for (std::size_t j = 0; j < LEC_names_.size(); j++)
   {
      LECs_[LEC_names_[j]] = 0.0;
   }
   gsl_m* mtx_const = get_matrix(0.0,chn,rel_correction); // Do not care about on-shell part
   // Save it
   saved_matrices_[chn]["const"] = mtx_const;

   // Set each lec=1 and the rest to 0
   for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
   {
      // Reset ALL lecs to zero, just to be sure
      // Would suffice to do it for all LECs in use
   
      for (std::size_t j = 0; j < LEC_names_.size(); j++)
      {
         LECs_[LEC_names_[j]] = 0.0;
      }
   
      // Set one LEC to 1.0;
      LECs_[LECs_in_use_[i]] = 1.0;
         
      // Compute potential matrix for these LECs
      // Note that the function get_matrix() will compute the matrix elements using the 
      // LECs in the private variable LECs_ 
      gsl_m* matrix = get_matrix(0.0,chn,rel_correction); // Do not care about on-shell part
      
      // Subtract the constant matrix
      ph::matrix_sub(matrix,mtx_const);

      #ifdef ENABLE_DEBUG
         std::cerr << std::endl << "populate_saved_mtx(): Printing matrix associated to LEC: " << LECs_in_use_[i] << "." << std::endl;
         print_m(matrix);
      #endif
      // Save the matrix as the matrix corresponding tho the non-zero LEC.
      saved_matrices_[chn][LECs_in_use_[i]] = matrix;
   }
   LECs_[LECs_in_use_[LECs_in_use_.size()-1]] = 0.0;
   // After the function has populated the matrices the default is that all LECs are set to zero, 
   // this is to alert the user if the LECs are not updater before getting the matrix.
}

/*
   This function returns the potential matrix in the given channel. It does not recompute all matrix elements
   but loades the saved part and just fills in the parts that depends on the on-shell momentum.
*/


template <class gsl_m>
gsl_m* Potential_mwpc<gsl_m>::get_saved_matrix(double q_on_shell, qs::quantum_channel chn, bool rel_correction)
{
   unsigned int J = chn.J;
   unsigned int S = chn.S;
   bool coupled = chn.coupled;

   double mu = ph::get_mN(chn.Tz)/2.0;
   
   // Note that the matrix is computed for the current LECs_!!
   // Note also that the LECs_ are screwed up by the act of saving the matrices!
   // This means that the 
   
   // Check if saved otherwise print error message
   if (saved_matrices_[chn].empty())
   {
      std::cerr << "get_saved_matrix(): There is no saved matrices in this channel, returning." << std::endl;
      return nullptr;
   }

   // Load the matrices and compute the necessary matrix elements on the diagonal
   
   gsl_m* matrix_saved_sum = nullptr;
   if (!coupled)
   {
      matrix_saved_sum  = ph::matrix_alloc((size_t)(mom_grid_size_ + 1),
              (size_t)(mom_grid_size_ + 1),matrix_saved_sum);
      gsl_m* tmp_matrix = ph::matrix_alloc((size_t)(mom_grid_size_ + 1),
              (size_t)(mom_grid_size_ + 1),matrix_saved_sum);
      ph::matrix_set_zero(matrix_saved_sum);
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
          ph::matrix_memcpy(tmp_matrix, saved_matrices_[chn][LECs_in_use_[i]]);
         //p_m(saved_matrices_[chn][LECs_in_use_[i]]);
         // Scale tmp_matrix by the correct LEC
          ph::matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);
         //std::cout << "LEC=" << LECs_[LECs_in_use_[i]] << std::endl;
         // Add this contribution to the sum

          ph::matrix_add(matrix_saved_sum,tmp_matrix);
      }
      // Add constant
      ph::matrix_add(matrix_saved_sum,saved_matrices_[chn]["const"]);
      // Delete the tmp matrix
      ph::matrix_free(tmp_matrix);

      /* The dashed part are filled in
         (     |)
         (     |)
         (------)
      */
      // Compute the neccesary elements, it suffice to just calculate the row
      for (std::size_t i = 0; i < mom_grid_size_+1; i++)
      {
         double tmp_arr[6];
         double p_in;
         if (i < mom_grid_size_) {
            p_in = p_grid_[i];
         } else {
            p_in = q_on_shell;
         }
         double p_out = q_on_shell; // row is fixed
         calc_element_V_arr(p_in,p_out, chn, &tmp_arr[0]);


         // Get the factor from the relativistic corrections
         // the cutoff and the grid
         double tot_fac = get_total_rel_cut_weight_factor(p_in,i,p_out,mom_grid_size_,mu,rel_correction);
         
         if (J == 0 && S == 1) // 3P_0 case
         {
             ph::matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[2]*tot_fac);
             ph::matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[2]*tot_fac);
         }
         else if (S==0)
         {
            // Take element zero and insert it into the matrix
             ph::matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[0]*tot_fac);
             ph::matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[0]*tot_fac);
            
         } else 
         {
             ph::matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[1]*tot_fac);
             ph::matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[1]*tot_fac);
         }
      }

   } else 
   {
      matrix_saved_sum       = ph::matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2,matrix_saved_sum);
      gsl_m* tmp_matrix = ph::matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2,matrix_saved_sum);
      ph::matrix_set_zero(matrix_saved_sum);
       for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
          ph::matrix_memcpy(tmp_matrix, saved_matrices_[chn][LECs_in_use_[i]]);

         // Scale tmp_matrix by the correct LEC
         //std::cout << "LEC=" << LECs_[LECs_in_use_[i]] << std::endl;
          ph::matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);

         // Add this contribution to the sum
          ph::matrix_add(matrix_saved_sum,tmp_matrix);
      }
      // Add constant
       ph::matrix_add(matrix_saved_sum,saved_matrices_[chn]["const"]);
     
      // Delete the tmp matrix
       ph::matrix_free(tmp_matrix);

      /* The dashed parts are filled in
         (     |)(     |)
         (     |)(     |)
         (------)(------)
         (     |)(     |)
         (     |)(     |)
         (------)(------)
      */
      // Compute the neccesary elements and add them
      for (std::size_t i = 0; i < mom_grid_size_+1; i++)
      {
         double tmp_arr[6];
         double p_in;
         if (i < mom_grid_size_) {
            p_in = p_grid_[i];
         } else {
            p_in = q_on_shell;
         }
         // column index
         double p_out = q_on_shell; // row is fixed
         calc_element_V_arr(p_in,p_out, chn, &tmp_arr[0]);

         // Get the factor from the relativistic corrections
         // the cutoff and the grid
         double tot_fac = get_total_rel_cut_weight_factor(p_in,i,p_out,mom_grid_size_,mu,rel_correction);
         
         
         // In this part of the code the symmetry of the potential in momentum is taken
         // advantage of. Note that this does NOT mean that every block is symmetric,
         // but rather the WHOLE matrix.

         // Take mm element and insert it into the matrix
         // mm
         ph::matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[3]*tot_fac); // Column
         ph::matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[3]*tot_fac); // Row

         // Take mp element one and insert it into the matrix
         // mp 
         ph::matrix_set(matrix_saved_sum,mom_grid_size_,i+mom_grid_size_+1,tmp_arr[4]*tot_fac); // Column
         ph::matrix_set(matrix_saved_sum,i,2*mom_grid_size_+1,tmp_arr[5]*tot_fac); // Row 

         // Take pm element one and insert it into the matrix
         // pm
         ph::matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i,tmp_arr[5]*tot_fac); // Column
         ph::matrix_set(matrix_saved_sum,i+mom_grid_size_+1,mom_grid_size_,tmp_arr[4]*tot_fac); // Row 

         // Take pp element one and insert it into the matrix
         // pp
         ph::matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i+mom_grid_size_+1,tmp_arr[2]*tot_fac); // Column
         ph::matrix_set(matrix_saved_sum,i+mom_grid_size_+1,2*mom_grid_size_+1,tmp_arr[2]*tot_fac); // Row 
      }
   }

   #ifdef ENABLE_DEBUG
      std::cout << std::endl << "get_saved_matrix(): Printing matrix." << std::endl;
      ph::print_m(matrix_saved_sum);
   #endif
     
   // Return pointer to the full potential matrix.

   //std::cout << gsl_matrix_get(matrix_saved_sum,3,3) << std::endl;
   return matrix_saved_sum;
}


template <class gsl_m>
gsl_m* Potential_mwpc<gsl_m>::get_matrix_no_onshell(qs::quantum_channel chn, 
        bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "get_matrix()" << std::endl;
   #endif
   double mu = ph::get_mN(chn.Tz)/2.0;

   // Allocate gsl matrices in the case of coupled and uncoupled channels.
   // The matrix becomes twise as large in the coupled case
   // In the construction of the saved matrix the desired on-shell momentum 
   // that the matrix element will be evaluated on are unknown. These matrix elements needs
   // to be computed at runtime
   gsl_m* matrix_data = nullptr;
   if (chn.coupled) {
      matrix_data = ph::matrix_alloc((2*mom_grid_size_),(2*mom_grid_size_),matrix_data);
   } else {
      matrix_data = ph::matrix_alloc(mom_grid_size_,mom_grid_size_,matrix_data);
   }

   // For each grid point
   
   // i: row index, j: column index
   // These loops populate the matrices everywhere except where
   // the on-shell part will go later.
   for (std::size_t i = 0; i < mom_grid_size_; i++)
   {
      for (std::size_t j = 0; j < mom_grid_size_; j++)
      {
         double V_arr[6]; // Array for data
         
         // Outgoing momentum is row index
         double p_in  = p_grid_[j];
         double p_out = p_grid_[i];
         
         // Get the factor from the relativistic corrections
         // the cutoff and the grid
         double tot_fac = get_total_rel_cut_weight_factor(p_in,j,p_out,i,mu,rel_correction);
         
         calc_element_V_arr(p_in, p_out, chn, &V_arr[0]);
     
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
                ph::matrix_set(matrix_data,j,i,V_arr[0]*tot_fac);
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               if (chn.J != 0)
               {
                   ph::matrix_set(matrix_data,j,i,V_arr[1]*tot_fac);
               } else // For J=0,S=1,L=1 case
               {
                   ph::matrix_set(matrix_data,j,i,V_arr[2]*tot_fac); // Take pp element to get L=1
               }
            }

         } else 
         {
            // The matrix is constructed as [[mm,mp],[pm,pp]]
             ph::matrix_set(matrix_data,j,i,V_arr[3]*tot_fac); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
             ph::matrix_set(matrix_data,j,i+(mom_grid_size_),V_arr[5]*tot_fac); //mp
             ph::matrix_set(matrix_data,j+(mom_grid_size_),i,V_arr[4]*tot_fac); //pm
             ph::matrix_set(matrix_data,j+(mom_grid_size_),i+(mom_grid_size_),V_arr[2]*tot_fac); //pp
         }
      }
   }
   return matrix_data;
}

template <class gsl_m>
double Potential_mwpc<gsl_m>::get_rel_cut(double p_in, 
        double p_out, double mu, bool rel_correction)
{
     // Compute relativistic factors
     double rel_fac = 1.0;
     if (rel_correction)
     {
        double E_rel_in = sqrt(4*mu*mu+p_in*p_in);
        double E_rel_out = sqrt(4*mu*mu+p_out*p_out);
        double rel_factor_in = sqrt(2*mu/E_rel_in);
        double rel_factor_out = sqrt(2*mu/E_rel_out);
        rel_fac = rel_factor_in*rel_factor_out;
     }

     // Compute the cutoff factor
     double cutoff_regulator = exp(-gsl_pow_uint(p_in/cutoff_Lambda_,cut_pow_))*
            exp(-gsl_pow_uint(p_out/cutoff_Lambda_,cut_pow_));
     //std::cout << cutoff_regulator << ", " << rel_fac << std::endl;
     return cutoff_regulator*rel_fac;
}

template <class gsl_m>
double Potential_mwpc<gsl_m>::get_total_rel_cut_weight_factor(double p_in, int j, 
        double p_out, int i, double mu, bool rel_correction)
{
     // Compute relativistic factors
     double rel_fac = 1.0;
     if (rel_correction)
     {
        double E_rel_in = sqrt(4*mu*mu+p_in*p_in);
        double E_rel_out = sqrt(4*mu*mu+p_out*p_out);
        double rel_factor_in = sqrt(2*mu/E_rel_in);
        double rel_factor_out = sqrt(2*mu/E_rel_out);
        rel_fac = rel_factor_in*rel_factor_out;
     }

     // Compute the cutoff factor
     double cutoff_regulator = 1.0;
     if (cut_on_shell_==true)
     {
        cutoff_regulator = exp(-gsl_pow_uint(p_in/cutoff_Lambda_,cut_pow_))*
            exp(-gsl_pow_uint(p_out/cutoff_Lambda_,cut_pow_));
     } else 
     {
         // If ingoing momenta is of shell
         if (j<(int)mom_grid_size_)
         {
            cutoff_regulator *= exp(-gsl_pow_uint(p_in/cutoff_Lambda_,cut_pow_));
            /*if (p_in > cutoff_Lambda_)
            {
                cutoff_regulator = 0;
            } else 
            {
                cutoff_regulator = 1;
            }*/
         }
        
         // If outgoing momenta is of shell
         if (i<(int)mom_grid_size_)
         {
            cutoff_regulator *= exp(-gsl_pow_uint(p_out/cutoff_Lambda_,cut_pow_));
            /*if (p_in > cutoff_Lambda_)
            {
                cutoff_regulator = 0;
            } else 
            {
                cutoff_regulator = 1;
            }*/
         }
     }
     
     if (sharp_cutoff_) {
        if (p_in > cutoff_Lambda_ + 300.0 || p_out > cutoff_Lambda_ + 300.0) {
           cutoff_regulator = 0.0;
        }
     }

     // Compute the fectors that come from including the weights and momenta
     double weights_momenta = 1.0;
     // Include only in the off shell and half of shell part
     if (inc_grid_weights_in_pot_)
     {
         // If ingoing momenta is of shell
         if (j<(int)mom_grid_size_)
         {
            weights_momenta *= std::sqrt(w_grid_[j])*p_grid_[j];
         }
        
         // If outgoing momenta is of shell
         if (i<(int)mom_grid_size_)
         {
            weights_momenta *= std::sqrt(w_grid_[i])*p_grid_[i];
         }
     }

     // Multiply everything together
     return rel_fac*cutoff_regulator*weights_momenta;
}

template <class gsl_m>
void Potential_mwpc<gsl_m>::print_LECs_and_params_info()
{
    std::cout << std::endl << 
        "#####################################################################"
        << std::endl;
    std::cout << 
        "#####################################################################"
        << std::endl;

    std::cout << "                     Printing potential info" << std::endl;

    std::cout << 
        "#####################################################################" 
        << std::endl;
    std::cout << 
        "#####################################################################"
        << std::endl;

    std::cout << "Terms in potential:" <<  std::endl;
    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;
    for (auto& T : terms_in_pot_)
    {
        std::cout << T.term_name_ << std::endl;
    }    

    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl << "The LECs in the same order as they must be set in the"
              << " set functions:" << std::endl;
    std::cout << "Current set LEC values are shown (All units in powers of MeV)" 
        << std::endl;
    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;

    for (auto s : LECs_in_use_)
    {
        std::cout << std::setw(20) << s << " = " << LECs_[s] << std::endl; 
    }
    
    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl << "The params in the same order as they must be set in the"
              << " set functions:" << std::endl;
    std::cout << "(All units in powers of MeV)" << std::endl;
    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;

    for (auto s : params_in_use_)
    {
        std::cout << std::setw(20) << s << " = " << params_[s] << std::endl; 
    }
    std::cout << 
        "---------------------------------------------------------------------"
        << std::endl;
    std::cout << 
        "#####################################################################"
        << std::endl;
    std::cout << 
        "#####################################################################" 
        << std::endl << std::endl;

}
