#include "pot_nn_mwpc.h"

// This struct is to define what set of parameters the function f_int want
struct my_f_params { double qi; double qo; int J; int l; Term* term; Potential_mwpc* this_pot;};


// Constructor
Potential_mwpc::Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,double* p_grid, double* w_grid, std::size_t mom_grid_size, unsigned int J_max)
{
   // Init constants
   N_GLI_PWA_ = N_GLI_PWA;
   p_grid_ = p_grid;
   w_grid_ = w_grid;
   mom_grid_size_ = mom_grid_size;
   J_max_ = J_max;
   
   // Construct terms and append them to terms_in_pot
   for (std::size_t i = 0; i < terms.size(); i++)
   {
      terms_in_pot_.push_back(Term(terms[i]));
      #ifdef ENABLE_DEBUG
      std::cout << "Added " << terms[i] << " to terms_in_pot_" << std::endl;
      #endif
   }
   
   // Make list of all LECs.
   LEC_names_.push_back("gA2");
   LEC_names_.push_back("C1S0");
   LEC_names_.push_back("C3S1");
   
   // Initialize all LECs to zero
   for (std::size_t i = 0; i < LEC_names_.size(); i++)
   {
      LECs_.insert( std::make_pair<std::string,double> (LEC_names_[i],0.0) );
   }

   // Update list of LECs that are in this instancs of the potential 
   // Use a set to not have duplicates
   std::set<std::string> tmp_set;
   for (std::size_t i = 0; i < terms_in_pot_.size(); i++)
   {
      std::vector<std::string> tmp_s = terms_in_pot_[i].get_lecs_in_term();
      for (std::size_t i = 0; i < tmp_s.size(); i++)
      {
         tmp_set.insert(tmp_s[i]);
      }
   }
   // Set the the vector equal to the set
   LECs_in_use_.assign(tmp_set.begin(),tmp_set.end());

   #ifdef ENABLE_DEBUG
      std::cout << "LECs in potential: ";
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         std::cout << LECs_in_use_[i] << ", ";
      }
      std::cout << std::endl;
   #endif

   // Make grid for angular integration
   const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
   int_ang_ = gsl_integration_fixed_alloc(T, N_GLI_PWA_, -1.0, 1.0, 0, 0);

   z_mesh   = gsl_integration_fixed_nodes(int_ang_);
   w_z_mesh = gsl_integration_fixed_weights(int_ang_);
  
   
   len_z_mesh = N_GLI_PWA_;

   // Store some Legendre Polynomials for J = 0,...,J_max

   stored_Legendre_polynomials_ = (double**) malloc( (J_max+1) * sizeof(double*));

   for (int i = 0; i < J_max_+1; i++)
   {
      stored_Legendre_polynomials_[i] = (double*) malloc(len_z_mesh*sizeof(double));
      for (int j = 0; j < len_z_mesh; j++)
      {
         stored_Legendre_polynomials_[i][j] = gsl_sf_legendre_Pl(i, z_mesh[j]);
      }
   }

   #ifdef ENABLE_DEBUG
      std::cout << "V object created" << std::endl;
   #endif
}

// Destructor
Potential_mwpc::~Potential_mwpc()
{
   // Free all memory allocations
   gsl_integration_fixed_free(int_ang_);
  
   for (int i = 0; i < J_max_; i++)
   {
      free(stored_Legendre_polynomials_[i]);
   }
   free(stored_Legendre_polynomials_);

   #ifdef ENABLE_DEBUG
      std::cout << "V object deleted" << std::endl;
   #endif
}


/*
   Computes the total isospin factor from \tau_1 \cdot \tau_2 from
   the constraint J+L+T = odd 
*/
int isoFac(int L,int S)
{
   int T = (1-L-S & 1); // L+S+T=odd
   return -3*(1-T) + T; // Return factor from \tau_i \cdot \tau_2 in terms of T
}


double Potential_mwpc::calc_element_JLS(double qi,double qo, int J, int L, int S, int Tz)
{
   // TODO implement
   return 0;
}

double Potential_mwpc::compute_A_integral(double qi, double qo, int J,int l, std::vector<double> v_alpha_arr)
{
   #ifdef ENABLE_DEBUG
   if (!(J<J_max_)) 
   {
      std::cerr << "Error in compute_A_integral(): J>=J_max, returning 0" << std::endl;
      return 0;
   }
   #endif
   double integral = 0;
   for (int i = 0; i < len_z_mesh; i++)
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

void Potential_mwpc::calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr)
{
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
         std::vector<double> v_alpha_arr = terms_in_pot_[i].my_v_alpha(qi,qo,z_mesh,len_z_mesh,LECs_);

         #ifdef ENABLE_DEBUG
            end = std::clock();
            std::cout << "Time taken to call get_v_alpha is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
            std::cout << " mu sec " << std::endl;
         #endif
         // Calculate the A_x integrals (don't calculate all)
         #ifdef ENABLE_DEBUG
            start = std::clock();
         #endif

         double A_0 = compute_A_integral(qi,qo,J,0,v_alpha_arr);
         double A_1 = 0;
         double A_P = compute_A_integral(qi,qo,J+1,0,v_alpha_arr);
         double A_M = 0;

         if (J!=0) {
            A_M = compute_A_integral(qi,qo,J-1,0,v_alpha_arr);
         }
         if (!coupled) {
            A_1 = compute_A_integral(qi,qo,J,1,v_alpha_arr);
         }

         #ifdef ENABLE_DEBUG
            end = std::clock();
            std::cout << "Time taken to compute A integrals is : " << 1000000.0*(double)(end-start)/(double)CLOCKS_PER_SEC; 
            std::cout << " mu sec " << std::endl;
         #endif
         // Call pwa with the correct spin structure from this term
         // This will fill up the array V_arr with the correct potential elements
         
         #ifdef ENABLE_DEBUG
            start = std::clock();
         #endif
         double tmp_arr[6];
         //std::cout << terms_in_pot_[i].get_spin_structure() << " " << terms_in_pot_[i].get_isovector() << std::endl;
         pwa(qi,qo,coupled,J,A_M,A_P,A_0,A_1,terms_in_pot_[i].get_spin_structure(),terms_in_pot_[i].get_isovector(),&tmp_arr[0]);
      
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
         // Find which ouput to append (if any)
         if (J==LS_term.J)
         {
            if (coupled)
            {
               if (LS_term.Li == J-1 && LS_term.Lo == J-1) // --
               {
                  V_coupled_mm += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
               } else if (LS_term.Li == J+1 && LS_term.Lo == J+1) // ++
               {
                  V_coupled_pp += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // -+
               {
                  V_coupled_mp += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // +-
               {
                  V_coupled_pm += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
               }
            } else
            {
               if (LS_term.S == 0) // S0
               {
                  V_uncoupled_S0 += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
                  //std::cout << terms_in_pot_[i].get_LEC_element(qi,qo,LECs_) << std::endl;
               } else if (LS_term.S == 1) // S1
               {
                  V_uncoupled_S1 += terms_in_pot_[i].get_v_alpha_well_def_pw(qi,qo,LECs_);
                  //std::cout << terms_in_pot_[i].get_LEC_element(qi,qo,LECs_) << std::endl;
               }
            }
         }
      } 
   }
   double fac = gsl_pow_int(2*M_PI,-3);

   V_arr[0] = V_uncoupled_S0*fac;
   V_arr[1] = V_uncoupled_S1*fac;
   V_arr[2] = V_coupled_pp*fac;
   V_arr[3] = V_coupled_mm*fac;
   V_arr[4] = V_coupled_pm*fac;
   V_arr[5] = V_coupled_mp*fac;
}

/*
   This funtion performs the partial wave decomposition, which depends on the tensor structure of the 
   respective term in the potential.
*/
void Potential_mwpc::pwa(double qi,double qo, bool coupled, int J,double A_M,double A_P,double A_0,double A_1,std::string spin_struct,bool isovector,double* V_arr)
{
   // Define some variables
   double V_uncoupled_S0 = 0;
   double V_uncoupled_S1 = 0;
   double V_coupled_mm   = 0;
   double V_coupled_pm   = 0;
   double V_coupled_mp   = 0;
   double V_coupled_pp   = 0;
   //std::cout << "A_0=" << A_0 << " A_P=" << A_P << " A_M=" << A_M << " A_1=" << A_1 << std::endl;
   if (spin_struct == "tensor")
   {
      // Check which elements that are non-zero by checking if the Channel is
      // coupled or not
      if (!coupled)
      {
         // OPEP uncoupled interactions
         V_uncoupled_S0    = 2.0 * (-(qo*qo+qi*qi)*A_0 + 2.0*qo*qi*A_1);
         
         if (J!=0) {
            V_uncoupled_S1 = 2.0 * ((qo*qo+qi*qi)*A_0 - 2.0*qo*qi*(1.0/(2.0*J+1.0))*(J*A_P + (J+1)*A_M));
         }
      } else
      {
         V_coupled_pp = (2.0/(2.0*J+1.0)) * (-(qo*qo+qi*qi)*A_P + 2.0*qo*qi*A_0);
         if (J!= 0)
         {
            V_coupled_mm = (2.0/(2*J+1)) * ((qo*qo+qi*qi)*A_M - 2*qo*qi*A_0);
            V_coupled_mp = (4.0*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_P + qo*qo*A_M - 2*qo*qi*A_0);
            V_coupled_pm = (4.0*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_M + qo*qo*A_P - 2*qo*qi*A_0);
         }
      }

      // Add isospin factor if the term is an isovector
      if (isovector)
      {
         V_uncoupled_S0 *= isoFac(J,0);
         V_uncoupled_S1 *= isoFac(J,1);
         V_coupled_mm   *= isoFac(J-1,1);
         V_coupled_pm   *= isoFac(J+1,1);
         V_coupled_mp   *= isoFac(J-1,1);
         V_coupled_pp   *= isoFac(J+1,1);
      }     
      V_arr[0] = V_uncoupled_S0;
      V_arr[1] = V_uncoupled_S1;
      V_arr[2] = V_coupled_pp;
      V_arr[3] = V_coupled_mm;
      V_arr[4] = V_coupled_pm;
      V_arr[5] = V_coupled_mp;
   } else
   {
      std::cerr << "Unknown <spin_struct> passed to function <pwa> in Potential_mwpc object." << std::endl;
      V_arr[0] = 0;
      V_arr[1] = 0;
      V_arr[2] = 0;
      V_arr[3] = 0;
      V_arr[4] = 0;
      V_arr[5] = 0;
   }
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
gsl_matrix* Potential_mwpc::get_matrix(double q_on_shell,qs::quantum_channel chn,bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "get_matrix()" << std::endl;
   #endif
   double mu;
   if (chn.tz == -1)
   {
      mu = constants::Mn/2.0; // nn
   } else if (chn.tz == 0)
   {
      mu = constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np
   } else if (chn.tz == 1)
   {
      mu = constants::Mp/2.0; // pp
   } else 
   {
      #ifdef ENABLE_DEBUG
         std::cerr << "Error in solve_in_chn(): Unknown isospin" << std::endl;
      #endif
   }

   // Allocate gsl matrices in the case of coupled and uncoupled channels.
   // The matrix becomes twise as large in the coupled case
   // In the construction of the saved matrix the desired on-shell momentum 
   // that the matrix element will be evaluated on are unknown. These matrix elements needs
   // to be computed at runtime
   gsl_matrix* matrix_data;
   if (chn.coupled) {
      matrix_data = gsl_matrix_alloc((2*mom_grid_size_ + 2),(2*mom_grid_size_ + 2));
   } else {
      matrix_data = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
   }

   // Set all elements to zero to be able to add the contribution of the on-shell 
   // part later. This ensures that the elements in these saved matrices are zero.
   gsl_matrix_set_zero(matrix_data);

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
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/450,6))*exp(-gsl_pow_uint(p_out/450,6));
         rel_fac *= cutoff_regulator;
      
         //std::cout << " LECS: " << LECs_["gA2"] << " " << LECs_["C1S0"] << " " << LECs_["C3S1"] << std::endl;
         calc_element_V_arr(p_in,p_out,chn.coupled,chn.J,&V_arr[0]);
         /*std::cout << "Rel fac: " << rel_fac << std::endl;
         for (int i= 0; i < 6; i++)
         {
            std::cout << V_arr[i] << " ";
         }*/
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
               gsl_matrix_set(matrix_data,j,i,V_arr[0]*rel_fac);
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               gsl_matrix_set(matrix_data,j,i,V_arr[1]*rel_fac);
            }
         } else 
         {
            // The matrix is constructed as [[mm,mp],[pm,pp]]
            gsl_matrix_set(matrix_data,j,i,V_arr[3]*rel_fac); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
            gsl_matrix_set(matrix_data,j,i+(mom_grid_size_+1),V_arr[5]*rel_fac); //mp
            gsl_matrix_set(matrix_data,j+(mom_grid_size_+1),i,V_arr[4]*rel_fac); //pm
            gsl_matrix_set(matrix_data,j+(mom_grid_size_+1),i+(mom_grid_size_+1),V_arr[2]*rel_fac); //pp
         }
      }
   }
   
   // Return the matrix
   return matrix_data;
}

void print_gsl_matrix(gsl_matrix* matrix)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < matrix->size1; i++)
   {
      for (std::size_t j = 0; j < matrix->size1; j++)
      {
         std::cout << gsl_matrix_get(matrix,i,j) << " ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}

void Potential_mwpc::clear_saved_matrices()
{
   #ifdef ENABLE_DEBUG
      std::cerr << "clear_saved_matrices()" << std::endl;
   #endif
   // Go through the map and remove the created matrices
   for (std::map<qs::quantum_channel, std::map<std::string, gsl_matrix*>, qs::comp>::iterator it1=saved_matrices_.begin(); it1!=saved_matrices_.begin(); ++it1)
   {
      for (std::map<std::string, gsl_matrix*>::iterator it=it1->second.begin(); it!=it1->second.end(); ++it)
      {
         gsl_matrix_free(it->second);
      }
   }
   #ifdef ENABLE_DEBUG
      std::cerr << "clear_saved_matrices() - end" << std::endl;
   #endif
   
}

/*
   This function popolates the saved matrices in the given channel.
*/
void Potential_mwpc::populate_saved_mtx(qs::quantum_channel chn, bool rel_correction)
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
   gsl_matrix* mtx_const = get_matrix(0.0,chn,rel_correction); // Do not care about on-shell part
   //print_gsl_matrix(mtx_const);
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
      gsl_matrix* matrix = get_matrix(0.0,chn,rel_correction); // Do not care about on-shell part
      
      // Subtract the constant matrix
      gsl_matrix_sub(matrix,mtx_const);

      #ifdef ENABLE_DEBUG
         std::cerr << std::endl << "populate_saved_mtx(): Printing matrix associated to LEC: " << LECs_in_use_[i] << "." << std::endl;
         print_gsl_matrix(matrix);
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

void p_m(gsl_matrix* matrix)
{
   std::cout << "---------" << std::endl;
   for (std::size_t i = 0; i < matrix->size1; i++)
   {
      for (std::size_t j = 0; j < matrix->size1; j++)
      {
         std::cout << gsl_matrix_get(matrix,i,j) << " ";
      }   
      std::cout << std::endl;
   }
   std::cout << "---------" << std::endl;
}

gsl_matrix* Potential_mwpc::get_saved_matrix(double q_on_shell, qs::quantum_channel chn, bool rel_correction)
{
   unsigned int J = chn.J;
   unsigned int S = chn.S;
   bool coupled = chn.coupled;

   double mu;
   if (chn.tz == -1)
   {
      mu = constants::Mn/2.0; // nn
   } else if (chn.tz == 0)
   {
      mu = constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np
   } else if (chn.tz == 1)
   {
      mu = constants::Mp/2.0; // pp
   } else 
   {
      #ifdef ENABLE_DEBUG
         std::cerr << "Error in solve_in_chn(): Unknown isospin" << std::endl;
      #endif
   }
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

   
   gsl_matrix* matrix_saved_sum;
   if (!coupled)
   {
      matrix_saved_sum       = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
      gsl_matrix* tmp_matrix = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
      gsl_matrix_set_zero(matrix_saved_sum);
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
         gsl_matrix_memcpy(tmp_matrix, saved_matrices_[chn][LECs_in_use_[i]]);
         //p_m(saved_matrices_[chn][LECs_in_use_[i]]);
         // Scale tmp_matrix by the correct LEC
         gsl_matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);
         //std::cout << "LEC=" << LECs_[LECs_in_use_[i]] << std::endl;
         // Add this contribution to the sum

         gsl_matrix_add(matrix_saved_sum,tmp_matrix);
      }
      // Add constant
      gsl_matrix_add(matrix_saved_sum,saved_matrices_[chn]["const"]);
      // Delete the tmp matrix
      gsl_matrix_free(tmp_matrix);

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
         calc_element_V_arr(p_in,p_out, coupled, J,&tmp_arr[0]);

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
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/450,6))*exp(-gsl_pow_uint(p_out/450,6));
         rel_fac *= cutoff_regulator;

         if (S==0)
         {
            // Take element zero and insert it into the matrix
            gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[0]*rel_fac);
            gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[0]*rel_fac);
            
         } else 
         {
            // Take element one and insert it into the matrix
            gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[1]*rel_fac);
            gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[1]*rel_fac);
         }
      }

   } else 
   {
      matrix_saved_sum       = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
      gsl_matrix* tmp_matrix = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
      gsl_matrix_set_zero(matrix_saved_sum);
       for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
         gsl_matrix_memcpy(tmp_matrix, saved_matrices_[chn][LECs_in_use_[i]]);

         // Scale tmp_matrix by the correct LEC
         //std::cout << "LEC=" << LECs_[LECs_in_use_[i]] << std::endl;
         gsl_matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);

         // Add this contribution to the sum
         gsl_matrix_add(matrix_saved_sum,tmp_matrix);
      }
      // Add constant
      gsl_matrix_add(matrix_saved_sum,saved_matrices_[chn]["const"]);
     
      // Delete the tmp matrix
      gsl_matrix_free(tmp_matrix);

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
         calc_element_V_arr(p_in,p_out, coupled, J,&tmp_arr[0]);

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
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/450.0,6))*exp(-gsl_pow_uint(p_out/450.0,6));
         rel_fac *= cutoff_regulator;

         // In this part of the code the symmetry of the potential in momentum is taken
         // advantage of. Note that this does NOT mean that every block is symmetric,
         // but rather the WHOLE matrix.

         // Take mm element and insert it into the matrix
         // mm
         gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[3]*rel_fac); // Column
         gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[3]*rel_fac); // Row

         // Take mp element one and insert it into the matrix
         // mp 
         gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i+mom_grid_size_+1,tmp_arr[4]*rel_fac); // Column
         gsl_matrix_set(matrix_saved_sum,i,2*mom_grid_size_+1,tmp_arr[5]*rel_fac); // Row 

         // Take pm element one and insert it into the matrix
         // pm
         gsl_matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i,tmp_arr[5]*rel_fac); // Column
         gsl_matrix_set(matrix_saved_sum,i+mom_grid_size_+1,mom_grid_size_,tmp_arr[4]*rel_fac); // Row 

         // Take pp element one and insert it into the matrix
         // pp
         gsl_matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i+mom_grid_size_+1,tmp_arr[2]*rel_fac); // Column
         gsl_matrix_set(matrix_saved_sum,i+mom_grid_size_+1,2*mom_grid_size_+1,tmp_arr[2]*rel_fac); // Row 
      }
   }

   #ifdef ENABLE_DEBUG
      std::cout << std::endl << "get_saved_matrix(): Printing matrix." << std::endl;
      print_gsl_matrix(matrix_saved_sum);
   #endif
     
   // Return pointer to the full potential matrix.

   //std::cout << gsl_matrix_get(matrix_saved_sum,3,3) << std::endl;
   return matrix_saved_sum;
}

gsl_matrix* Potential_mwpc::get_matrix_no_onshell(qs::quantum_channel chn, bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "get_matrix()" << std::endl;
   #endif
   double mu;
   if (chn.tz == -1)
   {
      mu = constants::Mn/2.0; // nn
   } else if (chn.tz == 0)
   {
      mu = constants::Mn*constants::Mp/(constants::Mn+constants::Mp); // np
   } else if (chn.tz == 1)
   {
      mu = constants::Mp/2.0; // pp
   } else 
   {
      #ifdef ENABLE_DEBUG
         std::cerr << "Error in solve_in_chn(): Unknown isospin" << std::endl;
      #endif
   }

   // Allocate gsl matrices in the case of coupled and uncoupled channels.
   // The matrix becomes twise as large in the coupled case
   // In the construction of the saved matrix the desired on-shell momentum 
   // that the matrix element will be evaluated on are unknown. These matrix elements needs
   // to be computed at runtime
   gsl_matrix* matrix_data;
   if (chn.coupled) {
      matrix_data = gsl_matrix_alloc((2*mom_grid_size_),(2*mom_grid_size_));
   } else {
      matrix_data = gsl_matrix_alloc(mom_grid_size_,mom_grid_size_);
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
         double p_in  = 0;
         double p_out = 0;

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
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/450,6))*exp(-gsl_pow_uint(p_out/450,6));
         rel_fac *= cutoff_regulator;
      
         calc_element_V_arr(p_in,p_out,chn.coupled,chn.J,&V_arr[0]);
     
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
               gsl_matrix_set(matrix_data,j,i,V_arr[0]*rel_fac);
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               gsl_matrix_set(matrix_data,j,i,V_arr[1]*rel_fac);
            }
         } else 
         {
            // The matrix is constructed as [[mm,mp],[pm,pp]]
            gsl_matrix_set(matrix_data,j,i,V_arr[3]*rel_fac); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
            gsl_matrix_set(matrix_data,j,i+(mom_grid_size_),V_arr[5]*rel_fac); //mp
            gsl_matrix_set(matrix_data,j+(mom_grid_size_),i,V_arr[4]*rel_fac); //pm
            gsl_matrix_set(matrix_data,j+(mom_grid_size_),i+(mom_grid_size_),V_arr[2]*rel_fac); //pp
         }
      }
   }
   // Return the matrix
   return matrix_data;
}
