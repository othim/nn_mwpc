#include "pot_nn_mwpc.h"

// This struct is to define what set of parameters the function f_int want
struct my_f_params { double qi; double qo; int J; int l; Term* term; Potential_mwpc* this_pot;};


// Constructor
Potential_mwpc::Potential_mwpc(std::vector<std::string> terms, unsigned int N_GLI_PWA,double* p_grid, double* w_grid, std::size_t mom_grid_size)
{
   // Init constants
   N_GLI_PWA_ = N_GLI_PWA;
   p_grid_ = p_grid;
   w_grid_ = w_grid;
   mom_grid_size_ = mom_grid_size;
   
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

   #ifdef ENABLE_DEBUG
      std::cout << "V object created" << std::endl;
   #endif
}

// Destructor
Potential_mwpc::~Potential_mwpc()
{
   // Free all memory allocations
   gsl_integration_fixed_free(int_ang_);
   #ifdef ENABLE_DEBUG
      std::cout << "V object deleted" << std::endl;
   #endif
}

/* 
   ** Some helper functions to opep_get_el **
   ******************************************
*/


/*
   V_\alpha from Erkelenz without the isospin factor
*/
double Potential_mwpc::pot_OPEP_mom(double qi, double qo, double z)
{
   double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	return -(LECs_["gA2"]/(4.0*constants::fpi*constants::fpi))*(1.0/(q2+constants::mpi*constants::mpi));
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

double f_int_helper(double z, void* p)
{
   // Decode the void* to the parameters
   struct my_f_params * params = (struct my_f_params *)p;
   double qi = (params->qi);
   double qo = (params->qo);
   int J = (params->J);
   int l = (params->l);
   
   return (params->term)->get_v_alpha(qi,qo,z,params->this_pot->LECs_)*gsl_pow_int(z,l)*gsl_sf_legendre_Pl(J, z);
}

double Potential_mwpc::compute_A_integral(double qi, double qo, int J,int l,Term* term)
{
   // Define function to integrate
   gsl_function F;
   F.function = &f_int_helper;
   struct my_f_params params = {qi, qo, J, l, term, this};
   F.params = &params;

   double result;
   // Perform the integration from -1 to 1
   gsl_integration_fixed(&F,&result,int_ang_);

   return result*M_PI;
}

void Potential_mwpc::calc_element_V_arr(double qi,double qo, bool coupled, int J, double* V_arr)
{
   for (std::size_t i = 0; i < terms_in_pot_.size(); i++)
   {  
      if (!terms_in_pot_[i].is_lec())
      {
         // Calculate the A_x integrals (don't calculate all)
         std::clock_t start, end;
         start = std::clock();

         double A_0 = compute_A_integral(qi,qo,J,0,&terms_in_pot_[i]);
         double A_1 = 0;
         double A_P = compute_A_integral(qi,qo,J+1,0,&terms_in_pot_[i]);
         double A_M = 0;

         if (J!=0) {
            A_M = compute_A_integral(qi,qo,J-1,0,&terms_in_pot_[i]);
         }
         if (!coupled) {
            A_1 = compute_A_integral(qi,qo,J,1,&terms_in_pot_[i]);
         }

         end = std::clock();
         std::cout << "Time taken to compute A integrals is : " << (double)(end-start)/(double)CLOCKS_PER_SEC; 
         std::cout << " sec " << std::endl;

         // Call pwa with the correct spin structure from this term
         // This will fill up the array V_arr with the correct potential elements
         
         start = std::clock();

         pwa(qi,qo,coupled,J,A_M,A_P,A_0,A_1,terms_in_pot_[i].get_spin_structure(),terms_in_pot_[i].get_isovector(),V_arr);

         end = std::clock();
         std::cout << "Time taken to call pwa is : " << (double)(end-start)/(double)CLOCKS_PER_SEC; 
         std::cout << " sec " << std::endl;
       
      } else 
      {
         // Define some variables and initialize V_arr to all zeros.
         double V_uncoupled_S0 = 0;
         double V_uncoupled_S1 = 0;
         double V_coupled_mm   = 0;
         double V_coupled_pm   = 0;
         double V_coupled_mp   = 0;
         double V_coupled_pp   = 0;

         // Get JLS from the Term to know where to apend the lec value 
         LS_term LS_term = terms_in_pot_[i].get_LS_term();
         // Find which ouput to append (if any)
         if (J==LS_term.J)
         {
            if (coupled)
            {
               if (LS_term.Li == J-1 && LS_term.Lo == J-1) // --
               {
                  V_coupled_mm = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J+1 && LS_term.Lo == J+1) // ++
               {
                  V_coupled_pp = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // -+
               {
                  V_coupled_mp = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.Li == J-1 && LS_term.Lo == J+1) // +-
               {
                  V_coupled_pm = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               }
            } else
            {
               if (LS_term.S == 0) // S0
               {
                  V_uncoupled_S0 = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               } else if (LS_term.S == 1) // S1
               {
                  V_uncoupled_S1 = terms_in_pot_[i].get_LEC_element(qi,qo,LECs_);
               }
            }
         }
         V_arr[0] = V_uncoupled_S0;
         V_arr[1] = V_uncoupled_S1;
         V_arr[2] = V_coupled_pp;
         V_arr[3] = V_coupled_mm;
         V_arr[4] = V_coupled_pm;
         V_arr[5] = V_coupled_mp;
      }
   }
}

void Potential_mwpc::pwa(double qi,double qo, bool coupled, int J,double A_M,double A_P,double A_0,double A_1,std::string spin_struct,bool isovector,double* V_arr)
{
   // Define some variables
   double V_uncoupled_S0 = 0;
   double V_uncoupled_S1 = 0;
   double V_coupled_mm   = 0;
   double V_coupled_pm   = 0;
   double V_coupled_mp   = 0;
   double V_coupled_pp   = 0;

   if (spin_struct == "tensor")
   {
      // Check which elements that are non-zero by checking if the Channel is
      // coupled or not
      if (!coupled)
      {
         // OPEP uncoupled interactions
         V_uncoupled_S0    = 2 * (-(qo*qo+qi*qi)*A_0 + 2*qo*qi*A_1);
         if (J!=0) {
            V_uncoupled_S1 = 2 * ((qo*qi+qi*qi)*A_0 - 2*qo*qi*(1.0/(2.0*J+1.0))*(J*A_P + (J+1)*A_M));
         }
      } else
      {
         V_coupled_pp = (2.0/(2.0*J+1.0)) * (-(qo*qo+qi*qi)*A_P + 2*qo*qi*A_0);
         if (J!= 0)
         {
            V_coupled_mm = (2./(2*J+1)) * ((qo*qo+qi*qi)*A_M - 2*qo*qi*A_0);
            V_coupled_mp = (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_P + qo*qo*A_M - 2*qo*qi*A_0);
            V_coupled_pm = (4*sqrt(J*(J+1))/(2.0*J+1.0)) * (qi*qi*A_M + qo*qo*A_P - 2*qo*qi*A_0);
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
    } else
    {
       std::cerr << "Unknown <spin_struct> passed to function <pwa> in Potential_mwpc object." << std::endl;
    }
      // Fill the output array with the correct values
      V_arr[0] = V_uncoupled_S0;
      V_arr[1] = V_uncoupled_S1;
      V_arr[2] = V_coupled_pp;
      V_arr[3] = V_coupled_mm;
      V_arr[4] = V_coupled_pm;
      V_arr[5] = V_coupled_mp;
}

gsl_matrix* Potential_mwpc::get_matrix(bool coupled, int J,int S, bool rel_correction, double regulator)
{
   // Allocate gsl matrices in the case of coupled and uncoupled channels.
   // The matrix becomes twise as large in the coupled case
   // In the construction of the saved matrix the desired on-shell momentum 
   // that the matrix element will be evaluated on are unknown. These matrix elements needs
   // to be computed at runtime
   gsl_matrix* matrix_data;
   if (coupled) {
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
   for (std::size_t i = 0; i < mom_grid_size_; i++)
   {
      for (std::size_t j = 0; j < mom_grid_size_; j++)
      {
         double V_arr[6]; // Array for data
         // Outgoing momentum is row index
         calc_element_V_arr(p_grid_[j],p_grid_[i],coupled,J,V_arr);

         if (!coupled)
         {
            if (S==0) 
            {
               // Take S=0 element of V_arr
               gsl_matrix_set(matrix_data,i,j,V_arr[0]);
            } else if (S==1)
            {
               // Take S=1 element of V_arr
               gsl_matrix_set(matrix_data,i,j,V_arr[0]);
            }
         } else 
         {
            // The matrix is constructed as [[mm,mp],[pm,pp]]
            gsl_matrix_set(matrix_data,i,j,V_arr[2]); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            gsl_matrix_set(matrix_data,i,j+(mom_grid_size_+1),V_arr[3]); //mp
            gsl_matrix_set(matrix_data,i+(mom_grid_size_+1),j,V_arr[4]); //pm
            gsl_matrix_set(matrix_data,i+(mom_grid_size_+1),j+(mom_grid_size_+1),V_arr[5]); //pp
         }
      }
   }
   
   // Add relativistic factors and cutoff

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
   // Go through the map and remove the created matrices
   for (std::map<std::string, gsl_matrix*>::iterator it=saved_matrices_.begin(); it!=saved_matrices_.end(); ++it)
   {
      gsl_matrix_free(it->second);
   }
}

void Potential_mwpc::populate_saved_mtx(double qi,double qo, bool coupled, int J, int S, bool rel_correction, bool cutoff_on)
{
   clear_saved_matrices(); // Clears the allocated pointers

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
      gsl_matrix* matrix = get_matrix(coupled, J, S,rel_correction, cutoff_on); 

      #ifdef ENABLE_DEBUG
         std::cout << std::endl << "populate_saved_mtx(): Printing matrix associated to LEC: " << LECs_in_use_[i] << "." << std::endl;
         print_gsl_matrix(matrix);
      #endif
      // Save the matrix as the matrix corresponding tho the non-zero LEC.
      saved_matrices_[LECs_in_use_[i]] = matrix;
   }
   LECs_[LECs_in_use_[LECs_in_use_.size()-1]] = 0.0;
   // After the function has populated the matrices the default is that all LECs are set to zero, 
   // this is to alert the user if the LECs are not updater before getting the matrix.
}

void Potential_mwpc::get_saved_matrix(double q_on_shell, bool coupled, int J, int S, bool rel_correction, bool cutoff_on,gsl_matrix* out_matrix)
{
   // Note that the matrix is computed for the current LECs_!!
   // Note also that the LECs_ are screwed up by the act of saving the matrices!
   // This means that the 
   
   // Check if saved otherwise print error message
   if (saved_matrices_.empty())
   {
      std::cerr << "get_saved_matrix(): There is no saved matrices, returning." << std::endl;
      return;
   }

   // Load the matrices and compute the necessary matrix elements on the diagonal

   
   gsl_matrix* matrix_saved_sum;
   if (!coupled)
   {
      matrix_saved_sum = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
      gsl_matrix* tmp_matrix       = gsl_matrix_alloc(mom_grid_size_ + 1,mom_grid_size_ + 1);
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
         gsl_matrix_memcpy(tmp_matrix, saved_matrices_[LECs_in_use_[i]]);

         // Scale tmp_matrix by the correct LEC
         gsl_matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);

         // Add this contribution to the sum
         gsl_matrix_add(matrix_saved_sum,tmp_matrix);
      }
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
         calc_element_V_arr(p_grid_[i],q_on_shell, coupled, J,tmp_arr);
         if (S==0)
         {
            // Take element zero and insert it into the matrix
            gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[0]);
            gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[0]);
         } else 
         {
            // Take element one and insert it into the matrix
            gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[1]);
            gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[1]);
         }
      }

   } else 
   {
      matrix_saved_sum       = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
      gsl_matrix* tmp_matrix = gsl_matrix_alloc(2*mom_grid_size_ + 2,2*mom_grid_size_ + 2);
      for (std::size_t i = 0; i < LECs_in_use_.size(); i++)
      {
         // Copy saved matrix to not mess it upp
         gsl_matrix_memcpy(tmp_matrix, saved_matrices_[LECs_in_use_[i]]);

         // Scale tmp_matrix by the correct LEC
         gsl_matrix_scale(tmp_matrix, LECs_[LECs_in_use_[i]]);

         // Add this contribution to the sum
         gsl_matrix_add(matrix_saved_sum,tmp_matrix);
      }
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
         calc_element_V_arr(p_grid_[i],q_on_shell, coupled, J,tmp_arr);

         // Take mm element and insert it into the matrix
         gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i,tmp_arr[2]); // Column
         gsl_matrix_set(matrix_saved_sum,i,mom_grid_size_,tmp_arr[2]); // Row

         // Take mp element one and insert it into the matrix
         gsl_matrix_set(matrix_saved_sum,mom_grid_size_,i+mom_grid_size_+1,tmp_arr[3]); // Column
         gsl_matrix_set(matrix_saved_sum,i,2*mom_grid_size_+1,tmp_arr[3]); // Row 

         // Take pm element one and insert it into the matrix
         gsl_matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i,tmp_arr[3]); // Column
         gsl_matrix_set(matrix_saved_sum,i+mom_grid_size_+1,mom_grid_size_,tmp_arr[3]); // Row 

         // Take pp element one and insert it into the matrix
         gsl_matrix_set(matrix_saved_sum,2*mom_grid_size_+1,i+mom_grid_size_+1,tmp_arr[3]); // Column
         gsl_matrix_set(matrix_saved_sum,i+mom_grid_size_+1,2*mom_grid_size_+1,tmp_arr[3]); // Row 
      }
   }

   #ifdef ENABLE_DEBUG
      std::cout << std::endl << "get_saved_matrix(): Printing matrix." << std::endl;
      print_gsl_matrix(matrix_saved_sum);
   #endif
     
   // Return pointer to the full potential matrix.
   out_matrix = matrix_saved_sum;
}