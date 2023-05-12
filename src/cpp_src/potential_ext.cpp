#include "potential_ext.h"

extern "C" {
    void cdbonn_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

// This function is not complete!!!
void cdbonn_correct_arg(double qi, double qo, bool coupled, int S, int J, int T, int Tz,  double* V_arr)
{
    int coup = (int)coupled;
    cdbonn_fort_interface(&qi, &qo, &coup, &S, &J, &T, &Tz, &V_arr[0]); 

}

extern "C" {
    void nijmegen_fort_interface(double *qi,
			  double *qo,
			  int *coup,
			  int *S,
			  int *J,
			  int *T,
			  int *Tz,
			  double *pot);
}

// This function is not complete!!!
void nijm_correct_arg(double qi, double qo, bool coupled, int S, int J, int T, int Tz,  double* V_arr)
{
    int coup = (int)coupled;
    nijmegen_fort_interface(&qi, &qo, &coup, &S, &J, &T, &Tz, &V_arr[0]); 

}



/* These declarations are neccessary to be able to have the declarations
 * in a separate .cpp file.
 */

template Potential_ext<gsl_matrix>::Potential_ext(double* p_grid, int p_grid_length, 
        double cutoff_Lambda, 
        void (*f)(double qi,double qo, bool coupled, int S, int J, int T, 
            int Tz, double* V_arr));

template Potential_ext<gsl_matrix>::~Potential_ext();

template gsl_matrix* Potential_ext<gsl_matrix>::get_matrix(
        double q_on_shell, qs::quantum_channel chn, bool rel_correction);

template gsl_matrix* Potential_ext<gsl_matrix>::get_matrix_no_onshell(
        qs::quantum_channel chn, bool rel_correction);


template void Potential_ext<gsl_matrix>::populate_saved_mtx(qs::quantum_channel chn, 
        bool rel_correction);

template gsl_matrix* Potential_ext<gsl_matrix>::get_saved_matrix(double q_on_shell, 
        qs::quantum_channel chn, bool rel_correction);

template void Potential_ext<gsl_matrix>::print_LECs_and_params_info();

template <class gsl_m>
Potential_ext<gsl_m>::Potential_ext(double* p_grid, int p_grid_length, double cutoff_Lambda, 
        void (*f)(double qi,double qo, bool coupled, int S, int J, int T, int Tz, double* V_arr))
{
    // Set the function opinter to the correct function
    my_element_V_arr = f; 

    p_grid_ = p_grid;
    mom_grid_size_ = p_grid_length;
    cutoff_Lambda_ = cutoff_Lambda;
}
template <class gsl_m>
Potential_ext<gsl_m>::~Potential_ext()
{
}

template <class gsl_m>
gsl_matrix* Potential_ext<gsl_m>::get_matrix(double q_on_shell, 
        qs::quantum_channel chn, bool rel_correction)
{
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
   for (int i = 0; i < (int)mom_grid_size_+1; i++)
   {
      for (int j = 0; j < (int)mom_grid_size_+1; j++)
      {
         double V_arr[6]; // Array for data
         
         // Outgoing momentum is row index
         double p_in  = 0;
         double p_out = 0;

         if (j < (int)mom_grid_size_) {
            p_in  = p_grid_[j];
         } else {
            p_in = q_on_shell;
         }
         if (i < (int)mom_grid_size_) {
            p_out = p_grid_[i];
         } else {
            p_out = q_on_shell;
         }
         // THESE LINES DIFFER FROM THE FUNCTION IN POT_NN_MWPC
         //
         // Compute cutoff
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/cutoff_Lambda_,6))*exp(-gsl_pow_uint(p_out/cutoff_Lambda_,6));
      
         //std::cout << " LECS: " << LECs_["gA2"] << " " << LECs_["C1S0"] << " " << LECs_["C3S1"] << std::endl;
         
         //std::cout << chn.coupled << " " << chn.J << " " <<  chn.S << " " << chn.T << std::endl << " ------- " << std::endl;
         my_element_V_arr(p_in,p_out,chn.coupled,chn.S, chn.J, chn.T, chn.Tz, &V_arr[0]);
         //std::cout << "-----" << std::endl;
         //---------------------------------------------------
         //std::cout << "Rel fac: " << rel_fac << std::endl;
         /*for (int i= 0; i < 6; i++)
         {
            std::cout << V_arr[i] << " ";
         }
         std::cout << std::endl;*/
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
               gsl_matrix_set(matrix_data,j,i,V_arr[0]*cutoff_regulator);
              // std::cout << "Pot el S0: " << std::endl;
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               if (chn.J != 0)
               {
                  gsl_matrix_set(matrix_data,j,i,V_arr[1]*cutoff_regulator);
                //  std::cout << "1" << std::endl;
               } else // For J=0,S=1,L=1 case
               {
                  gsl_matrix_set(matrix_data,j,i,V_arr[2]*cutoff_regulator); // Take pp element to get L=1

                  //std::cout << "2" << std::endl;
               }
               
            }
         } else 
         {
             
            //std::cout << "3" << std::endl;
            // The matrix is constructed as [[mm,mp],[pm,pp]]
            gsl_matrix_set(matrix_data,j,i,V_arr[3]*cutoff_regulator); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
            gsl_matrix_set(matrix_data,j,i+(mom_grid_size_+1),V_arr[5]*cutoff_regulator); //mp
            gsl_matrix_set(matrix_data,j+(mom_grid_size_+1),i,V_arr[4]*cutoff_regulator); //pm
            gsl_matrix_set(matrix_data,j+(mom_grid_size_+1),i+(mom_grid_size_+1),V_arr[2]*cutoff_regulator); //pp
         }
      }
   }
   
   // Return the matrix
   return matrix_data;
}
    
template <class gsl_m>
gsl_matrix* Potential_ext<gsl_m>::get_matrix_no_onshell(qs::quantum_channel chn,
        bool rel_correction)
{
   #ifdef ENABLE_DEBUG
      std::cerr << "get_matrix_no_onshell()" << std::endl;
   #endif

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
   for (int i = 0; i < (int)mom_grid_size_; i++)
   {
      for (int j = 0; j < (int)mom_grid_size_; j++)
      {
         double V_arr[6]; // Array for data
         
         // Outgoing momentum is row index
         double p_in  = p_grid_[j];
         double p_out = p_grid_[i];

         // Compute relativistic factors
         double cutoff_regulator = exp(-gsl_pow_uint(p_in/cutoff_Lambda_,6))
             *exp(-gsl_pow_uint(p_out/cutoff_Lambda_,6));
         
         my_element_V_arr(p_in,p_out,chn.coupled,chn.S, chn.J, chn.T, chn.Tz, &V_arr[0]);
        
         // This is the same as in pot_nn_mwpc
         // ---------------------------------- 
         if (!chn.coupled)
         {
            if (chn.S==0) 
            {
               // Take S=0 element of V_arr and multiply by the relativistic factor
               gsl_matrix_set(matrix_data,j,i,V_arr[0]*cutoff_regulator);
            } else if (chn.S==1)
            {
               // Take S=1 element of V_arr
               if (chn.J != 0)
               {
                  gsl_matrix_set(matrix_data,j,i,V_arr[1]*cutoff_regulator);
               } else // For J=0,S=1,L=1 case
               {
                  gsl_matrix_set(matrix_data,j,i,V_arr[2]*cutoff_regulator); // Take pp element to get L=1
               }
            }

         } else 
         {
            // The matrix is constructed as [[mm,mp],[pm,pp]]
            gsl_matrix_set(matrix_data,j,i,V_arr[3]*cutoff_regulator); //mm
            // Offsett with mom_grid_size_+1, sinze the one is for the
            // on-shell part of the matrix that will be added later
            
            //std::cout << "element=" << V_arr[5] << " rel_fac=" << rel_fac << std::endl;
            //std::cout << p_in << " " << p_out << std::endl;
            gsl_matrix_set(matrix_data,j,i+(mom_grid_size_),V_arr[5]*cutoff_regulator); //mp
            gsl_matrix_set(matrix_data,j+(mom_grid_size_),i,V_arr[4]*cutoff_regulator); //pm
            gsl_matrix_set(matrix_data,j+(mom_grid_size_),i+(mom_grid_size_),V_arr[2]*cutoff_regulator); //pp
         }
         // ---------------------------------
      }
   }
   return matrix_data;
}
 
template <class gsl_m>
void Potential_ext<gsl_m>::populate_saved_mtx(qs::quantum_channel chn, 
        bool rel_correction)
{
}

template <class gsl_m>
gsl_m* Potential_ext<gsl_m>::get_saved_matrix(double q_on_shell, 
        qs::quantum_channel chn, bool rel_correction)
{
    return get_matrix(q_on_shell, chn, rel_correction);
}

template <class gsl_m>
void Potential_ext<gsl_m>::print_LECs_and_params_info()
{
    std::cout << "Potential_ext do not have any LECs or params." << std::endl;
}
