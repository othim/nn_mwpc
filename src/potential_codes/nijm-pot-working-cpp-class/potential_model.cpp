
#include "potential_model.h"

//#include "../constants.h"

/* Include present potentials */
#include "nijmegen.h"


potential_model::potential_model(){
}

int find_Tz_from_system(std::string system){
	int Tz;
	if (system=="np"){
		Tz = 0;
	}
	else if (system=="pp"){
		Tz = 1;
	}
	else if (system=="nn"){
		Tz = -1;
	}
	else{
		std::cout << "Invalid system entered. Exiting ..." << std::endl;
		exit(0);
	}
	
	return Tz;
}

potential_model *potential_model::fetch_potential_ptr(std::string model, std::string system){
	
	int Tz = find_Tz_from_system(system);
	
	if (model=="nijmegen"){
		nijmegen* pot_ptr = new nijmegen();
		
		pot_ptr->setSystem(Tz);

		return pot_ptr;
	}
	else{
		std::cout << "Invalid potential model entered. Exiting ..." << std::endl;
		exit(-1);
	}
}
