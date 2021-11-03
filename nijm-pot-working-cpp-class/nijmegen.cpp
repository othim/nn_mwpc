

#include "nijmegen.h"

nijmegen::nijmegen(){
}

void nijmegen::setSystem(int Tz_in){
	Tz = -Tz_in;
}

void nijmegen::update_parameters(){}

void nijmegen::setParameter(int index, double value){
}

void nijmegen::retrieveSymmetric(bool &symmetric_in){
	symmetric_in = symmetric;
}

void nijmegen::retrieve_state_dep_labels(std::vector<std::string> &state_dep_labels_in){
	state_dep_labels_in = state_dep_labels;
}

void nijmegen::retrieve_state_dep_indices(std::vector<int> &state_dep_idx_lookup_in){
	state_dep_idx_lookup_in = state_dep_idx_lookup;
}

void nijmegen::retrieve_state_indep_labels(std::vector<std::string> &state_indep_labels_in){
	state_indep_labels_in = state_indep_labels;
}

void nijmegen::retrieve_state_indep_indices(std::vector<int> &state_indep_idx_lookup_in){
	state_indep_idx_lookup_in = state_indep_idx_lookup;
}

void nijmegen::retrieve_state_list(std::vector<int> &state_list_in){
	state_list_in = state_list;
}

void nijmegen::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, double *Varray){
	
	int coupledState = coupled;
	
	nijmegen_fort_interface(&qi, &qo, &coupledState, &S, &J, &T, &Tz, Varray);
}
