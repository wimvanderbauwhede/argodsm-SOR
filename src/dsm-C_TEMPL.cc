#include "dsm-api.h"

extern "C" {

void dsminitf_(int64_t* mem_sz, int64_t* cache_sz) {
	// std::cout << "argo::init("<<(std::size_t)(*mem_sz)/1024/1024<<","<<(std::size_t)(*cache_sz)/1024/1024<<")\n";
	argo::init((std::size_t)(*mem_sz),(std::size_t)(*cache_sz));
}

void dsmfinalisef_() {
	argo::finalize();
}

void dsmgetnnodesf_(int* n_nodes) {
	*n_nodes = argo::number_of_nodes();    
}
void dsmgetnidf_(int64_t* n_id) {
	*n_id = argo::node_id();
}

void dsmbarrierf_() {
	argo::barrier();
}

	// Initialize the lock
void dsminitlockf_(Ptr dsm_lock_ptr, Ptr dsm_lock_field_ptr)	{
	argo::globallock::global_tas_lock::internal_field_type* lock_field = fromWord<argo::globallock::global_tas_lock::internal_field_type*>(*dsm_lock_field_ptr);
	argo::globallock::global_tas_lock* lock = fromWord<argo::globallock::global_tas_lock*>(*dsm_lock_ptr);

	// argo::globallock::global_tas_lock::internal_field_type* 
	lock_field =
		argo::conew_<argo::globallock::global_tas_lock::internal_field_type>();
	// argo::globallock::global_tas_lock* 
	lock = new argo::globallock::global_tas_lock(lock_field);
	*dsm_lock_ptr = toWord<argo::globallock::global_tas_lock*>(lock);
	*dsm_lock_field_ptr = toWord<argo::globallock::global_tas_lock::internal_field_type*>(lock_field);
}

void dsmdeletelockf_(Ptr dsm_lock_ptr, Ptr dsm_lock_field_ptr)	{
	argo::globallock::global_tas_lock::internal_field_type* lock_field = fromWord<argo::globallock::global_tas_lock::internal_field_type*>(*dsm_lock_field_ptr);
	argo::globallock::global_tas_lock* lock = fromWord<argo::globallock::global_tas_lock*>(*dsm_lock_ptr);
	delete lock;
	argo::codelete_(lock_field);
}

void dsmlockf_(Ptr dsm_lock_ptr) {
	argo::globallock::global_tas_lock* lock = fromWord<argo::globallock::global_tas_lock*>(*dsm_lock_ptr);
	lock->lock();
}
	
void dsmunlockf_(Ptr dsm_lock_ptr) {
	argo::globallock::global_tas_lock* lock = fromWord<argo::globallock::global_tas_lock*>(*dsm_lock_ptr);
	lock->unlock();
}

#pragma GEN WrapperSubs

#ifndef GEN_WrapperSubs
void dsmallocfloatarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	// std::cout << "argo::conew_array<float>("<<(std::size_t)(*a_size)/1024/1024<<")\n";
	float* dsm_a_ptr_f = argo::conew_array<float>((std::size_t)(*a_size));
	// std::cout << "done:"<< dsm_a_ptr_f <<"\n";
	// 	void* vp = (void*)(dsm_a_ptr_f);
	// 	std::cout << "done2\n";
	// int64_t ivp = (int64_t)vp;
	// std::cout << "done3:"<< ivp << "\n";
	*dsm_a_ptr = toWord<float*>(dsm_a_ptr_f);
	// std::cout << "Boom!\n";
}

void dsmdeletefloatarrayf_(DSMArrayPtr dsm_a_ptr) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

// This is risky: are C floats and Fortran reals the same?
void dsmallocfloatf_(float* f_ptr, float* f_initval) {
	f_ptr = argo::conew_<float>(*f_initval);
};
void dsmdeletefloatf_(float* f_ptr) {
	argo::codelete_<float>(f_ptr);
};

void dsmwritefloatarrayf_(DSMArrayPtr dsm_a_ptr, int* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadfloatarrayf_(DSMArrayPtr dsm_a_ptr, int* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}
#endif



} // extern "C"
