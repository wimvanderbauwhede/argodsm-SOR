// !!! Don't edit this file!!! Edit dsm-C_TEMPL.f95 and run ./gen_C_dsmAPI_subs.pl !!!
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


#define GEN_WrapperSubs


// Generated API routines



// Scalar float

void dsmallocfloatf_(DSMScalarPtr dsm_s_ptr, float* f_initval) {
	float* dsm_s_ptr_f = argo::conew_<float>(*f_initval);
	*dsm_s_ptr	= toWord<float*>(dsm_s_ptr_f);    
};

void dsmdeletefloatf_(DSMScalarPtr dsm_s_ptr) {
	float* dsm_s_ptr_f = fromWord<float*>(*dsm_s_ptr);
	argo::codelete_<float>(dsm_s_ptr_f);
};

void dsmwritefloatf_(DSMScalarPtr dsm_s_ptr,float* val) {
	float* dsm_s_ptr_f = fromWord<float*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreadfloatf_(DSMScalarPtr dsm_s_ptr,float* val) {
	float* dsm_s_ptr_f = fromWord<float*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of float

void dsmallocfloatarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	float* dsm_a_ptr_f = argo::conew_array<float>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<float*>(dsm_a_ptr_f);
}

void dsmdeletefloatarrayf_(DSMArrayPtr dsm_a_ptr) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwritefloatarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadfloatarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpyfloatarrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, float* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(float));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromfloatarrayf_(float* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	float* dsm_buf_ptr_f = fromWord<float*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(float));
}
*/

// Local array of float

void dsmallocfloatlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	float* dsm_a_ptr_f = argo::new_array<float>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<float*>(dsm_a_ptr_f);
}

void dsmdeletefloatlocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwritefloatlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadfloatlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, float* val) {
	float* dsm_a_ptr_f = fromWord<float*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        

// Scalar double

void dsmallocdoublef_(DSMScalarPtr dsm_s_ptr, double* f_initval) {
	double* dsm_s_ptr_f = argo::conew_<double>(*f_initval);
	*dsm_s_ptr	= toWord<double*>(dsm_s_ptr_f);    
};

void dsmdeletedoublef_(DSMScalarPtr dsm_s_ptr) {
	double* dsm_s_ptr_f = fromWord<double*>(*dsm_s_ptr);
	argo::codelete_<double>(dsm_s_ptr_f);
};

void dsmwritedoublef_(DSMScalarPtr dsm_s_ptr,double* val) {
	double* dsm_s_ptr_f = fromWord<double*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreaddoublef_(DSMScalarPtr dsm_s_ptr,double* val) {
	double* dsm_s_ptr_f = fromWord<double*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of double

void dsmallocdoublearrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	double* dsm_a_ptr_f = argo::conew_array<double>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<double*>(dsm_a_ptr_f);
}

void dsmdeletedoublearrayf_(DSMArrayPtr dsm_a_ptr) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwritedoublearrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, double* val) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreaddoublearrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, double* val) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpydoublearrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, double* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(double));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromdoublearrayf_(double* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	double* dsm_buf_ptr_f = fromWord<double*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(double));
}
*/

// Local array of double

void dsmallocdoublelocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	double* dsm_a_ptr_f = argo::new_array<double>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<double*>(dsm_a_ptr_f);
}

void dsmdeletedoublelocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwritedoublelocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, double* val) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreaddoublelocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, double* val) {
	double* dsm_a_ptr_f = fromWord<double*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        

// Scalar char

void dsmalloccharf_(DSMScalarPtr dsm_s_ptr, char* f_initval) {
	char* dsm_s_ptr_f = argo::conew_<char>(*f_initval);
	*dsm_s_ptr	= toWord<char*>(dsm_s_ptr_f);    
};

void dsmdeletecharf_(DSMScalarPtr dsm_s_ptr) {
	char* dsm_s_ptr_f = fromWord<char*>(*dsm_s_ptr);
	argo::codelete_<char>(dsm_s_ptr_f);
};

void dsmwritecharf_(DSMScalarPtr dsm_s_ptr,char* val) {
	char* dsm_s_ptr_f = fromWord<char*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreadcharf_(DSMScalarPtr dsm_s_ptr,char* val) {
	char* dsm_s_ptr_f = fromWord<char*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of char

void dsmallocchararrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	char* dsm_a_ptr_f = argo::conew_array<char>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<char*>(dsm_a_ptr_f);
}

void dsmdeletechararrayf_(DSMArrayPtr dsm_a_ptr) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwritechararrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, char* val) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadchararrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, char* val) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpychararrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, char* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(char));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromchararrayf_(char* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	char* dsm_buf_ptr_f = fromWord<char*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(char));
}
*/

// Local array of char

void dsmalloccharlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	char* dsm_a_ptr_f = argo::new_array<char>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<char*>(dsm_a_ptr_f);
}

void dsmdeletecharlocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwritecharlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, char* val) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadcharlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, char* val) {
	char* dsm_a_ptr_f = fromWord<char*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        

// Scalar short

void dsmallocshortf_(DSMScalarPtr dsm_s_ptr, short* f_initval) {
	short* dsm_s_ptr_f = argo::conew_<short>(*f_initval);
	*dsm_s_ptr	= toWord<short*>(dsm_s_ptr_f);    
};

void dsmdeleteshortf_(DSMScalarPtr dsm_s_ptr) {
	short* dsm_s_ptr_f = fromWord<short*>(*dsm_s_ptr);
	argo::codelete_<short>(dsm_s_ptr_f);
};

void dsmwriteshortf_(DSMScalarPtr dsm_s_ptr,short* val) {
	short* dsm_s_ptr_f = fromWord<short*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreadshortf_(DSMScalarPtr dsm_s_ptr,short* val) {
	short* dsm_s_ptr_f = fromWord<short*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of short

void dsmallocshortarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	short* dsm_a_ptr_f = argo::conew_array<short>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<short*>(dsm_a_ptr_f);
}

void dsmdeleteshortarrayf_(DSMArrayPtr dsm_a_ptr) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwriteshortarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, short* val) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadshortarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, short* val) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpyshortarrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, short* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(short));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromshortarrayf_(short* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	short* dsm_buf_ptr_f = fromWord<short*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(short));
}
*/

// Local array of short

void dsmallocshortlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	short* dsm_a_ptr_f = argo::new_array<short>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<short*>(dsm_a_ptr_f);
}

void dsmdeleteshortlocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwriteshortlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, short* val) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadshortlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, short* val) {
	short* dsm_a_ptr_f = fromWord<short*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        

// Scalar int

void dsmallocintf_(DSMScalarPtr dsm_s_ptr, int* f_initval) {
	int* dsm_s_ptr_f = argo::conew_<int>(*f_initval);
	*dsm_s_ptr	= toWord<int*>(dsm_s_ptr_f);    
};

void dsmdeleteintf_(DSMScalarPtr dsm_s_ptr) {
	int* dsm_s_ptr_f = fromWord<int*>(*dsm_s_ptr);
	argo::codelete_<int>(dsm_s_ptr_f);
};

void dsmwriteintf_(DSMScalarPtr dsm_s_ptr,int* val) {
	int* dsm_s_ptr_f = fromWord<int*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreadintf_(DSMScalarPtr dsm_s_ptr,int* val) {
	int* dsm_s_ptr_f = fromWord<int*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of int

void dsmallocintarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	int* dsm_a_ptr_f = argo::conew_array<int>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<int*>(dsm_a_ptr_f);
}

void dsmdeleteintarrayf_(DSMArrayPtr dsm_a_ptr) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwriteintarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int* val) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadintarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int* val) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpyintarrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, int* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(int));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromintarrayf_(int* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	int* dsm_buf_ptr_f = fromWord<int*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(int));
}
*/

// Local array of int

void dsmallocintlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	int* dsm_a_ptr_f = argo::new_array<int>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<int*>(dsm_a_ptr_f);
}

void dsmdeleteintlocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwriteintlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int* val) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadintlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int* val) {
	int* dsm_a_ptr_f = fromWord<int*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        

// Scalar long

void dsmalloclongf_(DSMScalarPtr dsm_s_ptr, long* f_initval) {
	long* dsm_s_ptr_f = argo::conew_<long>(*f_initval);
	*dsm_s_ptr	= toWord<long*>(dsm_s_ptr_f);    
};

void dsmdeletelongf_(DSMScalarPtr dsm_s_ptr) {
	long* dsm_s_ptr_f = fromWord<long*>(*dsm_s_ptr);
	argo::codelete_<long>(dsm_s_ptr_f);
};

void dsmwritelongf_(DSMScalarPtr dsm_s_ptr,long* val) {
	long* dsm_s_ptr_f = fromWord<long*>(*dsm_s_ptr);
	*dsm_s_ptr_f = *val;
}

void dsmreadlongf_(DSMScalarPtr dsm_s_ptr,long* val) {
	long* dsm_s_ptr_f = fromWord<long*>(*dsm_s_ptr);
	*val = *dsm_s_ptr_f;
}

// Array of long

void dsmalloclongarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	long* dsm_a_ptr_f = argo::conew_array<long>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<long*>(dsm_a_ptr_f);
}

void dsmdeletelongarrayf_(DSMArrayPtr dsm_a_ptr) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwritelongarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, long* val) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadlongarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, long* val) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO! 
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpylongarrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, long* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(long));
}
// copy a dsm buffer into a local array
void dsmmemcpyfromlongarrayf_(long* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	long* dsm_buf_ptr_f = fromWord<long*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(long));
}
*/

// Local array of long

void dsmalloclonglocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	long* dsm_a_ptr_f = argo::new_array<long>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<long*>(dsm_a_ptr_f);
}

void dsmdeletelonglocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwritelonglocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, long* val) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadlonglocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, long* val) {
	long* dsm_a_ptr_f = fromWord<long*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        
// Array of pointer

void dsmallocpointerarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	int64_t* dsm_a_ptr_f = argo::conew_array<int64_t>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<int64_t*>(dsm_a_ptr_f);
}

void dsmdeletepointerarrayf_(DSMArrayPtr dsm_a_ptr) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	argo::codelete_array(dsm_a_ptr_f);
}

void dsmwritepointerarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int64_t* val) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadpointerarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int64_t* val) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

// TODO!
/*
// Copy a local buffer into the dsm array, call this dsmmemcpyto?
void dsmmemcpypointerarrayf_(DSMArrayPtr dsm_a_ptr, int64_t * start_idx, int64_t* buf_ptr, int64_t * buf_start_idx, int64_t* buf_sz) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
    std::memcpy((void*)&dsm_a_ptr_f[*start_idx],(void*)&buf_ptr[*buf_start_idx],(*buf_sz)*sizeof(int64_t));
}
// copy a dsm buffer into a local array
void dsmmemcpyfrompointerarrayf_(int64_t* loc_a_ptr, int64_t * start_idx, DSMArrayPtr dsm_buf_ptr, int64_t * buf_start_idx, int64_t* dsm_buf_sz) {
	int64_t* dsm_buf_ptr_f = fromWord<int64_t*>(*dsm_buf_ptr);
    std::memcpy((void*)&loc_a_ptr_f[*start_idx],(void*)&dsm_buf_ptr_f[buf_start_idx],(*dsm_buf_sz)*sizeof(int64_t));
}
*/

// Local array of pointer

void dsmallocpointerlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* a_size) {
	int64_t* dsm_a_ptr_f = argo::new_array<int64_t>((std::size_t)(*a_size));
	*dsm_a_ptr = toWord<int64_t*>(dsm_a_ptr_f);
}

void dsmdeletepointerlocalarrayf_(DSMArrayPtr dsm_a_ptr) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	argo::delete_array(dsm_a_ptr_f);
}

void dsmwritepointerlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int64_t* val) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	dsm_a_ptr_f[*idx] = *val;
}
void dsmreadpointerlocalarrayf_(DSMArrayPtr dsm_a_ptr, int64_t* idx, int64_t* val) {
	int64_t* dsm_a_ptr_f = fromWord<int64_t*>(*dsm_a_ptr);
	*val = dsm_a_ptr_f[*idx];
}

        
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
