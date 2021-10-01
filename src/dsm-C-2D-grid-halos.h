#include "dsm-api.h"

void indexing_case1(DSM3DArray& a, int* i,int* j,int* k,int* lin_idx) ;

void indexing_case2(DSM3DArray& a, int i, int j, int k, int64_t* part_idx, int64_t* lin_idx);

void indexing_case3(DSM3DArrayH& a, int i, int j, int k, int64_t* part_idx, int64_t* sel, int64_t* comb_lin_idx);

// These are the routines to write to and read from a Case-2 3-D DSM array
void dsmWrite3DReal4Array(DSM3DArray& a, int64_t i, int64_t j, int64_t k, float* val);
float dsmRead3DReal4Array(DSM3DArray& a, int64_t i, int64_t j, int64_t k) ;

void dsmWriteArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k, float* val) ;
float dsmReadArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k) ;
