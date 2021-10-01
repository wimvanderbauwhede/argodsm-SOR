// !!! Don't edit this file!!! Edit dsmArray_TEMPL.f95 and run ./gen_C_dsmAPI_subs_header.pl !!!
/*
   This is a function-based ArgoDSM API intended for use in Fortran programs

   The approach is to cast the pointer to the array to a 64-bit integer and return it.
   This word gets passed around in the Fortran code (as integer(8)).
   Every call to the API takes this word as its first argument, casts it back to the object and so on.

   So we create some low-level casting functions first
*/
#ifndef _DSM_ARRAY_C_H_
#define _DSM_ARRAY_C_H_
#include <cstdint>
#include <cstring>
typedef int64_t* DSMArrayPtr;
typedef int64_t* DSMScalarPtr;
typedef int64_t* Ptr; 
// The minimal API is as follows:
void dsmgetnnodesf_(int*);
void dsmgetnidf_(int64_t*);
void dsminitf_(int64_t* , int64_t* ) ;
void dsmfinalisef_() ;

void dsmbarrierf_();
void dsminitlockf_(Ptr, Ptr );
void dsmdeletelockf_(Ptr, Ptr );
void dsmlockf_(Ptr);
void dsmunlockf_(Ptr);

void dsmmemcpyf_(DSMArrayPtr, int64_t*, void*, int64_t);


#define GEN_WrapperSubs


// Generated API routines


// Scalar float

void dsmallocfloatf_(DSMScalarPtr,float*) ;
void dsmdeletefloatf_(DSMScalarPtr);
void dsmwritefloatf_(DSMScalarPtr,float*);
void dsmreadfloatf_(DSMScalarPtr,float*);

// Array of float

void dsmallocfloatarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletefloatarrayf_(DSMArrayPtr) ;
void dsmwritefloatarrayf_(DSMArrayPtr,int64_t*,float*);
void dsmreadfloatarrayf_(DSMArrayPtr,int64_t*,float*);
// TODO! void dsmmemcpyfloatarrayf_(DSMArrayPtr, int64_t*, float*, int64_t*);

// Local array of float

void dsmallocfloatlocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletefloatlocalarrayf_(DSMArrayPtr) ;
void dsmwritefloatlocalarrayf_(DSMArrayPtr,int64_t*,float*);
void dsmreadfloatlocalarrayf_(DSMArrayPtr,int64_t*,float*);
// TODO! void dsmmemcpyfloatlocalarrayf_(DSMArrayPtr, int64_t*, float*, int64_t*);

        
// Scalar double

void dsmallocdoublef_(DSMScalarPtr,double*) ;
void dsmdeletedoublef_(DSMScalarPtr);
void dsmwritedoublef_(DSMScalarPtr,double*);
void dsmreaddoublef_(DSMScalarPtr,double*);

// Array of double

void dsmallocdoublearrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletedoublearrayf_(DSMArrayPtr) ;
void dsmwritedoublearrayf_(DSMArrayPtr,int64_t*,double*);
void dsmreaddoublearrayf_(DSMArrayPtr,int64_t*,double*);
// TODO! void dsmmemcpydoublearrayf_(DSMArrayPtr, int64_t*, double*, int64_t*);

// Local array of double

void dsmallocdoublelocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletedoublelocalarrayf_(DSMArrayPtr) ;
void dsmwritedoublelocalarrayf_(DSMArrayPtr,int64_t*,double*);
void dsmreaddoublelocalarrayf_(DSMArrayPtr,int64_t*,double*);
// TODO! void dsmmemcpydoublelocalarrayf_(DSMArrayPtr, int64_t*, double*, int64_t*);

        
// Scalar char

void dsmalloccharf_(DSMScalarPtr,char*) ;
void dsmdeletecharf_(DSMScalarPtr);
void dsmwritecharf_(DSMScalarPtr,char*);
void dsmreadcharf_(DSMScalarPtr,char*);

// Array of char

void dsmallocchararrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletechararrayf_(DSMArrayPtr) ;
void dsmwritechararrayf_(DSMArrayPtr,int64_t*,char*);
void dsmreadchararrayf_(DSMArrayPtr,int64_t*,char*);
// TODO! void dsmmemcpychararrayf_(DSMArrayPtr, int64_t*, char*, int64_t*);

// Local array of char

void dsmalloccharlocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletecharlocalarrayf_(DSMArrayPtr) ;
void dsmwritecharlocalarrayf_(DSMArrayPtr,int64_t*,char*);
void dsmreadcharlocalarrayf_(DSMArrayPtr,int64_t*,char*);
// TODO! void dsmmemcpycharlocalarrayf_(DSMArrayPtr, int64_t*, char*, int64_t*);

        
// Scalar short

void dsmallocshortf_(DSMScalarPtr,short*) ;
void dsmdeleteshortf_(DSMScalarPtr);
void dsmwriteshortf_(DSMScalarPtr,short*);
void dsmreadshortf_(DSMScalarPtr,short*);

// Array of short

void dsmallocshortarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeleteshortarrayf_(DSMArrayPtr) ;
void dsmwriteshortarrayf_(DSMArrayPtr,int64_t*,short*);
void dsmreadshortarrayf_(DSMArrayPtr,int64_t*,short*);
// TODO! void dsmmemcpyshortarrayf_(DSMArrayPtr, int64_t*, short*, int64_t*);

// Local array of short

void dsmallocshortlocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeleteshortlocalarrayf_(DSMArrayPtr) ;
void dsmwriteshortlocalarrayf_(DSMArrayPtr,int64_t*,short*);
void dsmreadshortlocalarrayf_(DSMArrayPtr,int64_t*,short*);
// TODO! void dsmmemcpyshortlocalarrayf_(DSMArrayPtr, int64_t*, short*, int64_t*);

        
// Scalar int

void dsmallocintf_(DSMScalarPtr,int*) ;
void dsmdeleteintf_(DSMScalarPtr);
void dsmwriteintf_(DSMScalarPtr,int*);
void dsmreadintf_(DSMScalarPtr,int*);

// Array of int

void dsmallocintarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeleteintarrayf_(DSMArrayPtr) ;
void dsmwriteintarrayf_(DSMArrayPtr,int64_t*,int*);
void dsmreadintarrayf_(DSMArrayPtr,int64_t*,int*);
// TODO! void dsmmemcpyintarrayf_(DSMArrayPtr, int64_t*, int*, int64_t*);

// Local array of int

void dsmallocintlocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeleteintlocalarrayf_(DSMArrayPtr) ;
void dsmwriteintlocalarrayf_(DSMArrayPtr,int64_t*,int*);
void dsmreadintlocalarrayf_(DSMArrayPtr,int64_t*,int*);
// TODO! void dsmmemcpyintlocalarrayf_(DSMArrayPtr, int64_t*, int*, int64_t*);

        
// Scalar long

void dsmalloclongf_(DSMScalarPtr,long*) ;
void dsmdeletelongf_(DSMScalarPtr);
void dsmwritelongf_(DSMScalarPtr,long*);
void dsmreadlongf_(DSMScalarPtr,long*);

// Array of long

void dsmalloclongarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletelongarrayf_(DSMArrayPtr) ;
void dsmwritelongarrayf_(DSMArrayPtr,int64_t*,long*);
void dsmreadlongarrayf_(DSMArrayPtr,int64_t*,long*);
// TODO! void dsmmemcpylongarrayf_(DSMArrayPtr, int64_t*, long*, int64_t*);

// Local array of long

void dsmalloclonglocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletelonglocalarrayf_(DSMArrayPtr) ;
void dsmwritelonglocalarrayf_(DSMArrayPtr,int64_t*,long*);
void dsmreadlonglocalarrayf_(DSMArrayPtr,int64_t*,long*);
// TODO! void dsmmemcpylonglocalarrayf_(DSMArrayPtr, int64_t*, long*, int64_t*);

        
// Array of pointer

void dsmallocpointerarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletepointerarrayf_(DSMArrayPtr) ;
void dsmwritepointerarrayf_(DSMArrayPtr,int64_t*,int64_t*);
void dsmreadpointerarrayf_(DSMArrayPtr,int64_t*,int64_t*);
// TODO! void dsmmemcpypointerarrayf_(DSMArrayPtr, int64_t*, int64_t*, int64_t*);

// Local array of pointer

void dsmallocpointerlocalarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdeletepointerlocalarrayf_(DSMArrayPtr) ;
void dsmwritepointerlocalarrayf_(DSMArrayPtr,int64_t*,int64_t*);
void dsmreadpointerlocalarrayf_(DSMArrayPtr,int64_t*,int64_t*);
// TODO! void dsmmemcpypointerlocalarrayf_(DSMArrayPtr, int64_t*, int64_t*, int64_t*);

#ifndef GEN_WrapperSubs
void dsmallocfloatarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmallocfloatf_(DSMScalarPtr,float*) ;
void dsmdeletefloatf_(DSMScalarPtr);
void dsmdeletefloatarrayf_(DSMArrayPtr) ;
void dsmwritefloatarrayf_(DSMArrayPtr,int*,float*);
void dsmreadfloatarrayf_(DSMArrayPtr,int*,float*);

void dsmallocdoublearrayf_(DSMArrayPtr, int64_t* ) ;
void dsmallocdoublef_(double*,double*) ;
void dsmdeletedoublef_(double*);
void dsmdeletedoublearrayf_(DSMArrayPtr) ;
void dsmwritedoublearrayf_(DSMArrayPtr,int*,double*);
void dsmreaddoublearrayf_(DSMArrayPtr,int*,double*);

void dsmallocchararrayf_(DSMArrayPtr, int64_t* ) ;
void dsmalloccharf_(char*,char*) ;
void dsmdeletecharf_(char*);
void dsmdeletechararrayf_(DSMArrayPtr) ;
void dsmwritechararrayf_(DSMArrayPtr,int*,char*);
void dsmreadchararrayf_(DSMArrayPtr,int*,char*);

void dsmallocshortarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmallocshortf_(short*,short*) ;
void dsmdeleteshortf_(short*);
void dsmdeleteshortarrayf_(DSMArrayPtr) ;
void dsmwriteshortarrayf_(DSMArrayPtr,int*,short*);
void dsmreadshortarrayf_(DSMArrayPtr,int*,short*);

void dsmallocintarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmallocintf_(int*,int*) ;
void dsmdeleteintf_(int*);
void dsmdeleteintarrayf_(DSMArrayPtr) ;
void dsmwriteintarrayf_(DSMArrayPtr,int*,int*);
void dsmreadintarrayf_(DSMArrayPtr,int*,int*);

void dsmalloclongarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmalloclongf_(long*,long*) ;
void dsmdeletelongf_(long*);
void dsmdeletelongarrayf_(DSMArrayPtr) ;
void dsmwritelongarrayf_(DSMArrayPtr,int*,long*);
void dsmreadlongarrayf_(DSMArrayPtr,int*,long*);
#endif

#endif // _DSM_ARRAY_C_H_
