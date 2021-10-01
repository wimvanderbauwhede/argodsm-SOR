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

#pragma GEN WrapperSubs

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
