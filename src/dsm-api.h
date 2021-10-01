/*
   This is a function-based ArgoDSM API intended for use in Fortran programs

   The approach is to cast the pointer to the array to a 64-bit integer and return it.
   This word gets passed around in the Fortran code (as INTEGER*8).
   Every call to the API takes this word as its first argument, casts it back to the object and so on.

   So we create some low-level casting functions first
*/
#ifndef _DSM_ARRAY_F_H_
#define _DSM_ARRAY_F_H_


// #include <cassert>
// #include <limits>
// #include <iostream>
// #include <vector>

// #include <string>
// #include <cctype>
// #include <algorithm>

#include "argo/argo.hpp"

template<typename TPtr> TPtr fromWord(int64_t ivp) {
	int64_t* ip=(int64_t*)ivp;
	void* vp=(void*)ip;
	TPtr tp = (TPtr)vp;
	return tp;
}

template<typename TPtr> int64_t toWord(TPtr tp) {
	void* vp = reinterpret_cast<void*>(tp);
	int64_t ivp = (int64_t)vp;
	return ivp;
}

struct DSM3DArray {
	int64_t isz,jsz,ksz; // Size including the halos
	int64_t ioff = 0,joff = 0,koff = 0; // Offset (0 means core starts at 1)
	int ilh = 0,jlh = 0,klh = 0;
	int ihh = 0,jhh = 0,khh = 0;
	int64_t* coll_ptr = 0; // Pointer to the collective (argo::conew_array) array of argo::new_array<int64_t> pointer
	int64_t** ptr = 0; // Pointer to a local (new int64_t[]) array with those argo::new_array<int64_t> pointers copied from the collective one
	bool shared = true; // Collective 
};
#ifdef NEW_DSM3
struct DSM3DArrayH {
	int64_t isz,jsz,ksz; // Size including the halos
	int64_t ioff = 0,joff = 0,koff = 0; // Offset (0 means core starts at 1)
	int ilh = 0,jlh = 0,klh = 0;
	int ihh = 0,jhh = 0,khh = 0;
	int64_t* core_ptr = 0; // Pointer to the core
	int64_t* coll_ptr = 0; // Pointer to the collective (argo::conew_array) array of argo::new_array<int64_t> pointer
	int64_t* ptr = 0; // Pointer to a local (new int64_t[]) array with those argo::new_array<int64_t> pointers copied from the collective one
	bool shared = true; // Collective 
};
#else
struct DSM3DArrayH {
	int64_t isz,jsz,ksz; // Size including the halos
	int64_t ioff = 0,joff = 0,koff = 0; // Offset (0 means core starts at 1)
	int ilh = 0,jlh = 0,klh = 0;
	int ihh = 0,jhh = 0,khh = 0;
	int64_t* coll_ptr = 0; // Pointer to the collective (argo::conew_array) array of argo::new_array<int64_t> pointer	
	int64_t** ptr = 0; // Pointer to a local (new int64_t[]) array with those argo::new_array<int64_t> pointers copied from the collective one
	bool shared = true; // Collective 
};
#endif

extern "C" {
#include "dsm-C.h"
}
#endif // _DSM_ARRAY_F_H_
