#pragma once

// true -> stores the data of the optimal estimator contiguously in one big buffer, slightly faster
// false -> stores the data of the optimal estimator in different buffers, slightly slower
#define OPTIMIS_ONE_CONTIGUOUS_ARRAY true

// true -> enable debug functions of the direct estimator
#define OPTIMIS_ENABLE_DEBUG false


#ifdef _MSC_VER 
#define MIS_forceinline __forceinline
#endif
#ifdef __GNUG__
#define MIS_forceinline __attribute__((always_inline))
#endif
#ifndef MIS_forceinline
#define MIS_forceinline inline
#endif