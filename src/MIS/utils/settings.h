#pragma once

// true -> stores the data of the optimal estimator contiguously in one big buffer, slightly faster
// false -> stores the data of the optimal estimator in different buffers, slightly slower
#ifndef OPTIMIS_ONE_CONTIGUOUS_ARRAY
	#define OPTIMIS_ONE_CONTIGUOUS_ARRAY true
#endif

// 1 -> old version
// 2 -> new version, probably faster.
#ifndef OPTIMIS_IMPL
#define OPTIMIS_IMPL 2
#endif

#ifdef _MSC_VER 
#define MIS_forceinline __forceinline
#endif
#ifdef __GNUG__
#define MIS_forceinline __attribute__((always_inline))
#endif
#ifndef MIS_forceinline
#define MIS_forceinline inline
#endif