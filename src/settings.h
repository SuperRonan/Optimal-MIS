#pragma once

// true -> stores the data of the optimal estimator contiguously in one big buffer, slightly faster
// false -> stores the data of the optimal estimator in different buffers, slightly slower
#define OPTIMIS_ONE_CONTIGUOUS_ARRAY true

// true -> enable debug functions of the direct estimator
#define OPTIMIS_ENABLE_DEBUG false