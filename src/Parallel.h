#pragma once

#include <omp.h>
#include <vector>



class Parallel
{
protected:


public:

	static void init()
	{}

	static void setNumThreads(int n)
	{
		omp_set_num_threads(std::max(1, std::min(getCPUThreads(), n)));
	}

	static int getNumThreads()
	{
		return omp_get_max_threads();
	}

	static int getCPUThreads()
	{
		return omp_get_num_procs();
	}

	static int getCurrentNumThreads()
	{
		return omp_get_num_threads();
	}

	template <class T>
	__forceinline static std::vector<T> preAllocate(T const& t=T())
	{
		return std::vector<T>(getNumThreads(), t);
	}

	static int tid()
	{
		return omp_get_thread_num();
	}

	template <class Function>
	__forceinline static void ParallelFor(int min, int max, Function const& function)
	{
#pragma omp parallel for schedule(dynamic)
			for (int i = min; i < max; ++i)
			{
				function(i);
			}
	}
};