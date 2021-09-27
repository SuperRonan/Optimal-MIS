#pragma once

#include <omp.h>
#include <vector>
#include "settings.h"



namespace MIS
{
	/// <summary>
	/// Static class for parallelization, based on open MP
	/// </summary>
	class Parallel
	{
	protected:
	public:

		static void init()
		{}

		/// <summary>
		/// Sets the number of threads. If n is greater than the number of threads that can be provided, it will be clamped
		/// </summary>
		static void setNumThreads(int n)
		{
			omp_set_num_threads(std::max(1, std::min(getCPUThreads(), n)));
		}

		/// <summary>
		/// Returns the current number of threads setted by "setNumThreads"
		/// </summary>
		static int getNumThreads()
		{
			return omp_get_max_threads();
		}

		/// <summary>
		/// Returns the maximum number of (physical) threads of the CPU
		/// </summary>
		static int getCPUThreads()
		{
			return omp_get_num_procs();
		}

		/// <summary>
		/// Returns the number of OMP threads running at the time of calling this function
		/// </summary>
		/// <returns></returns>
		static int getCurrentNumThreads()
		{
			return omp_get_num_threads();
		}

		/// <summary>
		/// Pre allocate a vector of T the size of the wished number of threads
		/// </summary>
		template <class T>
		MIS_forceinline static std::vector<T> preAllocate(T const& t=T())
		{
			return std::vector<T>(getNumThreads(), t);
		}

		/// <summary>
		/// Return the OMP thread id 
		/// </summary>
		static int tid()
		{
			return omp_get_thread_num();
		}

		/// <summary>
		/// Execute in parallel of for loop from min (included) to max (excluded), and calls function at each iteration
		/// </summary>
		/// <param name="min">: included</param>
		/// <param name="max">: excluded</param>
		/// <param name="function">: A function taking an int between min and max</param>
		template <class Function>
		MIS_forceinline static void parallelFor(int min, int max, Function const& function)
		{
	#pragma omp parallel for schedule(dynamic)
				for (int i = min; i < max; ++i)
				{
					function(i);
				}
		}
	};
}