#pragma once

namespace MIS
{
	enum class Heuristic : int {
		Balance = 1, 
		Power = 2,
		Naive = 4,
		CutOff = 8,
		Maximum = 16,
		Direct = -1,
		Progressive = -2,
	};

// This can still be inlined
	//inline bool isOptimal(Heuristic h)
	//{
	//	return (int)h < 0;
	//}

// This can't not be inlined
	constexpr inline bool isOptimal(Heuristic h) { 
		return (int)h < 0;
	}

	template <class Float>
	struct EstimatorCreateInfo
	{
		Heuristic heuristic;
		int N;
		// Unfortunately, cannot use a union and default values		
		Float power_beta = 2; 
		Float cutoff_alpha = 0.8; 
		int progressive_step = 4;
	};
}