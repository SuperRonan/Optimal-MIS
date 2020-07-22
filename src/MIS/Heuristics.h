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
		//Progressive = -2,
	};

	__forceinline bool iseOptimal(Heuristic h)
	{
		return (int)h < 0;
	}
}