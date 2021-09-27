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
#define isOptimal(h) ((int h) < 0)
}