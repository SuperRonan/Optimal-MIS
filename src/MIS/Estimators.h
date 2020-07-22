#pragma once

#include "Estimator.h"
#include "BalanceEstimator.h"
#include "DirectEstimator.h"
#include "PowerEstimator.h"
#include "CutOffEstimator.h"
#include "MaximumEstimator.h"
#include "NaiveEstimator.h"

#include "Heuristics.h"

#include <exception>

namespace MIS
{
	template <class Spectrum, class Float>
	Estimator<Spectrum, Float>* createEstimator(Heuristic h, int N)
	{
		Estimator<Spectrum, Float>* res = nullptr;
		switch (h)
		{
			case Heuristic::Balance:
				res = new BalanceEstimator<Spectrum, Float>(N);
			break;
			case Heuristic::Power:
				res = new PowerEstimator<Spectrum, Float>(N);
			break;
			case Heuristic::Naive:
				res = new NaiveEstimator<Spectrum, Float>(N);
			break;
			case Heuristic::CutOff:
				res = new CutOffEstimator<Spectrum, Float>(N);
			break;
			case Heuristic::Maximum:
				res = new MaximumEstimator<Spectrum, Float>(N);
			break;
			case Heuristic::Direct:
				res = new DirectEstimator<Spectrum, Float>(N);
			break;
			default:
				throw std::logic_error("Heuristic not recognized! Cannot create the estimator");
			break;
		}
		return res;
	}
}