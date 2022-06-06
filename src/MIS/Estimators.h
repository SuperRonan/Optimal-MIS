#pragma once

#include "Estimator.h"
#include "BalanceEstimator.h"
#include "DirectEstimator.h"
#include "PowerEstimator.h"
#include "CutOffEstimator.h"
#include "MaximumEstimator.h"
#include "NaiveEstimator.h"
#include "ProgressiveEstimator.h"
#include "AlphaEstimator.h"

#include "Heuristics.h"

#include <exception>

namespace MIS
{
	/// <summary>
	/// Integer can be signed or unsigned, it does not matter.
	/// It is here to precise the width of the integers to use (we recomand 64 bits)
	/// <returns></returns>
	template <class Spectrum, class Float, class Integer=int64_t>
	Estimator<Spectrum, Float>* createEstimator(EstimatorCreateInfo<Float> const& estimator_create_info=EstimatorCreateInfo<Float>())
	{
		Estimator<Spectrum, Float>* res = nullptr;
		switch (estimator_create_info.heuristic)
		{
			case Heuristic::Balance:
				res = new BalanceEstimator<Spectrum, Float>(estimator_create_info.N);
			break;
			case Heuristic::Power:
				res = new PowerEstimator<Spectrum, Float>(estimator_create_info.N, estimator_create_info.power_beta);
			break;
			case Heuristic::Naive:
				res = new NaiveEstimator<Spectrum, Float>(estimator_create_info.N);
			break;
			case Heuristic::CutOff:
				res = new CutOffEstimator<Spectrum, Float>(estimator_create_info.N, estimator_create_info.cutoff_alpha);
			break;
			case Heuristic::Maximum:
				res = new MaximumEstimator<Spectrum, Float>(estimator_create_info.N);
			break;
			case Heuristic::Direct:
				res = new DirectEstimator<Spectrum, Float, Integer>(estimator_create_info.N);
			break;
			case Heuristic::Progressive:
				res = new ProgressiveEstimator<Spectrum, Float, Integer>(estimator_create_info.N, estimator_create_info.progressive_step);
			break;
			case Heuristic::Alpha:
				res = new AlphaEstimator<Spectrum, Float>(estimator_create_info.N);
			break;
			default:
				throw std::logic_error("Heuristic not recognized! Cannot create the estimator");
			break;
		}
		return res;
	}
}