#pragma once

#include "ImageEstimator.h"
#include "ImageBalanceEstimator.h"
#include "ImageDirectEstimator.h"
#include "ImageCutOffEstimator.h"
#include "ImageMaximumEstimator.h"
#include "ImagePowerEstimator.h"
#include "ImageNaiveEstimator.h"

#include "Heuristics.h"

#include <exception>

namespace MIS
{
	/// <summary>
	/// Integer can be signed or unsigned, it is here to precise the width of the integers to use (we recomand 64 bits)
	/// <returns></returns>
	template <class Spectrum, class Float, bool ROW_MAJOR, class Integer = int64_t>
	ImageEstimator<Spectrum, Float, ROW_MAJOR>* createImageEstimator(int width, int height, EstimatorCreateInfo<Float> const& estimator_create_info = EstimatorCreateInfo<Float>())
	{
		ImageEstimator<Spectrum, Float, ROW_MAJOR>* res = nullptr;
		switch (estimator_create_info.heuristic)
		{
		case Heuristic::Balance:
			res = new ImageBalanceEstimator<Spectrum, Float, ROW_MAJOR>(estimator_create_info.N, width, height);
			break;
		case Heuristic::Power:
			res = new ImagePowerEstimator<Spectrum, Float, ROW_MAJOR>(estimator_create_info.N, width, height, estimator_create_info.power_beta);
			break;
		case Heuristic::Naive:
			res = new ImageNaiveEstimator<Spectrum, Float, ROW_MAJOR>(estimator_create_info.N, width, height);
			break;
		case Heuristic::CutOff:
			res = new ImageCutOffEstimator<Spectrum, Float, ROW_MAJOR>(estimator_create_info.N, width, height, estimator_create_info.cutoff_alpha);
			break;
		case Heuristic::Maximum:
			res = new ImageMaximumEstimator<Spectrum, Float, ROW_MAJOR>(estimator_create_info.N, width, height);
			break;
		case Heuristic::Direct:
			res = new ImageDirectEstimator<Spectrum, Float, Integer, ROW_MAJOR>(estimator_create_info.N, width, height);
			break;
		default:
			throw std::logic_error("Heuristic not recognized! Cannot create the estimator");
			break;
		}
		return res;
	}
}