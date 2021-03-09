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
	template <class Spectrum, class Float, bool ROW_MAJOR, class Uint = size_t>
	ImageEstimator<Spectrum, Float, ROW_MAJOR>* createImageEstimator(Heuristic h, int N, int width, int height)
	{
		ImageEstimator<Spectrum, Float, ROW_MAJOR>* res = nullptr;
		switch (h)
		{
		case Heuristic::Balance:
			res = new ImageBalanceEstimator<Spectrum, Float, ROW_MAJOR>(N, width, height);
			break;
		case Heuristic::Power:
			res = new ImagePowerEstimator<Spectrum, Float, ROW_MAJOR>(N, width, height);
			break;
		case Heuristic::Naive:
			res = new ImageNaiveEstimator<Spectrum, Float, ROW_MAJOR>(N, width, height);
			break;
		case Heuristic::CutOff:
			res = new ImageCutOffEstimator<Spectrum, Float, ROW_MAJOR>(N, width, height);
			break;
		case Heuristic::Maximum:
			res = new ImageMaximumEstimator<Spectrum, Float, ROW_MAJOR>(N, width, height);
			break;
		case Heuristic::Direct:
			res = new ImageDirectEstimator<Spectrum, Float, Uint, ROW_MAJOR>(N, width, height);
			break;
		default:
			throw std::logic_error("Heuristic not recognized! Cannot create the estimator");
			break;
		}
		return res;
	}
}