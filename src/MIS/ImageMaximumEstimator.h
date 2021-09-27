#pragma once

#include "ImageSimpleEstimator.h"

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageMaximumEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	public:

		ImageMaximumEstimator(int N, int width, int height) :
			ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>::ImageSimpleEstimator(N, width, height, Heuristic::Maximum)
		{}

		ImageMaximumEstimator(ImageMaximumEstimator const& other) = default;

		ImageMaximumEstimator(ImageMaximumEstimator&& other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index) const override
		{
			Float const& wt = balance_weights[tech_index];
			for (int i = 0; i < this->m_numtechs; ++i)
				if (balance_weights[i] > wt)
					return Float(0);
			return Float(1);
		}
	};
}