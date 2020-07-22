#pragma once

#include "ImageSimpleEstimator.h"

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageNaiveEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	public:

		ImageNaiveEstimator(int N, int width, int height) :
			ImageSimpleEstimator(N, width, height, Heuristic::Naive)
		{}

		ImageNaiveEstimator(ImageNaiveEstimator const& other) = default;

		ImageNaiveEstimator(ImageNaiveEstimator&& other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index) const override
		{
			int non_zero_techs = 0;
			for (int i = 0; i < m_numtechs; ++i)
				if (balance_weights[i] > 0)
					++non_zero_techs;
			return Float(1.0) / Float(non_zero_techs);
		}
	};
}