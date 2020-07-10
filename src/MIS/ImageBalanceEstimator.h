#pragma once
#include "ImageSimpleEstimator.h"

namespace MIS
{
    template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageBalanceEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	public:

		ImageBalanceEstimator(int N, int width, int height) :
			ImageSimpleEstimator(N, width, height)
		{}

		ImageBalanceEstimator(ImageBalanceEstimator const& other) = default;

		ImageBalanceEstimator(ImageBalanceEstimator && other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index)const override
		{
			return balance_weights[tech_index];
		}
	};
}