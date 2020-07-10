#pragma once

#include <ImageSimpleEstimator.h>

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageMaximumEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	public:

		ImageMaximumEstimator(int N, int width, int height) :
			ImageSimpleEstimator(N, width, height)
		{}

		ImageMaximumEstimator(ImageMaximumEstimator const& other) = default;

		ImageMaximumEstimator(ImageMaximumEstimator&& other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index) const override
		{
			Float const& wt = balance_weights[tech_index];
			if (std::all_of(balance_weights, balance_weights + m_numtechs, [wt](Float wi) {return wt > wi; }))
				return Float(1.0);
			return Float(0);
		}
	};
}