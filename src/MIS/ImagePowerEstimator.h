#pragma once

#include "ImageSimpleEstimator.h"

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImagePowerEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	protected:

		Float m_beta;

	public:

		ImagePowerEstimator(int N, int width, int height, Float beta=2.0) :
			ImageSimpleEstimator(N, width, height),
			m_beta(beta)
		{}

		ImagePowerEstimator(ImagePowerEstimator const& other) = default;

		ImagePowerEstimator(ImagePowerEstimator&& other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index) const override
		{
			Float sum = 0;
			for (int i = 0; i < m_numtechs; ++i)
				sum += std::pow(balance_weights[i], m_beta);
			Float power_weight = std::pow(balance_weights[tech_index], m_beta) / sum;
			return power_weight;
		}
	};
}