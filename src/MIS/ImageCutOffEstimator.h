#pragma once

#include <ImageSimpleEstimator.h>

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageCutOffEstimator : public ImageSimpleEstimator<Spectrum, Float, ROW_MAJOR>
	{
	protected:

		Float m_alpha;

	public:

		ImageCutOffEstimator(int N, int width, int height, Float alpha=0.5):
			ImageSimpleEstimator(N, width, height),
			m_alpha(alpha)
		{}

		ImageCutOffEstimator(ImageCutOffEstimator const& other) = default;

		ImageCutOffEstimator(ImageCutOffEstimator&& other) = default;

		virtual Float weight(const Float* balance_weights, int tech_index) const override
		{
			int best_index = 0;
			Float max_w = balance_weights[0];
			for (int i = 1; i < m_numtechs; ++i)
				if (balance_weights[i] > max_w)
				{
					best_index = i;
					max_w = balance_weights[i];
				}
			if (balance_weights[tech_index] >= m_alpha * max_w)
			{
				Float sum = 0;
				for (int i = 0; i < m_numtechs; ++i)
					if (balance_weights[i] >= m_alpha * max_w)
						sum += balance_weights[i];
				return balance_weights[tech_index] / sum;
			}
			return Float(0);
		}
	};
}