#pragma once

#include "SimpleEstimator.h"

namespace MIS
{
	/// <summary>
	/// Implementation of Veach's Cut Off Heuristic estimator
	/// </summary>
	/// <typeparam name="Spectrum"></typeparam>
	/// <typeparam name="Float"></typeparam>
	template <class Spectrum, class Float>
	class CutOffEstimator : public SimpleEstimator<Spectrum, Float>
	{
	protected:

		Float m_alpha;

	public:

		CutOffEstimator(int N, Float alpha=0.5):
			SimpleEstimator(N),
			m_alpha(alpha)
		{}

		CutOffEstimator(CutOffEstimator const&) = default;

		CutOffEstimator(CutOffEstimator &&) = default;

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
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
				Float weight = balance_weights[tech_index] / sum;
				m_result += estimate * weight;
			}
			// else weight is 0
		}

	};
}