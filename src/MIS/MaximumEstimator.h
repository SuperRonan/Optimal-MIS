#pragma once

#include "SimpleEstimator.h"
#include <algorithm>

namespace MIS
{
	/// <summary>
	/// Implementation of Veach's Maximum Heuristic estimator
	/// </summary>
	/// <typeparam name="Spectrum"></typeparam>
	/// <typeparam name="Float"></typeparam>
	template <class Spectrum, class Float>
	class MaximumEstimator : public SimpleEstimator<Spectrum, Float>
	{
	public:

		MaximumEstimator(int N):
			SimpleEstimator<Spectrum, Float>::SimpleEstimator(N, Heuristic::Maximum)
		{}

		MaximumEstimator(MaximumEstimator const&) = default;

		MaximumEstimator(MaximumEstimator &&) = default;

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
		{
			Float const& wt = balance_weights[tech_index];
			for (int i = 0; i < this->m_numtechs; ++i)
				if (balance_weights[i] > wt)
					return;
			this->m_result += estimate;
		}
	};
}