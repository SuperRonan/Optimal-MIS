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
			SimpleEstimator(N, Heuristic::Maximum)
		{}

		MaximumEstimator(MaximumEstimator const&) = default;

		MaximumEstimator(MaximumEstimator &&) = default;

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
		{
			Float const& wt = balance_weights[tech_index];
			if (std::all_of(balance_weights, balance_weights + m_numtechs, [wt](Float wi) {return wt > wi; }))
				m_result += estimate;
		}
	};
}