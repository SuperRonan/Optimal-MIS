#pragma once

#include "SimpleEstimator.h"

namespace MIS
{
	/// <summary>
	/// Implementation of the Naive MIS estimator
	/// </summary>
	/// <typeparam name="Spectrum"></typeparam>
	/// <typeparam name="Float"></typeparam>
	template <class Spectrum, class Float>
	class NaiveEstimator : public SimpleEstimator<Spectrum, Float>
	{
	public:

		NaiveEstimator(int N):
			SimpleEstimator(N, Heuristic::Naive)
		{}

		NaiveEstimator(NaiveEstimator const&) = default;

		NaiveEstimator(NaiveEstimator &&) = default;

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
		{
			int non_zero_techs = 0;
			for (int i = 0; i < this->m_numtechs; ++i)
				if (balance_weights[i] > 0)
					++non_zero_techs;
			assert(non_zero_techs != 0);
			Float weight = Float(1.0) / Float(non_zero_techs);
			this->m_result += estimate * weight;
		}
	};
}