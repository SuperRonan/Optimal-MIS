#pragma once

#include <SimpleEstimator.h>
#include <cmath>
#include <numeric>

namespace MIS
{
	/// <summary>
	/// Implementation of Veach's Power Heuristic estimator
	/// </summary>
	/// <typeparam name="Spectrum"></typeparam>
	/// <typeparam name="Float"></typeparam>
	template <class Spectrum, class Float>
	class PowerEstimator : public SimpleEstimator<Spectrum, Float>
	{
	protected:
		
		Float m_beta;

	public:

		PowerEstimator(int N, Float beta=2.0):
			SimpleEstimator(N),
			m_beta(beta)
		{}

		PowerEstimator(PowerEstimator const&) = default;

		PowerEstimator(PowerEstimator &&) = default;
		
		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
		{
			Float sum = 0;
			for (int i = 0; i < m_numtechs; ++i)
				sum += std::pow(balance_weights[i], m_beta);
			Float power_weight = std::pow(balance_weights[tech_index], m_beta) / sum;
			m_result += estimate * power_weight;
		}
	};
}