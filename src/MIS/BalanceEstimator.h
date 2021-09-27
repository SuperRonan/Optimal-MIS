#pragma once

#include "SimpleEstimator.h"

namespace MIS
{
    /// <summary>
    /// Implementation of Veach's Balance Heuristic estimator
    /// </summary>
    /// <typeparam name="Spectrum"></typeparam>
    /// <typeparam name="Float"></typeparam>
    template <class Spectrum, class Float = double>
    class BalanceEstimator: public SimpleEstimator<Spectrum, Float>
    {
    public:

        BalanceEstimator(int N):
            SimpleEstimator(N, Heuristic::Balance)
        {}

        BalanceEstimator(BalanceEstimator const&) = default;

        BalanceEstimator(BalanceEstimator &&) = default;

        virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
        {
            this->m_result += estimate * balance_weights[tech_index];
        }
    };
}