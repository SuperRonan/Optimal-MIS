#pragma once

#include <SimpleEstimator.h>

namespace MIS
{
    template <class Spectrum, class Float = double>
    class BalanceEstimator: public SimpleEstimator<Spectrum, Float>
    {
    public:

        BalanceEstimator(int N):
            SimpleEstimator(N)
        {}

        BalanceEstimator(BalanceEstimator const& other) = default;

        virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
        {
            m_result += estimate * balance_weights[tech_index];
        }

        virtual void addOneTechniqueEstimate(Spectrum const& estimate, int) override
        {
            m_result += estimate;
        }
    };
}