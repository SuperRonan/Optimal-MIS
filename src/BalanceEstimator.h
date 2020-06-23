#pragma once

#include <Estimator.h>

namespace MIS
{
template <class Spectrum, class Float = double>
class BalanceEstimator: public Estimator<Spectrum, Float>
{
protected:
    Spectrum m_res = 0;
public:

    BalanceEstimator(int N):
        Estimator(N)
    {}

    BalanceEstimator(BalanceEstimator const& other) = default;

    virtual void addEstimate(Spectrum const& balance_estimate, const Float*, int) override
    {
        m_res += balance_estimate;
    }

    virtual void addOneTechniqueEstimate(Spectrum const& estimate, int) override
    {
        m_res += estimate;
    }

    virtual Spectrum solve(int iterations) override
    {
        return m_res / (Float)iterations;
    }

    virtual void reset() override
    {
        m_res = 0;
    }
};
}