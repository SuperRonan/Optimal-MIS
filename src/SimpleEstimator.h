#pragma once

#include <Estimator.h>

namespace MIS
{
    /// <summary>
    /// Partial implementation for simple MIS Heuristics.
    /// </summary>
    /// <typeparam name="Spectrum"></typeparam>
    /// <typeparam name="Float"></typeparam>
    template <class Spectrum, class Float>
    class SimpleEstimator : public Estimator<Spectrum, Float>
    {
    protected:

        /// <summary>
        /// The result of the estimator is accumulated in this variable.
        /// </summary>
        Spectrum m_result;

    public:

        SimpleEstimator(int N):
            Estimator(N)
        {}

        SimpleEstimator(SimpleEstimator const&) = default;

        SimpleEstimator(SimpleEstimator &&) = default;

        virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index) override
        {
            m_result += estimate;
        }

        virtual Spectrum solve(int iterations) override
        {
            return m_result / Float(iterations);
        }

        virtual void reset() override
        {
            m_result = 0;
        }
    };
}