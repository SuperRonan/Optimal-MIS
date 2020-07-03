#include <iostream>
#include "Spectrum.h"
#include <Estimator.h>
#include <BalanceEstimator.h>
#include <DirectEstimator.h>
#include <PowerEstimator.h>
#include <CutOffEstimator.h>
#include <MaximumEstimator.h>
#include <NaiveEstimator.h>
#include <string>

template <class Stream, class Float, int N>
Stream& operator<<(Stream& stream, MISExample::Spectrum<Float, N> const& spec)
{
    stream << '[';
    for (int i = 0; i < spec.size(); ++i)
    {
        stream << std::to_string(spec[i]);
        if (i < spec.size() - 1)
            stream << ", ";
    }
    stream << "]";
    return stream;
}

template <class Estimator>
void testEstimator()
{
    // Just a function to check that it compiles well
    const int N = 2;
    Estimator estimator(N);
    Estimator::Spectrum_Type estimate;
    Estimator::Float_Type weights[N];

    estimator.setSampleForTechnique(1, 2);

    estimate = 1;
    weights[0] = 0.5;
    weights[1] = 0.7;
    estimator.addEstimate(estimate, weights, 0);

    estimate = 2.0;
    weights[0] = 0.3;
    weights[1] = 0.8;
    estimator.addEstimate(estimate, weights, 1);

    estimator.addOneTechniqueEstimate(estimate * 2, 1);

    std::cout << "Estimator result: " << estimator.solve(1)<<std::endl;
    estimator.reset();
}

template <class Spectrum, class Float = double>
void testEstimators()
{
    testEstimator<MIS::BalanceEstimator<Spectrum, Float>>();
    testEstimator<MIS::NaiveEstimator<Spectrum, Float>>();
    testEstimator<MIS::PowerEstimator<Spectrum, Float>>();
    testEstimator<MIS::CutOffEstimator<Spectrum, Float>>();
    testEstimator<MIS::MaximumEstimator<Spectrum, Float>>();
    testEstimator<MIS::DirectEstimator<Spectrum, Float>>();
}

int main(int argc, char ** argv)
{
    using Float = double;
    using RGBColor = MISExample::Spectrum<Float, 3>;

    testEstimators<Float, Float>();
    testEstimators<RGBColor, Float>();

}