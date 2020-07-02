#include <iostream>
#include "Spectrum.h"
#include <Estimator.h>
#include <BalanceEstimator.h>
#include <DirectEstimator.h>
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

int main(int argc, char ** argv)
{
    using Float = double;
    using RGBColor = MISExample::Spectrum<Float, 3>;
    //MIS::BalanceEstimator<RGBColor, double> estimator(1);
    MIS::DirectEstimator<RGBColor, Float> estimator(1);
    Float weights[] = { 1.0 };
    estimator.addEstimate(1.0, weights, 0);
    std::cout << "Estimator result: " << estimator.solve(1) << std::endl;
}