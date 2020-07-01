#include <iostream>
#include "Spectrum.h"
#include <Estimator.h>
#include <BalanceEstimator.h>
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
    using RGBColor = MISExample::Spectrum<double, 3>;
    MIS::BalanceEstimator<RGBColor, double> estimator(1);
    estimator.addEstimate(1.0, nullptr, 0);
    std::cout<<"Hello World!"<<std::endl;
    std::cout << "Estimator result: " << estimator.solve(1) << std::endl;
}