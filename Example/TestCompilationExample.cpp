#include <iostream>
#include "Spectrum.h"

#include <MIS/Estimators.h>
          
#include <MIS/ImageEstimators.h>

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
    Estimator _estimator(N);
    Estimator estimator = _estimator;
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

template <class ImageEstimator>
void testImageEstimator()
{
    using Float = ImageEstimator::Float_Type;
    using Spectrum = ImageEstimator::Spectrum_Type;
    const int N = 2;
    int w = 10;
    int h = 10;
    
    ImageEstimator _estimator(N, w, h);
    ImageEstimator estimator = _estimator;
    estimator.setSampleForTechnique(1, 2);
    
    for (int i = 0; i < w; ++i)  for (int j = 0; j < h; ++j)
    {
        Float u = Float(i) / Float(w);
        Float v = Float(j) / Float(h);
        Spectrum estimate = u + v;

        Float weights[N];
        weights[0] = u + v;
        weights[1] = 0.7;
        estimator.addEstimate(estimate, weights, 0, u, v);

        estimate = u + v;
        weights[0] = u + v;
        weights[1] = 0.8;
        estimator.addEstimate(estimate, weights, 1, u, v);

        estimator.addOneTechniqueEstimate(estimate * 2, 1, u, v);
    }
    estimator.loop();
    std::vector<Spectrum> res(w * h, 0);
    estimator.solve(res.data(), 1);


}

template <class Spectrum, class Float = double>
void testImageEstimators()
{
    testImageEstimator<MIS::ImageBalanceEstimator<Spectrum, Float, true>>();
    testImageEstimator<MIS::ImageDirectEstimator<Spectrum, Float, true>>();
    testImageEstimator<MIS::ImagePowerEstimator<Spectrum, Float, true>>();
    testImageEstimator<MIS::ImageNaiveEstimator<Spectrum, Float, true>>();
    testImageEstimator<MIS::ImageCutOffEstimator<Spectrum, Float, true>>();
    testImageEstimator<MIS::ImageMaximumEstimator<Spectrum, Float, true>>();
}

int main(int argc, char ** argv)
{
    using Float = double;
    using RGBColor = MISExample::Spectrum<Float, 3>;

    testEstimators<Float, Float>();
    testEstimators<RGBColor, Float>();

    testImageEstimators<Float, Float>();
    testImageEstimators<RGBColor, Float>();
}