
#include <iostream>
#include "Spectrum.h"

#include <MIS/Estimators.h>
		  
#include <MIS/ImageEstimators.h>

#include <string>

#include <random>
#include <memory>


const std::vector<MIS::Heuristic> heuristics = { MIS::Heuristic::Balance,
		MIS::Heuristic::CutOff, MIS::Heuristic::Maximum, MIS::Heuristic::Power,
		MIS::Heuristic::Naive, MIS::Heuristic::Direct, MIS::Heuristic::Progressive };

const std::vector<MIS::Heuristic> image_heuristics = { MIS::Heuristic::Balance,
		MIS::Heuristic::CutOff, MIS::Heuristic::Maximum, MIS::Heuristic::Power,
		MIS::Heuristic::Naive, MIS::Heuristic::Direct};

template <class Estimator>
void testEstimator()
{
	const auto printSystem = [&](MIS::LinearSystem< Estimator::Float_Type> const& system)
	{
		std::cout << "matrix: \n";
		std::cout << system.tech_matrix << "\n";
		std::cout << "vectors: \n";
		std::cout << system.contrib_vectors << "\n";
		std::cout << "alphas: \n";
		std::cout << system.alphas << "\n";
	};
	// Just a function to check that it compiles well
	const int N = 2;
	Estimator _estimator(N);
	Estimator estimator = _estimator;
	typename Estimator::Spectrum_Type estimate;
	typename Estimator::Float_Type weights[N];

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

	std::cout << "Estimator " << (int)estimator.m_heuristic << " result: " << estimator.solve(1) << std::endl;

	if constexpr (std::is_same<Estimator, MIS::DirectEstimator< Estimator::Spectrum_Type, Estimator::Float_Type>>::value)
	{
		MIS::LinearSystem< Estimator::Float_Type> system = estimator.getLinearSystem(1);
		printSystem(system);
	}

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
	testEstimator<MIS::ProgressiveEstimator<Spectrum, Float>>();
}

template <class Spectrum, class Float>
void testVirtualEstimators()
{
	const int N = 2;
	
	for (MIS::Heuristic h : heuristics)
	{
		MIS::EstimatorCreateInfo<Float> eci;
		eci.heuristic = h;
		eci.N = 2;
		MIS::Estimator<Spectrum, Float>* estimator = MIS::createEstimator<Spectrum, Float>(eci);
		Spectrum estimate;
		Float weights[N];
		
		estimator->setSampleForTechnique(1, 2);

		estimate = 1;
		weights[0] = 0.5;
		weights[1] = 0.7;
		estimator->addEstimate(estimate, weights, 0);

		estimate = 2.0;
		weights[0] = 0.3;
		weights[1] = 0.8;
		estimator->addEstimate(estimate, weights, 1);

		estimator->addOneTechniqueEstimate(estimate * 2, 1);

		std::cout << "Estimator " << (int)estimator->m_heuristic << " result: " << estimator->solve(1) << std::endl;
		estimator->reset();
	}   
}

template <class ImageEstimator>
void testImageEstimator()
{
	using Float = typename ImageEstimator::Float_Type;
	using Spectrum = typename ImageEstimator::Spectrum_Type;
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
	testImageEstimator<MIS::ImageDirectEstimator<Spectrum, Float, size_t, true>>();
	testImageEstimator<MIS::ImagePowerEstimator<Spectrum, Float, true>>();
	testImageEstimator<MIS::ImageNaiveEstimator<Spectrum, Float, true>>();
	testImageEstimator<MIS::ImageCutOffEstimator<Spectrum, Float, true>>();
	testImageEstimator<MIS::ImageMaximumEstimator<Spectrum, Float, true>>();
}

template <class Spectrum, class Float = double>
void testVirtualImageEstimators()
{
	const int N = 2;
	int width = 10;
	int height = 10;

	for (MIS::Heuristic h : image_heuristics)
	{
		MIS::EstimatorCreateInfo<Float> eci;
		eci.heuristic = h;
		eci.N = N;
		MIS::ImageEstimator<Spectrum, Float, true>* estimator = MIS::createImageEstimator<Spectrum, Float, true>(width, height, eci);
		estimator->setSampleForTechnique(1, 2);

		for (int i = 0; i < estimator->width(); ++i)  for (int j = 0; j < estimator->height(); ++j)
		{
			Float u = Float(i) / Float(width);
			Float v = Float(j) / Float(height);
			Spectrum estimate = u + v;

			Float weights[N];
			weights[0] = u + v;
			weights[1] = 0.7;
			estimator->addEstimate(estimate, weights, 0, u, v);

			estimate = u + v;
			weights[0] = u + v;
			weights[1] = 0.8;
			estimator->addEstimate(estimate, weights, 1, u, v);

			estimator->addOneTechniqueEstimate(estimate * 2, 1, u, v);
		}
		estimator->loop();
		std::vector<Spectrum> res(width * height, 0);
		estimator->solve(res.data(), 1);

		if (h == MIS::Heuristic::Direct)
		{
			MIS::ImageDirectEstimator<Spectrum, Float, size_t, true>* direct_estimator = (MIS::ImageDirectEstimator<Spectrum, Float, size_t, true> *) estimator;
			MIS::LinearSystem<Float> system = direct_estimator->getPixelLinearSystem(1, 5, 5);
		}

		delete estimator;
	}
}

int main(int argc, char ** argv)
{
	using Float = double;
	using RGBColor = MISExample::Spectrum<Float, 3>;

	testEstimators<Float, Float>();
	testEstimators<RGBColor, Float>();

	testVirtualEstimators<Float, Float>();
	testVirtualEstimators<RGBColor, Float>();

	testImageEstimators<Float, Float>();
	testImageEstimators<RGBColor, Float>();

	testVirtualImageEstimators<Float, Float>();
	testVirtualImageEstimators<RGBColor, Float>();

}