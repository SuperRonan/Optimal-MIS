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

	std::cout << "Estimator " << (int)estimator.m_heuristic << " result: " << estimator.solve(1) << std::endl;
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
		MIS::Estimator<Spectrum, Float>* estimator = MIS::createEstimator<Spectrum, Float>(h, N);
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
		MIS::ImageEstimator<Spectrum, Float, true>* estimator = MIS::createImageEstimator<Spectrum, Float, true>(h, N, width, height);
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
	}
}

template <class Float, class Estimator>
void testConservativeHypothesis(bool conservative, bool use_pdf=false, bool verbose=false)
{
	struct Sampler
	{
		std::mt19937_64 m_rng;

		Sampler(std::mt19937_64::result_type seed = std::mt19937_64::default_seed) : m_rng(seed) {};

		Float uniform(Float min = 0, Float max = 1)
		{
			Float u = Float(m_rng()) / Float(m_rng.max());
			return u * (max - min) + min;
		};
	};

	struct Sample 
	{
		Float x;
		Float p;
	};

	struct Tech
	{
		int n = 1;
		Tech(int n = 1) :n(n) {};
		virtual Sample sample(Sampler& sampler) const = 0;
		virtual Float pdf(Float x) const = 0;
	};

	struct UniTech : public Tech
	{
		Float min, max;
		UniTech(int n=1, Float min = 0, Float max = 1) : Tech(n), min(min), max(max) {};

		virtual Sample sample(Sampler& sampler) const override
		{
			Sample res;
			res.x = sampler.uniform(min, max);
			res.p = pdf(res.x);
			return res;
		};

		virtual Float pdf(Float x) const override
		{
			if (x >= min && x <= max)	return 1.0 / std::abs(max - min);
			else						return 0;
		};
	};

	std::vector<std::shared_ptr<Tech>> techs = {
		std::make_shared<UniTech>(1, 0, 3), // Worst
		std::make_shared<UniTech>(1, 0, 2), // half bad
		std::make_shared<UniTech>(1, 1, 3), // half bad
		std::make_shared<UniTech>(1, 1, 2), // good 
		//std::make_shared<UniTech>(1, 1.25, 1.75), // optimal
	};

	std::vector<Float> Wb(techs.size());
	
	
	size_t iterations = 64;

	const auto f = [](Float x) {
		if (x > 1.25 && x < 1.75)	return Float(2);
		else						return Float(0);
	};
	// Integral of f
	const Float F = 1;

	Estimator estimator(techs.size());
	for (int i = 0; i < techs.size(); ++i)	estimator.setSampleForTechnique(i, techs[i]->n);

	std::cout << "Ground Truth: " << F << std::endl;
	Sampler sampler;
	Float avg_error = 0, avg_abs_error=0;
	int L = 16*16*16;
	for (int n = 0; n < L; ++n)
	{
		estimator.reset();

		for (size_t iter = 0; iter < iterations; ++iter)
		{
			for (int i = 0; i < techs.size(); ++i)
			{
				const Tech& ti = *techs[i];
				for (int j = 0; j < ti.n; ++j)
				{
					Sample s = ti.sample(sampler);
					Float fx = f(s.x) / (ti.n * s.p);
					Float sum = 0;
					if (fx == 0 && conservative)
						continue;
					for (int l = 0; l < techs.size(); ++l)
					{
						Float ql = techs[l]->pdf(s.x) * (use_pdf ? 1 : techs[l]->n);
						sum += ql;
						Wb[l] = ql;
					}
					for (int l = 0; l < techs.size(); ++l)	Wb[l] /= sum;
					estimator.addEstimate(fx, Wb.data(), i);
				}
			}
			estimator.loop();
		}

		Float estimation = estimator.solve(iterations);
		Float error = F - estimation;
		if(verbose)
			std::cout << "Conservative estimation: " << estimation <<", Error: " << error << std::endl;
		avg_abs_error += std::abs(error);
		avg_error += error;

		if (n == 0)
		{
			if constexpr (std::is_same<Estimator, MIS::DirectEstimator<Float, Float>>::value)
			{
				MIS::DirectEstimator<Float, Float>& d = estimator;
				MIS::DirectEstimator<Float, Float>::LinearSystem system = d.getLinearSystem(iterations);
				for (int i = 0; i < techs.size(); ++i)	system.alpha[i] *= techs[i]->n;
				std::cout << "A: \n" << system.tech_matrix / iterations << std::endl;
				std::cout << "b: \n" << system.contrib_vector / iterations << std::endl;
				std::cout << "a: \n" << system.alpha << std::endl;
			}
		}
	}
	std::cout << "Avg  error : " << avg_error / L << std::endl;
	std::cout << "Avg |error|: " << avg_abs_error / L << std::endl;
}

int main(int argc, char ** argv)
{
	using Float = double;
	using RGBColor = MISExample::Spectrum<Float, 3>;

	//testEstimators<Float, Float>();
	//testEstimators<RGBColor, Float>();

	//testVirtualEstimators<Float, Float>();
	//testVirtualEstimators<RGBColor, Float>();

	//testImageEstimators<Float, Float>();
	//testImageEstimators<RGBColor, Float>();

	//testVirtualImageEstimators<Float, Float>();
	//testVirtualImageEstimators<RGBColor, Float>();


	std::cout << "Balance estimator: " << std::endl;
	testConservativeHypothesis<Float, MIS::BalanceEstimator<Float, Float>>(false);
	std::cout << "Direct  estimator: " << std::endl;
	testConservativeHypothesis<Float, MIS::DirectEstimator<Float, Float>>(true, false);
	//std::cout << "PDF instead of effective PDFs (Not expected to work)" << std::endl;
	//testConservativeHypothesis<Float, MIS::DirectEstimator<Float, Float>>(true, true);
	std::cout << "Non-Conservative Hypothesis: " << std::endl;
	testConservativeHypothesis<Float, MIS::DirectEstimator<Float, Float>>(false, false);
}