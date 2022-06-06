
#include <iostream>
#include "Spectrum.h"

#include <MIS/Estimators.h>
		  
#include <MIS/ImageEstimators.h>

#include <string>
#include <random>
#include <memory>


const std::vector<MIS::Heuristic> heuristics = { MIS::Heuristic::Balance,
		MIS::Heuristic::CutOff, MIS::Heuristic::Maximum, MIS::Heuristic::Power,
		MIS::Heuristic::Naive, MIS::Heuristic::Direct, MIS::Heuristic::Progressive, MIS::Heuristic::Alpha };

const std::vector<MIS::Heuristic> image_heuristics = { MIS::Heuristic::Balance,
		MIS::Heuristic::CutOff, MIS::Heuristic::Maximum, MIS::Heuristic::Power,
		MIS::Heuristic::Naive, MIS::Heuristic::Direct};


template <class Float>
class SamplingTechnique
{
public:
	// Assume s in [0, 1]
	virtual Float sample(Float s, Float& pdf) const = 0;

	virtual Float pdf(Float x) const = 0;
};

template <class Float>
class UniformTechnique : public SamplingTechnique<Float>
{
	Float _min = 0, _max = 1;

public:

	UniformTechnique(Float min=0, Float max=1):
		_min(min),
		_max(max)
	{}

	virtual Float sample(Float s, Float& pdf) const final override
	{
		Float res = _min + s * (_max - _min); 
		pdf = this->pdf(res);
		return res;
	}

	virtual Float pdf(Float x) const final override
	{
		bool in = x <= _max && x >= _min;
		return in ? Float(1) / (_max - _min) : 0;
	}
};

template <class Float>
class LinearTechnique : public SamplingTechnique<Float>
{
	Float _x1, _y1, _x2, _y2;

public:

	//    |
	// y2 |            x
	//    |           / 
	//    |          /
	//    |         /
	//    |        /
	//    |       /
	//    |      /
	//    |     /
	// y1 |    x 
	//    |
	//----+------------------
	//    |    x1     x2

	LinearTechnique(Float x1=0, Float y1=0, Float x2=1, Float y2=1):
		_x1(x1),
		_y1(y1),
		_x2(x2),
		_y2(y2)
	{}

	virtual Float sample(Float s, Float& pdf) const final override
	{
		Float a = (_y2 - _y1) / (_x2 - _x1);
		Float b = _y1 - a * _x1;
		Float c = 0.5 * _x1 * _x1 + b * _x1;
		Float integral = 0.5 * a * (_x2 * _x2 - _x1 * _x1) + b * (_x2 - _x1);
		a /= integral;
		b /= integral;
		c /= integral;
		Float res = (std::sqrt(-2 * a * c + 2 * a * s + b * b) - b) / a;
		pdf = this->pdf(res);
		return res;
	}

	virtual Float pdf(Float x) const final override
	{
		bool in = x >= _x1 && x <= _x2;
		if (!in)	return 0;
		Float a = (_y2 - _y1) / (_x2 - _x1);
		Float b = _y1 - a * _x1;
		Float integral = 0.5 * a * (_x2 * _x2 - _x1 * _x1) + b * (_x2 - _x1);
		return (a * x + b) / integral;
	}
};

template <class Spectrum, class Float>
Spectrum f(Float x)
{
	Spectrum res;
	MIS::SpectrumWrapper<Spectrum> _res(res);
	for (int i = 0; i < _res.size(); ++i)
	{
		_res[i] = std::pow(x, i);
	}
	return res;
}

template <class Spectrum, class Float>
Spectrum computeFGroundTruth(int n=4096)
{
	std::uniform_real_distribution<Float> sampler(0, 1);
	std::mt19937_64 rng;
	Spectrum res = 0;
	for (int i = 0; i < n; ++i)
	{
		res += f<Spectrum>(sampler(rng));
	}
	return res / Float(n);
}

template <class Estimator>
void testEstimator()
{
	using Float = typename Estimator::Float_Type;
	using Spectrum = typename Estimator::Spectrum_Type;
	std::cout << "Ground truth: " << computeFGroundTruth<Spectrum, Float>() << std::endl;
	const auto printSystem = [&](MIS::LinearSystem<Float> const& system)
	{
		std::cout << "matrix: \n";
		std::cout << system.tech_matrix << "\n";
		std::cout << "vectors: \n";
		std::cout << system.contrib_vectors << "\n";
		std::cout << "alphas: \n";
		std::cout << system.alphas << "\n";
	};
	// Just a function to check that it compiles well
	const int N = 3;
	Estimator _estimator(N);
	Estimator estimator = _estimator;
	if(estimator.m_heuristic == MIS::Heuristic::Alpha)
	{
		MIS::AlphaEstimator<Spectrum, Float>* ae = (MIS::AlphaEstimator<Spectrum, Float>*)&estimator;
		Eigen::Matrix<Float, Eigen::Dynamic, Estimator::Wrapper_Type::size()> alpha(N, Estimator::Wrapper_Type::size());
		alpha.fill(0);
		ae->setAlpha(std::move(alpha));
	}

	std::vector<std::shared_ptr<SamplingTechnique<Float>>> techniques = {
		std::make_shared<UniformTechnique<Float>>(),
		std::make_shared<LinearTechnique<Float>>(0, 0, 0.5, 0.2),
		std::make_shared<LinearTechnique<Float>>(0.5, 0.2, 1, 1),
	};

	const int samples_per_tech[] = { 1, 2, 2 };
	for (int i = 0; i < N; ++i)
	{
		estimator.setSamplesForTechnique(i, samples_per_tech[i]);
	}

	std::uniform_real_distribution<Float> sampler(0, 1);
	std::mt19937_64 rng;

	const int iterations = 16*16;
	
	for (int n = 0; n < iterations; ++n)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < samples_per_tech[i]; ++j)
			{
				Float wb[N];
				Float* qdf = wb;
				Float sample = techniques[i]->sample(sampler(rng), qdf[i]);
				qdf[i] *= samples_per_tech[i];
				Float sum = qdf[i];
				Spectrum estimate = f<Spectrum>(sample) / (qdf[i]);
				for (int ii = 0; ii < N; ++ii)
				{
					if (i != ii)
					{
						qdf[ii] = techniques[ii]->pdf(sample) * samples_per_tech[ii];
						sum += qdf[ii];
					}
				}
				if (sum == 0)	__debugbreak();
				for (int ii = 0; ii < N; ++ii)
				{
					wb[ii] = qdf[ii] / sum;
				}
				estimator.addEstimate(estimate, wb, i);
			}
		}
		estimator.loop();
	}

	std::cout << "Estimator " << (int)estimator.m_heuristic << " result: " << estimator.solve(iterations) << std::endl;

	if constexpr (std::is_same<Estimator, MIS::DirectEstimator< Spectrum, Float>>::value)
	{
		MIS::LinearSystem<Float> system = estimator.getLinearSystem(1);
		//printSystem(system);
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
	testEstimator<MIS::AlphaEstimator<Spectrum, Float>>();
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
		if (estimator->m_heuristic == MIS::Heuristic::Alpha)
		{
			MIS::AlphaEstimator<Spectrum, Float>* ae = (MIS::AlphaEstimator<Spectrum, Float>*)estimator;
			Eigen::Matrix<Float, Eigen::Dynamic, MIS::Estimator<Spectrum, Float>::Wrapper_Type::size()> alpha(eci.N, MIS::Estimator<Spectrum, Float>::Wrapper_Type::size());
			alpha.fill(1);
			ae->setAlpha(std::move(alpha));
		}

		Spectrum estimate;
		Float weights[N];
		
		estimator->setSamplesForTechnique(1, 2);

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
	estimator.setSamplesForTechnique(1, 2);
	
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
		estimator->setSamplesForTechnique(1, 2);

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

	constexpr bool a = MIS::isOptimal(MIS::Heuristic::Balance);
	constexpr bool b = MIS::isOptimal(MIS::Heuristic::Direct);
}