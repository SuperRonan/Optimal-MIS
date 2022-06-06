#pragma once

#include "Estimator.h"
#include "DirectCommons.h"
#include "utils/SpectrumWrapper.h"
#include <cassert>

namespace MIS
{
	template <class Spectrum, class Float=double>
	class AlphaEstimator : public Estimator<Spectrum, Float>
	{
	protected:

		using Wrapper = SpectrumWrapper<Spectrum>;
		using CWrapper = SpectrumWrapper<const Spectrum>;

	public:

		using VectorsT = Eigen::Matrix<Float, Eigen::Dynamic, Wrapper::size()>;

	protected:

		using VectorT = Eigen::Matrix<Float, Eigen::Dynamic, 1>;
		using VectorL = Eigen::Matrix<Float, 1, Eigen::Dynamic>;
		using VectorSpectrum = Eigen::Matrix<Float, 1, Wrapper::size()>;

		VectorsT _alpha;
		VectorL _samples_per_tech;

		Spectrum _estimate;

	public:

		AlphaEstimator(int N) : 
			Estimator<Spectrum, Float>(N, Heuristic::Alpha),
			_alpha(N, Wrapper::size()),
			_samples_per_tech(N),
			_estimate(0)
		{
			this->_samples_per_tech.fill(1);
			this->_alpha.fill(0);
		}

		AlphaEstimator(AlphaEstimator const&) = default;
		AlphaEstimator(AlphaEstimator&&) = default;

		void setAlpha(VectorsT const& alpha)
		{
			assert(alpha.rows() == this->m_numtechs);
			assert(alpha.cols() == Wrapper::size());
			this->_alpha = alpha;
		}

		void setAlpha(VectorsT&& alpha)
		{
			assert(alpha.rows() == this->m_numtechs);
			assert(alpha.cols() == Wrapper::size());
			this->_alpha = std::move(alpha);
		}

		virtual void setSamplesForTechnique(int tech_index, int n)override
		{
			this->_samples_per_tech[tech_index] = n;
		}

		virtual void addEstimate(Spectrum const& estimate, const Float* wb, int tech_index) override
		{
			Spectrum tmp = estimate * wb[tech_index];
			Wrapper _tmp(tmp);

			for (int k = 0; k < Wrapper::size(); ++k)
			{
				for (int i = 0; i < this->m_numtechs; ++i)
				{
					_tmp[k] += wb[i] * this->_alpha(i, k);
				}
			}

			this->_estimate += _tmp;
		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index) override
		{
			Spectrum tmp = estimate;
			Wrapper _tmp(tmp);

			for (int k = 0; k < Wrapper::size(); ++k)
			{
				_tmp[k] += this->_alpha(tech_index, k);
			}

			this->_estimate += _tmp;
		}

		Spectrum alphaSum()const
		{
			Spectrum res;
			Wrapper _res(res);
			VectorSpectrum doot = _samples_per_tech * _alpha;
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				_res[k] = doot[k];
			}
			return res;
		}

		virtual Spectrum solve(int iterations)override
		{
			return this->_estimate / float(iterations) + this->alphaSum();
		}

		virtual void reset()override
		{
			this->_estimate = 0;
		}

	};
}