#pragma once

#include <Parallel.h>
#include <cassert>
#include <mutex>

namespace MIS
{
	// Major:
	// true -> row major
	// false -> col major
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageEstimator
	{
	protected:

		const int m_numtechs;

		const int m_width, m_height;

		constexpr static int spectrumDim()
		{
			return Spectrum::size();
		}

		int PixelTo1D(int i, int j)const
		{
			if constexpr (ROW_MAJOR) // row major
				return i * m_height + j;
			else // col major
				return j * m_width + i;
		}

	public:

		ImageEstimator(int N, int width, int height) :
			m_numtechs(N),
			m_width(width),
			m_height(height)
		{}

		ImageEstimator(ImageEstimator const& other) = default;

		virtual void setSampleForTechnique(int techIndex, int n)
		{}

		virtual void addEstimate(Spectrum const& balance_estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		virtual void addOneTechniqueEstimate(Spectrum const& balance_estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		virtual void loop() = 0;

		virtual void solve(Image::Image<Spectrum, ROW_MAJOR>& res, int iterations) = 0;

		virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const
		{}
	};

} // namespace MIS