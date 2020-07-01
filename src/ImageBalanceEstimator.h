#pragma once
#include <ImageEstimator.h>

namespace MIS
{
    template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageBalanceEstimator : public ImageEstimator<Spectrum, Float, ROW_MAJOR>
	{
	protected:

		Image::Image<Spectrum, MAJOR> m_image;

		std::mutex m_mutex;

	public:

		ImageBalanceEstimator(int N, int width, int height) :
			ImageEstimator(N, width, height),
			m_image(width, height)
		{
			m_image.fill(0);
		}

		ImageBalanceEstimator(ImageBalanceEstimator const& other) :
			ImageEstimator(other)
		{
			m_image = other.m_image;
		}

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			addOneTechniqueEstimate(estimate * balance_weights[tech_index], tech_index, u, v, thread_safe_update);
		}

		virtual void loop() override
		{}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			Math::Vector<int, 2> pixel = { u * m_width, v * m_height };
			if (thread_safe_update)
				m_mutex.lock();
			m_image(pixel) += estimate;
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void solve(Image::Image<Spectrum, ROW_MAJOR>& res, int iterations) override
		{
			Parallel::ParallelFor(0, m_width * m_height,
				[&](int pixel)
				{
					res[pixel] += m_image[pixel] / iterations;
				});
		}
	};
}