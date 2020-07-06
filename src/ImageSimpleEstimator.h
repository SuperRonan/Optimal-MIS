#pragma once

#include <ImageEstimator.h>
#include <vector>

namespace MIS
{
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageSimpleEstimator : public ImageEstimator<Spectrum, Float, ROW_MAJOR>
	{
	protected:

		std::vector<Spectrum> m_buffer;

		std::mutex m_mutex;

		Spectrum& pixel(Float u, Float v)
		{
			int i = u * m_width;
			int j = v * m_height;
			int id = pixelTo1D(i, j);
			Spectrum& p = m_buffer[id];
			return p;
		}

		virtual Float weight(const Float* balance_weights, int tech_index)const = 0;

	public:

		ImageSimpleEstimator(int N, int width, int height):
			ImageEstimator(N, width, height),
			m_buffer(width * height, Spectrum(0))
		{}

		ImageSimpleEstimator(ImageSimpleEstimator const& other) = default;

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			Spectrum& p = pixel(u, v);
			if (thread_safe_update)
				m_mutex.lock();
			p += estimate * weight(balance_weights, tech_index);
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index, Float u, Float v, bool thread_safe_update = false)override
		{
			Spectrum& p = pixel(u, v);
			if (thread_safe_update)
				m_mutex.lock();
			p += estimate;
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void loop() override
		{}

		virtual void solve(Spectrum* res, int iterations)override
		{
			loopThroughImage([&](int i, int j)
			{
				int id = pixelTo1D(i, j);
				res[id] += m_buffer[id] / iterations;
			});
		}

	};
}