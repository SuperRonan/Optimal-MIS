#pragma once

#include <Parallel.h>
#include <cassert>
#include <mutex>
#include <SpectrumWrapper.h>
#include <Parallel.h>

namespace MIS
{
	/// <summary>
	/// Base interface for the Image estimator type. 
	/// It is faster than an Image of estimator (in the sense of type like "Image&lt;tEstimator&gt;"), mostly for the memory management. 
	/// Storing a whole image of estimators is also necessary when using sampling techniques such as light tracing.
	/// </summary>
	/// <typeparam name="Spectrum">: The spectrum type. See Examples/Spectrum.h for the requirements of this type.</typeparam>
	/// <typeparam name="Float">: The Floating point type used for the computation. Expecting float or double. 
	/// If other custom types are used, it should be trivially copyable. 
	/// In the reset function, it is also assumed that the 0 of Float is represented by sizeof(Float) 0 bytes. 
	/// </typeparam>
	/// <typeparam name="ROW_MAJOR">: Indicates whether the memory array of the image should be in row or col major.</typeparam>
	template <class Spectrum, class Float, bool ROW_MAJOR>
	class ImageEstimator
	{
	protected:

		using Wrapper = SpectrumWrapper<Spectrum>;

		const int m_numtechs;

		const int m_width, m_height;

		constexpr static int spectrumDim()
		{
			return Wrapper::size();
		}

		/// <summary>
		/// Returns the 1D index of a pixel in the image array
		/// </summary>
		/// <param name="i">: Index of the col: usualy from left to right: -</param>
		/// <param name="j">: Index of the row: usualy from top to bottom: |</param>
		/// <returns>: The 1D index</returns>
		int pixelTo1D(int i, int j)const
		{
			if constexpr (ROW_MAJOR) // row major
				return i * m_height + j;
			else // col major
				return j * m_width + i;
		}

		int pixelTo1D(Float u, Float v)const
		{
			int i = u * m_width;
			int j = v * m_height;
			return pixelTo1D(i, j);

		}
		
		template <class Function>
		__forceinline void loopThroughImage(const Function& function)const
		{
			if constexpr (ROW_MAJOR)
			{
				Parallel::ParallelFor(0, m_width, [&](int i) {
					for (int j = 0; j < m_height; ++j)
						function(i, j);
					});
			}
			else
			{
				Parallel::ParallelFor(0, m_height, [&](int j) {
					for (int i = 0; i < m_width; ++i)
						function(i, j);
					});
			}
		}

	public:

		using Spectrum_Type = Spectrum;
		using Float_Type = Float;

		/// <param name="N">: The number of techniques.</param>
		ImageEstimator(int N, int width, int height) :
			m_numtechs(N),
			m_width(width),
			m_height(height)
		{}

		ImageEstimator(ImageEstimator const& other) = default;

		/// <summary>
		/// Sets the number of sample per iteration for a technique. 
		/// </summary>
		/// <param name="tech_index">: The index of the technique.</param>
		/// <param name="n">: The number of sample per iteration.</param>
		virtual void setSampleForTechnique(int tech_index, int n)
		{}

		/// <summary>
		/// Updates the estimator with the information of the drawn sample at a position (u, v) in the image.
		/// </summary>
		/// <param name="estimate">: f(x) / p(x) , with f(x) the raw contribtion of the sample and p(x) the pdf of x being sampled by technique "tech_index".</param>
		/// <param name="balance_weights">: An array of the balance weights for each techniques.</param>
		/// <param name="tech_index">: The index of the technique that produced the sample being added.</param>
		/// <param name="u">: x axis (left to right) coordinate normalized in [0, 1[</param>
		/// <param name="v">: y axis (top to bottom) coordinate normalized in [0, 1[</param>
		/// <param name="thread_safe_update">: In case of potential concurent access to the same pixel (from the light tracer for instance). It is not necessary.</param>
		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		/// <summary>
		/// Specialization of the addEstimate function in the case only one technique can produce the sample. 
		/// Then the balance weight of tech_index is 1 and 0 for all other techniques. 
		/// The update of the estimator is then faster for the optimal estimators.
		/// </summary>
		/// <param name="estimate">: f(x) / p(x) , with f(x) the raw contribtion of the sample and p(x) the pdf of x being sampled by technique "tech_index".</param>
		/// <param name="tech_index">: The index of the technique that produced the sample being added.</param>
		/// <param name="u">: x axis (left to right) coordinate normalized in [0, 1[</param>
		/// <param name="v">: y axis (top to bottom) coordinate normalized in [0, 1[</param>
		/// <param name="thread_safe_update">: In case of potential concurent access to the same pixel (from the light tracer for instance). It is not necessary.</param>
		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		/// <summary>
		/// Used by the progressive estimator only. Should be called once after each iteration. 
		/// </summary>
		virtual void loop() = 0;

		/// <summary>
		/// The estimator fills an array of results with its estimation for every pixel. 
		/// </summary>
		/// <param name="res">: A 1D array of size width * height to be add-filled with the result of the estimator. For each pixel, the estimator's estimate is added to the pixel in the res array.</param>
		/// <param name="iterations">: The total number of iterations the algorithm ran.</param>
		virtual void solve(Spectrum* res, int iterations) = 0;

		/// <summary>
		/// Function to launch the debug, only for the direct estimator. This function is currently deprecated.
		/// </summary>
		/// <param name="iterations">: The number of iterations the algorithm ran.</param>
		/// <param name="col_sum">: Output an image of the difference between the estimated sum of the cols / rows of the tech matrix and its expected value (number of sample for the technique). 
		/// If the image does not average to 0, there is probably a bias.</param>
		/// <param name="matrix">: Output images of the tech matices.</param>
		/// <param name="vec">: Output images of the contrib vectors.</param>
		/// <param name="alpha">: Output images of the alpha vectors.</param>
		virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const
		{}
	};

} // namespace MIS