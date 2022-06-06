#pragma once

#include <cassert>
#include <mutex>
#include "utils/SpectrumWrapper.h"
#include "utils/Parallel.h"
#include "Heuristics.h"

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

		int m_numtechs;

		int m_width, m_height;

		constexpr static int spectrumDim()
		{
			return Wrapper::size();
		}

	public:

		const Heuristic m_heuristic;

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
			int i = int(u * m_width);
			int j = int(v * m_height);
			return pixelTo1D(i, j);

		}
		
		template <class Function>
		MIS_forceinline void loopThroughImage(const Function& function)const
		{
			if constexpr (ROW_MAJOR)
			{
				Parallel::parallelFor(0, m_width, [&](int i) {
					for (int j = 0; j < m_height; ++j)
						function(i, j);
					});
			}
			else
			{
				Parallel::parallelFor(0, m_height, [&](int j) {
					for (int i = 0; i < m_width; ++i)
						function(i, j);
					});
			}
		}

		int numTechs()const { return m_numtechs; }

		int width()const { return m_width; }

		int height()const { return m_height; }



		using Spectrum_Type = Spectrum;
		using Float_Type = Float;

		/// <param name="N">: The number of techniques.</param>
		ImageEstimator(int N, int width, int height, Heuristic h) :
			m_numtechs(N),
			m_width(width),
			m_height(height),
			m_heuristic(h)
		{}

		ImageEstimator(ImageEstimator const& other) = default;

		ImageEstimator(ImageEstimator&& other) = default;

		/// <summary>
		/// Sets the number of sample per iteration for a technique. 
		/// </summary>
		/// <param name="tech_index">: The index of the technique.</param>
		/// <param name="n">: The number of sample per iteration.</param>
		virtual void setSamplesForTechnique(int tech_index, int n)
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

	};

} // namespace MIS