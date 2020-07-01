#pragma once

namespace MIS
{
	/// <summary>
	/// Interface for MIS Estimators.
	/// </summary>
	/// <typeparam name="Spectrum">: The spectrum type. See Examples/Spectrum.h for the requirements of this type.</typeparam>
	/// <typeparam name="Float">: The Floating point type used for the computation. Expecting float or double. 
	/// If other custom types are used, it should be trivially copyable. 
	/// In the reset function, it is also assumed that the 0 of Float is represented by sizeof(Float) 0 bytes. 
	/// </typeparam>
	template <class Spectrum, class Float=double>
	class Estimator
	{
	protected:

		const int m_numtechs;

	public:

		/// <param name="N"> The number of techniques.</param>
		Estimator(int N):
			m_numtechs(N)
		{}

		Estimator(Estimator const& other) = default;

		/// <summary>
		/// Sets the number of sample per iteration for a technique. 
		/// </summary>
		/// <param name="tech_index">: The index of the technique.</param>
		/// <param name="n">: The number of sample per iteration.</param>
		virtual void setSampleForTechnique(int tech_index, int n) {}

		/// <summary>
		/// Updates the estimator with the information of the drawn sample.
		/// </summary>
		/// <param name="estimate">: f(x) / p(x) , with f(x) the raw contribtion of the sample and p(x) the pdf of x being sampled by technique "tech_index".</param>
		/// <param name="balance_weights">: An array of the balance weights for each techniques.</param>
		/// <param name="tech_index">: The index of the technique that produced the sample being added.</param>
		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) = 0;

		/// <summary>
		/// Specialization of the addEstimate function in the case only one technique can produce the sample. 
		/// Then the balance weight of tech_index is 1 and 0 for all other techniques. 
		/// The update of the estimator is then faster for the optimal estimators.
		/// </summary>
		/// <param name="estimate">: f(x) / p(x) , with f(x) the raw contribtion of the sample and p(x) the pdf of x being sampled by technique "tech_index".</param>
		/// <param name="tech_index">: The index of the technique that produced the sample being added.</param>
		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index) = 0;

		/// <summary>
		/// Used by the progressive estimator only. Should be called once after each iteration. 
		/// </summary>
		virtual void loop() {}

		/// <summary>
		/// The estimator returns its estimation with this function. 
		/// </summary>
		/// <param name="iterations">: The total number of iterations the algorithm ran.</param>
		/// <returns>The estimator's estimation.</returns>
		virtual Spectrum solve(int iterations) = 0;

		/// <summary>
		/// Resets the estimator back to a "zero" state. Should be called before starting to estimate an integral. 
		/// </summary>
		virtual void reset() = 0;
	};
}