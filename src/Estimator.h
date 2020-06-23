#pragma once

namespace MIS
{
	template <class Spectrum, class Float=double>
	class Estimator
	{
	protected:

		const int m_numtechs;

	public:

		Estimator(int N):
			m_numtechs(N)
		{}

		Estimator(Estimator const& other) = default;

		virtual void setSampleForTechnique(int techIndex, int n) {}

		virtual void addEstimate(Spectrum const& balance_estimate, const Float* balance_weights, int tech_index) = 0;

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index) = 0;

		// For the Progressive estimator
		virtual void loop() {}

		virtual Spectrum solve(int iterations) = 0;

		virtual void reset() = 0;
	};
}