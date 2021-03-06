#pragma once

#include <vector>
#include "Estimator.h"
#include <Eigen/Dense>
#include "utils/SpectrumWrapper.h"
#include <algorithm>

namespace MIS
{
	template <class Spectrum, class Float = double>
	class ProgressiveEstimator : public Estimator<Spectrum, Float>
	{
	protected:

		using StorageFloat = Float;
		using SolvingFloat = Float;
		using StorageUInt = unsigned int;

		using MatrixT = Eigen::Matrix<SolvingFloat, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorT = Eigen::Matrix<SolvingFloat, Eigen::Dynamic, 1>;
		using SolverT = Eigen::ColPivHouseholderQR<MatrixT>;

		using Wrapper = SpectrumWrapper<Spectrum>;
		using CWrapper = SpectrumWrapper<const Spectrum>;

		// msize => matrix size
		// also the vector offset
		const int msize;

		int matTo1D(int row, int col) const {
			// return row * numTechs + col;
			if (col > row) std::swap(row, col);
			return (row * (row + 1)) / 2 + col;
		}

		// Re-solve the system every U iteration
		// U: Update step
		const int U;
		int m_loop_counter;

		Spectrum m_result;

		// One contiguous array to store all the data of the estimator
		// So usually all the data can be stored in one line of cache
		std::vector<char> m_data;

		StorageFloat* m_matrix_data;
		StorageFloat* m_vectors_data;
		StorageUInt* m_sample_count;

		std::vector<unsigned int> m_sample_per_technique;

		// Pre allocation for the solving step
		MatrixT m_matrix;
		VectorT m_vector;
		SolverT m_solver;
		VectorT m_MVector;
		VectorT m_alphas[Wrapper::size()];

	public:

		ProgressiveEstimator(int N, int U=4) :
			Estimator(N, Heuristic::Progressive),
			msize(N* (N + 1) / 2),
			U(U),
			m_loop_counter(0),
			m_result(0),
			m_data((msize + Wrapper::size() * N) * sizeof(StorageFloat) + N * sizeof(StorageUInt), 0),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * N) * sizeof(StorageFloat)))),
			m_sample_per_technique(N, 1),
			m_matrix(N, N),
			m_vector(N),
			m_solver(N, N),
			m_MVector(N)
		{
			const VectorT zero = VectorT(N).setZero();
			std::fill_n(m_alphas, Wrapper::size(), zero);
			m_MVector.fill(1);
		}

		ProgressiveEstimator(ProgressiveEstimator const& other) :
			Estimator(other),
			msize(other.msize),
			U(other.U),
			m_loop_counter(other.m_loop_counter),
			m_result(other.m_result),
			m_data(other.m_data),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * m_numtechs) * sizeof(StorageFloat)))),
			m_sample_per_technique(other.m_sample_per_technique),
			m_matrix(other.m_matrix),
			m_vector(other.m_vector),
			m_solver(other.m_solver),
			m_MVector(other.m_MVector)
		{
			std::copy_n(other.m_alphas, Wrapper::size(), m_alphas);
		}

		ProgressiveEstimator(ProgressiveEstimator&& other) :
			Estimator(std::move(other)),
			msize(other.msize),
			U(other.U),
			m_loop_counter(other.m_loop_counter),
			m_result(std::move(other.m_result)),
			m_data(std::move(other.m_data)),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * m_numtechs) * sizeof(StorageFloat)))),
			m_sample_per_technique(std::move(other.m_sample_per_technique)),
			m_matrix(std::move(other.m_matrix)),
			m_vector(std::move(other.m_vector)),
			m_solver(std::move(other.m_solver)),
			m_MVector(std::move(other.m_MVector))
		{
			for (int i = 0; i < Wrapper::size(); ++i)	m_alpha[i] = std::move(other.m_alpha[i]);
		}

		virtual void setSampleForTechnique(int tech_index, int n)override
		{
			m_sample_per_technique[tech_index] = n;
			m_MVector[tech_index] = n;
		}

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index) override
		{
			assert(balance_weights != nullptr);
			const Spectrum _balance_estimate = estimate * balance_weights[tech_index];
			const CWrapper balance_estimate = _balance_estimate;
			++m_sample_count[tech_index];
			for (int i = 0; i < m_numtechs; ++i)
			{
				for (int j = 0; j <= i; ++j) // Exploit the symmetry of the matrix
				{
					const int mat_index = matTo1D(i, j);
					Float tmp = balance_weights[i] * balance_weights[j];
					m_matrix_data[mat_index] += tmp;
				}
			}
			if (!balance_estimate.isZero())
			{
				// Update the contrib vector
				for (int k = 0; k < Wrapper::size(); ++k)
				{
					StorageFloat* vector = m_vectors_data + m_numtechs * k;
					for (int i = 0; i < m_numtechs; ++i)
					{
						Float tmp = balance_estimate[k] * balance_weights[i];
						vector[i] += tmp;
					}
				}

				// Evaluate partial Fo (equation 10) for this sample
				// Only the right hand side of the equation (the part in the sum)
				Spectrum Fo = 0;
				Wrapper _Fo(Fo);
				for (int k = 0; k < Wrapper::size(); ++k)
				{
					Float& fo = _Fo[k];
					VectorT& alpha = m_alphas[k];
					// Dot(alpha, Wb)
					Float doot = 0;
					for (int i = 0; i < m_numtechs; ++i)	doot += balance_weights[i] * alpha[i];
					fo = balance_estimate[k] - doot;
				}
				m_result += Fo;
			}

		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index)override
		{
			const CWrapper _estimate = estimate;
			const int mat_index = matTo1D(tech_index, tech_index);
			++m_sample_count[tech_index];
			m_matrix_data[mat_index] += 1.0; // this is not really necessary, it could be done implicitely during the solving function, but it is cleaner.
			Spectrum Fo = 0;
			Wrapper _Fo(Fo);
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				(m_vectors_data + m_numtechs * k)[tech_index] += _estimate[k];
				_Fo[k] = (_estimate[k] - (m_alphas[k][tech_index]));
			}
			m_result += Fo;
		}

		inline void fillMatrix(int iterations)
		{
			for (int i = 0; i < m_numtechs; ++i)
			{
				for (int j = 0; j < i; ++j)
				{
					const int mat_id = matTo1D(i, j);
					Float elem = m_matrix_data[mat_id];
					assert(elem >= 0);
					if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
					m_matrix(i, j) = elem;
					m_matrix(j, i) = elem;
				}
				const int mat_id = matTo1D(i, i);
				Float elem = m_matrix_data[mat_id];
				assert(elem >= 0);
				if (std::isnan(elem) || std::isinf(elem))	elem = 0;
				size_t expected = m_sample_per_technique[i] * (size_t)iterations;
				size_t actually = m_sample_count[i];
				m_matrix(i, i) = elem + (Float)(expected - actually); // Unsampled samples
			}
		}

		virtual Spectrum solve(int iterations)override
		{
			return m_result / Float(iterations);
		}

		virtual void loop() override
		{
			Wrapper result = m_result;
			// Partial estimate of Fo: right hand side: sum of alpha N
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				result[k] += m_alphas[k].dot(m_MVector);
			}
			
			// Update the estimate of alpha 
			if (m_loop_counter != 0 && m_loop_counter % U == 0)
			{
				bool matrix_solved = false;
				for (int k = 0; k < Wrapper::size(); ++k)
				{
					VectorT& alpha = m_alphas[k];
					bool is_zero;
					const StorageFloat* cvector = m_vectors_data + k * m_numtechs;
					for (int i = 0; i < m_numtechs; ++i)
					{
						Float elem = cvector[i];
						if (std::isnan(elem) || std::isinf(elem))	elem = 0;
						m_vector[i] = elem;
						is_zero = is_zero & (elem == 0);
					}
					if (!is_zero)
					{
						if (!matrix_solved)
						{
							fillMatrix(m_loop_counter+1);
							m_solver = m_matrix.colPivHouseholderQr();
							matrix_solved = true;
						}
						alpha = m_solver.solve(m_vector);
					}
					else
						alpha.setZero();
				}
			}
			++m_loop_counter;
		}

		virtual void reset()override
		{
			std::fill(m_data.begin(), m_data.end(), 0);
			for (int k = 0; k < Wrapper::size(); ++k)	m_alphas[k].setZero();
			m_result = 0;
		}
	};
}