#pragma once
#include <vector>
#include <Estimator.h>
#include <Eigen/Dense>
#include <utils/SpectrumWrapper.h>
namespace MIS
{
	template <class Spectrum, class Float = double>
	class DirectEstimator : public Estimator<Spectrum, Float>
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

	public:

		DirectEstimator(int N) :
			Estimator(N),
			msize(N* (N + 1) / 2),
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
			m_MVector.fill(1);
		}

		DirectEstimator(DirectEstimator const& other):
			Estimator(other),
			msize(other.msize),
			m_data(other.m_data),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * m_numtechs) * sizeof(StorageFloat)))),
			m_sample_per_technique(other.m_sample_per_technique),
			m_matrix(other.m_matrix),
			m_vector(other.m_vector),
			m_solver(other.m_solver),
			m_MVector(other.m_MVector)
		{}

		DirectEstimator(DirectEstimator&& other) :
			Estimator(std::move(other)),
			msize(other.msize),
			m_data(std::move(other.m_data)),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * m_numtechs) * sizeof(StorageFloat)))),
			m_sample_per_technique(std::move(other.m_sample_per_technique)),
			m_matrix(std::move(other.m_matrix)),
			m_vector(std::move(other.m_vector)),
			m_solver(std::move(other.m_solver)),
			m_MVector(std::move(other.m_MVector))
		{}


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
				for (int k = 0; k < Wrapper::size(); ++k)
				{
					StorageFloat* vector = m_vectors_data + m_numtechs * k;
					for (int i = 0; i < m_numtechs; ++i)
					{
						Float tmp = balance_estimate[k] * balance_weights[i];
						vector[i] += tmp;
					}
				}
			}
		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index)override
		{
			const CWrapper _estimate = estimate;
			const int mat_index = matTo1D(tech_index, tech_index);
			++m_sample_count[tech_index];
			m_matrix_data[mat_index] += 1.0; // this is not really necessary, it could be done implicitely during the solving function, but it is cleaner.
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				(m_vectors_data + m_numtechs * k)[tech_index] += _estimate[k];
			}
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
				if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
				size_t expected = m_sample_per_technique[i] * (size_t)iterations;
				size_t actually = m_sample_count[i];
				m_matrix(i, i) = elem + (Float)(expected - actually); // Unsampled samples
			}
		}

		virtual Spectrum solve(int iterations)override
		{
			Spectrum _res = 0;
			Wrapper res = _res;
			bool matrix_solved = false;
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				bool is_zero;
				const StorageFloat* cvector = m_vectors_data + k * m_numtechs;
				for (int i = 0; i < m_numtechs; ++i)
				{
					Float elem = cvector[i];
					if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
					m_vector[i] = elem;
					is_zero = is_zero & (elem == 0);
				}
				if (!is_zero)
				{
					if (!matrix_solved)
					{
						fillMatrix(iterations);
						m_solver = m_matrix.colPivHouseholderQr();
						matrix_solved = true;
					}
					VectorT& alpha = m_vector;
					alpha = m_solver.solve(m_vector);
					SolvingFloat estimate = alpha.dot(m_MVector);
					res[k] = estimate;
				}
			}
			return res;
		}

		virtual void reset()override
		{
			std::fill(m_data.begin(), m_data.end(), 0);
		}
    };
}