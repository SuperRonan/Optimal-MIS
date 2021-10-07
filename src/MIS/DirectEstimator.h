#pragma once
#include <vector>
#include "Estimator.h"
#include "DirectCommons.h"
#include "utils/SpectrumWrapper.h"
#include "utils/settings.h"

namespace MIS
{
	/// <summary>
	/// Integer can be signed or unsigned, it does not matter.
	/// It is here to precise the width of the integers to use (we recomand 64 bits)
	/// <returns></returns>
	template <class Spectrum, class Float = double, class Integer = int64_t>
	class DirectEstimator_V1 : public Estimator<Spectrum, Float>
	{
	protected:
		
		using SINT = typename std::make_signed<Integer>::type;
		using UINT = typename std::make_unsigned<Integer>::type;

		using StorageFloat = Float;
		using SolvingFloat = Float;
		using StorageUInt = UINT;

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

		std::vector<UINT> m_sample_per_technique;

        // Pre allocation for the solving step
		mutable MatrixT m_matrix;
		mutable VectorT m_vector;
		mutable SolverT m_solver;
		mutable VectorT m_MVector;

	public:

		DirectEstimator_V1(int N) :
			Estimator<Spectrum, Float>::Estimator(N, Heuristic::Direct),
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
			reset();
		}

		DirectEstimator_V1(DirectEstimator_V1 const& other):
			Estimator<Spectrum, Float>::Estimator(other),
			msize(other.msize),
			m_data(other.m_data),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * this->m_numtechs) * sizeof(StorageFloat)))),
			m_sample_per_technique(other.m_sample_per_technique),
			m_matrix(other.m_matrix),
			m_vector(other.m_vector),
			m_solver(other.m_solver),
			m_MVector(other.m_MVector)
		{}

		DirectEstimator_V1(DirectEstimator_V1&& other) :
			Estimator<Spectrum, Float>::Estimator(std::move(other)),
			msize(other.msize),
			m_data(std::move(other.m_data)),
			m_matrix_data((StorageFloat*)m_data.data()),
			m_vectors_data(((StorageFloat*)m_data.data()) + msize),
			m_sample_count((StorageUInt*)(m_data.data() + ((msize + Wrapper::size() * this->m_numtechs) * sizeof(StorageFloat)))),
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
			for (int i = 0; i < this->m_numtechs; ++i)
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
					StorageFloat* vector = m_vectors_data + this->m_numtechs * k;
					for (int i = 0; i < this->m_numtechs; ++i)
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
				(m_vectors_data + this->m_numtechs * k)[tech_index] += _estimate[k];
			}
		}

		inline void fillMatrix(int iterations) const
		{
			for (int i = 0; i < this->m_numtechs; ++i)
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
				SINT expected = m_sample_per_technique[i] * (SINT)iterations;
				SINT actually = m_sample_count[i];
				// Use signed integer in case the actual number of samples given overflows the expected one
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
				bool is_zero = true;
				const StorageFloat* cvector = m_vectors_data + k * this->m_numtechs;
				for (int i = 0; i < this->m_numtechs; ++i)
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

		LinearSystem<Float> getLinearSystem(int iterations) 
		{
			LinearSystem<Float> res(this->m_numtechs, Wrapper::size());
			
			fillMatrix(iterations);
			m_solver = m_matrix.colPivHouseholderQr();
			res.tech_matrix = m_matrix;
			for (int k = 0; k < CWrapper::size(); ++k)
			{
				const StorageFloat* cvector = m_vectors_data + k * this->m_numtechs;
				for (int i = 0; i < this->m_numtechs; ++i)
				{
					Float elem = cvector[i];
					m_vector[i] = elem;
					res.contrib_vectors(i, k) = elem;
				}
				m_vector = m_solver.solve(m_vector);
				for (int i = 0; i < this->m_numtechs; ++i)
				{
					res.alphas(i, k) = m_vector(i);
				}
			}

			return res;
		}
	};

	/// <summary>
	/// Integer can be signed or unsigned, it does not matter.
	/// It is here to precise the width of the integers to use (we recomand 64 bits)
	/// <returns></returns>
	template <class Spectrum, class Float = double, class Integer = int64_t>
	class DirectEstimator_V2 : public Estimator<Spectrum, Float>
	{
	protected:

		using SINT = typename std::make_signed<Integer>::type;
		using UINT = typename std::make_unsigned<Integer>::type;

		using StorageFloat = Float;
		using SolvingFloat = Float;
		using StorageUInt = UINT;

		using Wrapper = SpectrumWrapper<Spectrum>;
		using CWrapper = SpectrumWrapper<const Spectrum>;

		using MatrixT = Eigen::Matrix<SolvingFloat, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorT = Eigen::Matrix<SolvingFloat, Eigen::Dynamic, 1>;
		using VectorL = Eigen::Matrix<SolvingFloat, 1, Eigen::Dynamic>;
		using VectorsT = Eigen::Matrix<SolvingFloat, Eigen::Dynamic, Wrapper::size()>;
		using VectorSpectrum = Eigen::Matrix<SolvingFloat, 1, Wrapper::size()>;
		using VectorI = Eigen::Matrix<SINT, Eigen::Dynamic, 1>;
		using SolverT = Eigen::ColPivHouseholderQR<MatrixT>;


		// msize => matrix size
		// also the vector offset
		const int msize;

		int matTo1D(int row, int col) const {
			// return row * numTechs + col;
			if (col > row) std::swap(row, col);
			return (row * (row + 1)) / 2 + col;
		}

		VectorI m_sample_per_technique;
		VectorI m_sample_count;
		VectorL m_MVector;

		MatrixT m_matrix;
		SolverT m_solver;
		// One column per spectrum channel
		VectorsT m_vectors;
		VectorsT m_alphas;

	public:

		DirectEstimator_V2(int N) :
			Estimator<Spectrum, Float>::Estimator(N, Heuristic::Direct),
			msize(N* (N + 1) / 2),
			m_sample_per_technique(N, 1),
			m_sample_count(N),
			m_MVector(N),
			m_matrix(N, N),
			m_solver(N, N),
			m_vectors(N, Wrapper::size()),
			m_alphas(N, Wrapper::size())
		{
			m_sample_per_technique.fill(1);
			m_MVector.fill(1);
			reset();
		}

		DirectEstimator_V2(DirectEstimator_V2 const& other) = default;

		DirectEstimator_V2(DirectEstimator_V2&& other) = default;

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
			// Only update the upper left triangle of the matrix
			for (int i = 0; i < this->m_numtechs; ++i)
			{
				for (int j = 0; j <= i; ++j) // Exploit the symmetry of the matrix
				{
					Float tmp = balance_weights[i] * balance_weights[j];
					m_matrix(i, j) += tmp;
				}
			}
			if (!balance_estimate.isZero())
			{
				for (int i = 0; i < this->m_numtechs; ++i)
				{
					for (int k = 0; k < Wrapper::size(); ++k)
					{
						Float tmp = balance_estimate[k] * balance_weights[i];
						m_vectors(i, k) += tmp;
					}
				}
			}
		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index)override
		{
			const CWrapper _estimate = estimate;
			// Only update one row of the vectors.
			// The matrix will be updated with unsampled samples
			for (int k = 0; k < Wrapper::size(); ++k)
			{
				m_vectors(tech_index, k) += _estimate[k];
			}
		}

		inline void completeMatrix(int iterations)
		{
			for (int i = 0; i < this->m_numtechs; ++i)
			{
				// Complete the down right triangle
				for (int j = 0; j < i; ++j)
				{
					Float elem = m_matrix(i, j);
					assert(elem >= 0);
					if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
					m_matrix(j, i) = m_matrix(i, j) = elem;
				}
				// Check Diagonal + unsampled samples
				Float elem = m_matrix(i, i);
				assert(elem >= 0);
				if (std::isnan(elem) || std::isinf(elem))	elem = 0;
				SINT expected = m_sample_per_technique[i] * (SINT)iterations;
				SINT actually = m_sample_count[i];
				// Use signed integer in case the actual number of samples given overflows the expected one
				elem += (Float)(expected - actually); // Unsampled samples
				m_matrix(i, i) = elem;
			}
			m_sample_count = m_sample_per_technique;
		}

		virtual Spectrum solve(int iterations)override
		{
			Spectrum _res = 0;
			Wrapper res = _res;
			bool is_zero = m_vectors.isZero(0);
			// Hope there are no inf nor nan
			if (!is_zero)
			{
				completeMatrix(iterations);
				m_solver = m_matrix.colPivHouseholderQr();
				m_alphas = m_solver.solve(m_vectors);
				VectorSpectrum doot = m_MVector * m_alphas;
				for (int k = 0; k < Wrapper::size(); ++k)
					res[k] = doot[k];
			}
			return res;
		}

		virtual void reset()override
		{
			m_sample_count.fill(SINT(0));
			m_matrix.fill(Float(0));
			m_vectors.fill(Float(0));
		}

		LinearSystem<Float> getLinearSystem(int iterations)
		{
			LinearSystem<Float> res(this->m_numtechs, Wrapper::size());

			completeMatrix(iterations);
			m_solver = m_matrix.colPivHouseholderQr();
			res.tech_matrix = m_matrix;

			res.contrib_vectors = m_vectors;
			res.alphas = m_solver.solve(m_vectors);

			return res;
		}
	};


#if OPTIMIS_IMPL == 1
	template <class Spectrum, class Float = double, class Integer = int64_t>
	using DirectEstimator = DirectEstimator_V1<Spectrum, Float, Integer>;
#elif OPTIMIS_IMPL == 2
	template <class Spectrum, class Float = double, class Integer = int64_t>
	using DirectEstimator = DirectEstimator_V2<Spectrum, Float, Integer>;
#else
	static_assert(false, "Bad OPTIMIS_IMPL definition (Should be 1 or 2)");
#endif
}