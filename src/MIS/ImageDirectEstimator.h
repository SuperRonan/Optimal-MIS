#pragma once

#include "ImageEstimator.h"
#include "DirectCommons.h"
#include "utils/settings.h"
#include "utils/SpectrumWrapper.h"

namespace MIS
{
    template <class Spectrum, class Float, class Uint, bool ROW_MAJOR>
	class ImageDirectEstimator : public ImageEstimator<Spectrum, Float, ROW_MAJOR>
	{
	protected:

		using solvingFloat = Float;
		using MatrixT = Eigen::Matrix<solvingFloat, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorT = Eigen::Matrix<solvingFloat, Eigen::Dynamic, 1>;
		
		using StorageUInt = Uint;
		using StorageFloat = Float;

		using Wrapper = SpectrumWrapper<Spectrum>;
		using CWrapper = SpectrumWrapper<const Spectrum>;

		const int msize;

		MIS_forceinline int matTo1D(int row, int col) const {
			// return row * numTechs + col;
			if (col > row) std::swap(row, col);
			return (row * (row + 1)) / 2 + col;
		}

		struct PixelData
		{
			StorageFloat* techMatrix;
			// The three channels are one after the other [RRRRRRRR GGGGGGGG BBBBBBBBB]
			StorageFloat* contribVector;
			StorageUInt* sampleCount;

			PixelData(StorageFloat* tm, StorageFloat* cv, StorageUInt* sc) :
				techMatrix(tm),
				contribVector(cv),
				sampleCount(sc)
			{}

			PixelData(const StorageFloat* tm, const StorageFloat* cv, const StorageUInt* sc) :
				techMatrix((StorageFloat*)tm),
				contribVector((StorageFloat*)cv),
				sampleCount((StorageUInt*)sc)
			{}
		};

		PixelData getPixelData(int index)const
		{
#if OPTIMIS_ONE_CONTIGUOUS_ARRAY
			const size_t offset = index * m_pixel_data_size;
			const char* address = (m_data.data() + offset);
			PixelData res((StorageFloat*)address, (StorageFloat*)(address + m_vector_ofsset), (StorageUInt*)+(address + m_counter_offset));
			return res;
#else
			const StorageFloat* matrix = m_matrices.data() + msize * index;
			const StorageFloat* vector = m_vectors.data() + m_numtechs * spectrumDim() * index;
			const StorageUInt* counts = m_sampleCounts.data() + m_numtechs * index;
			PixelData res(matrix, vector, counts);
			return res;
#endif
		}


		PixelData getPixelData(Float u, Float v)const
		{
			const int index = this->pixelTo1D(u, v);
			return getPixelData(index);
		}

		PixelData getPixelData(int i, int j)const
		{
			const int index = this->pixelTo1D(i, j);
			return getPixelData(index);
		}

		

#if OPTIMIS_ONE_CONTIGUOUS_ARRAY
		// In Byte
		const Uint m_pixel_data_size;
		const Uint m_vector_ofsset;
		const Uint m_counter_offset;

		std::vector<char> m_data;
#else
		std::vector<StorageFloat> m_matrices;
		std::vector<StorageFloat> m_vectors;
		std::vector<StorageUInt> m_sampleCounts;
#endif

		//describes how many times more samples each technique has
		std::vector<Uint> m_sample_per_technique;

		std::mutex m_mutex;

	public:

#if OPTIMIS_ONE_CONTIGUOUS_ARRAY
		ImageDirectEstimator(int N, int width, int height) :
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(N, width, height, Heuristic::Direct),
			msize(N* (N + 1) / 2),
			m_pixel_data_size(msize * sizeof(StorageFloat) + this->spectrumDim() * this->m_numtechs * sizeof(StorageFloat) + this->m_numtechs * sizeof(StorageUInt)),
			m_vector_ofsset(msize * sizeof(StorageFloat)),
			m_counter_offset(msize * sizeof(StorageFloat) + this->spectrumDim() * this->m_numtechs * sizeof(StorageFloat)),
			m_data(std::vector<char>(width* height* m_pixel_data_size, (char)0)),
			m_sample_per_technique(std::vector<Uint>(this->m_numtechs, 1))
		{}

		ImageDirectEstimator(ImageDirectEstimator const& other) :
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(other),
			msize(other.msize),
			m_pixel_data_size(other.m_pixel_data_size),
			m_vector_ofsset(other.m_vector_ofsset),
			m_counter_offset(other.m_counter_offset),
			m_data(other.m_data),
			m_sample_per_technique(other.m_sample_per_technique)
		{}

		ImageDirectEstimator(ImageDirectEstimator && other):
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(other),
			msize(other.msize),
			m_pixel_data_size(other.m_pixel_data_size),
			m_vector_ofsset(other.m_vector_ofsset),
			m_counter_offset(other.m_counter_offset),
			m_data(std::move(other.m_data)),
			m_sample_per_technique(std::move(other.m_sample_per_technique))
		{}
#else
		ImageDirectEstimator(int N, int width, int height) :
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(N, width, height, Heuristic::Direct),
			msize(N* (N + 1) / 2),
			m_sample_per_technique(std::vector<Uint>(this->m_numtechs, 1))
		{
			int res = width * height;
			m_matrices = std::vector<StorageFloat>(res * msize, (StorageFloat)0.0);
			m_vectors = std::vector<StorageFloat>(res * this->m_numtechs * this->spectrumDim(), (StorageFloat)0.0);
			m_sampleCounts = std::vector<StorageUInt>(res * this->m_numtechs, (StorageUInt)0);
		}

		ImageDirectEstimator(ImageDirectEstimator const& other) :
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(other),
			msize(other.msize),
			m_matrices(other.m_matrices),
			m_vectors(other.m_vectors),
			m_sampleCounts(other.m_sampleCounts),
			m_sample_per_technique(other.m_sample_per_technique)
		{}

		ImageDirectEstimator(ImageDirectEstimator && other) :
			ImageEstimator<Spectrum, Float, ROW_MAJOR>::ImageEstimator(other),
			msize(other.msize),
			m_matrices(std::move(other.m_matrices)),
			m_vectors(std::move(other.m_vectors)),
			m_sampleCounts(std::move(other.m_sampleCounts)),
			m_sample_per_technique(std::move(other.m_sample_per_technique))
		{}
#endif

		virtual void setSampleForTechnique(int techIndex, int n) override
		{
			m_sample_per_technique[techIndex] = n;
		}

		void checkSample(Spectrum const& balance_estimate, Float* balance_weights, int tech_index)
		{
			for (int i = 0; i < this->m_numtechs; ++i)
			{
				Float weight = balance_weights[i];
				if (weight < 0 || std::isnan(weight) || std::isinf(weight))
				{
					std::cout << weight << std::endl;
				}
			}
		}

		virtual void addEstimate(Spectrum const& estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			const Spectrum _balance_estimate = estimate * balance_weights[tech_index];
			const CWrapper balance_estimate = _balance_estimate;
			PixelData data = this->getPixelData(u, v);
			if (thread_safe_update)
				m_mutex.lock();
			++data.sampleCount[tech_index];
			for (int i = 0; i < this->m_numtechs; ++i)
			{
				for (int j = 0; j <= i; ++j)
				{
					const int mat_index = matTo1D(i, j);
					Float tmp = balance_weights[i] * balance_weights[j];
					data.techMatrix[mat_index] += tmp;
				}
			}
			if (!balance_estimate.isZero())
			{
				for (int k = 0; k < this->spectrumDim(); ++k)
				{
					StorageFloat* vector = data.contribVector + k * this->m_numtechs;
					for (int i = 0; i < this->m_numtechs; ++i)
					{
						Float tmp = balance_weights[i] * balance_estimate[k];
						vector[i] += tmp;
					}
				}
			}
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void addOneTechniqueEstimate(Spectrum const& estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			const CWrapper _estimate = estimate;
			PixelData data = this->getPixelData(u, v);
			int mat_index = matTo1D(tech_index, tech_index);
			if (thread_safe_update)
				m_mutex.lock();
			++data.sampleCount[tech_index];
			data.techMatrix[mat_index] += 1.0;
			for (int k = 0; k < this->spectrumDim(); ++k)
			{
				StorageFloat* vector = data.contribVector + k * this->m_numtechs;
				vector[tech_index] += _estimate[k];
			}
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void loop() override
		{}

		inline void fillMatrix(MatrixT& matrix, PixelData const& data, int iterations)const
		{
			for (int i = 0; i < this->m_numtechs; ++i)
			{
				for (int j = 0; j < i; ++j)
				{
					Float elem = data.techMatrix[matTo1D(i, j)];
					if (std::isnan(elem) || std::isinf(elem))	elem = 0;
					matrix(i, j) = elem;
					matrix(j, i) = elem;
				}
				matrix(i, i) = data.techMatrix[matTo1D(i, i)];
				// Unsampled samples
				Uint expected = m_sample_per_technique[i] * (Uint)iterations;
				Uint actually = data.sampleCount[i];
				matrix(i, i) += (Float)(expected - actually);
			}
		}


		virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const override
		{
#if OPTIMIS_ENABLE_DEBUG
			if (col_sum)
				saveColSum(iterations);
			if (matrix)
				saveMatrices(iterations);
			if (vec)
				saveVectors(iterations);
			if (alpha)
				saveAlphas(iterations);
#else 
			std::cout << "Warning! Image Direct Estimator Debug disabled. Nothing happens!" << std::endl;
#endif
		}

		virtual void solve(Spectrum * res, int iterations) override
		{
			using Solver = Eigen::ColPivHouseholderQR<MatrixT>;
			VectorT MVector(this->m_numtechs);
			for (int i = 0; i < this->m_numtechs; ++i)	MVector(i) = m_sample_per_technique[i];
			std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(this->m_numtechs, this->m_numtechs));
			std::vector<VectorT> vectors = Parallel::preAllocate(VectorT(this->m_numtechs));
			std::vector<Solver> solvers = Parallel::preAllocate(Solver(this->m_numtechs, this->m_numtechs));
			int resolution = this->m_width * this->m_height;
			this->loopThroughImage([&](int i, int j)
				{
					int tid = Parallel::tid();
					MatrixT& matrix = matrices[tid];
					VectorT& vector = vectors[tid];
					Solver& solver = solvers[tid];

					Spectrum _color = 0;
					Wrapper color = _color;

					int pixel_id = this->pixelTo1D(i, j);
					PixelData data = this->getPixelData(pixel_id);

					bool matrix_done = false;

					for (int k = 0; k < this->spectrumDim(); ++k)
					{
						bool isZero = true;
						const StorageFloat* contribVector = data.contribVector + k * this->m_numtechs;
						for (int i = 0; i < this->m_numtechs; ++i)
						{
							Float elem = contribVector[i];
							if (std::isnan(elem) || std::isinf(elem))	elem = 0;
							vector[i] = elem;
							isZero = isZero & (contribVector[i] == 0);
						}
						if (!isZero)
						{
							if (!matrix_done)
							{
								fillMatrix(matrix, data, iterations);
								solver = matrix.colPivHouseholderQr();
								matrix_done = true;
							}
							// Now vector is alpha
							vector = solver.solve(vector);

							solvingFloat estimate = vector.dot(MVector);
							color[k] = estimate;
						}
					}
					res[pixel_id] += color;
				});
		}

		LinearSystem<Float> getPixelLinearSystem(int iterations, int x, int y)
		{
			using Solver = Eigen::ColPivHouseholderQR<MatrixT>;
			LinearSystem<Float> res;
			res.tech_matrix = MatrixT(this->m_numtechs, this->m_numtechs);
			res.contrib_vectors = MatrixT(this->m_numtechs, this->spectrumDim());
			res.alphas = MatrixT(this->m_numtechs, this->spectrumDim());
			const PixelData data = this->getPixelData(this->pixelTo1D(x, y));
			fillMatrix(res.tech_matrix, data, iterations);
			Solver solver = res.tech_matrix.colPivHouseholderQr();
			for (int k = 0; k < this->spectrumDim(); ++k)
			{
				const StorageFloat* contrib_vector = data.contribVector + k * this->m_numtechs;
				for (int i = 0; i < this->m_numtechs; ++i)
				{
					Float elem = contrib_vector[i];
					res.contrib_vectors(i, k) = elem;
				}
			}
			res.alphas = solver.solve(res.contrib_vectors);
			return res;
		}
	};
}