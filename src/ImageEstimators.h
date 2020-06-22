#pragma once

#include <Image/Image.h>
#include <System/Parallel.h>
#include <cassert>
#include <mutex>
#include <Eigen/Dense>

#define PRINT(var) std::cout << #var << ": " << var << std::endl;

namespace MIS
{
	// Major:
	// true -> row major
	// false -> col major
	template <class Spectrum, class Float, bool MAJOR>
	class ImageEstimator
	{
	protected:

		const int m_numtechs;

		const int m_width, m_height;

		constexpr static int spectrumDim()
		{
			return Spectrum::nSamples;
		}

		int PixelTo1D(int i, int j)const
		{
			if constexpr (MAJOR) // row major
				return i * m_height + j;
			else // col major
				return j * m_width + i;
		}

	public:

		ImageEstimator(int N, int width, int height) :
			m_numtechs(N),
			m_width(width),
			m_height(height)
		{}

		ImageEstimator(ImageEstimator const& other) = default;

		virtual void setSampleForTechnique(int techIndex, int n)
		{}

		virtual void addEstimate(Spectrum const& balance_estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		virtual void addOneTechniqueEstimate(Spectrum const& balance_estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) = 0;

		virtual void loop() = 0;

		virtual void solve(Image::Image<Spectrum, MAJOR>& res, int iterations) = 0;

		virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const
		{}
	};

	template <class Spectrum, class Float, bool MAJOR>
	class BalanceEstimatorImage : public ImageEstimator<Spectrum, Float, MAJOR>
	{
	protected:

		Image::Image<Spectrum, MAJOR> m_image;

		std::mutex m_mutex;

	public:

		BalanceEstimatorImage(int N, int width, int height) :
			ImageEstimator(N, width, height),
			m_image(width, height)
		{
			m_image.fill(0);
		}

		BalanceEstimatorImage(BalanceEstimatorImage const& other) :
			ImageEstimator(other)
		{
			m_image = other.m_image;
		}


		virtual void addEstimate(Spectrum const& balance_estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			addOneTechniqueEstimate(balance_estimate, tech_index, u, v, thread_safe_update);
		}

		virtual void loop() override
		{}

		virtual void addOneTechniqueEstimate(Spectrum const& balance_estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			Math::Vector<int, 2> pixel = { u * m_width, v * m_height };
			if (thread_safe_update)
				m_mutex.lock();
			m_image(pixel) += balance_estimate;
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void solve(Image::Image<Spectrum, MAJOR>& res, int iterations) override
		{
			Parallel::ParallelFor(0, m_width * m_height,
				[&](int pixel)
				{
					res[pixel] += m_image[pixel] / iterations;
				});
		}
	};

	template <class Spectrum, class Float, bool MAJOR>
	class DirectEstimatorImage : public ImageEstimator<Spectrum, Float, MAJOR>
	{
	protected:

#define ONE_CONTIGUOUS_ARRAY
#define ENABLE_DEBUG
		using solvingFloat = Float;
		using MatrixT = Eigen::Matrix<solvingFloat, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorT = Eigen::Matrix<solvingFloat, Eigen::Dynamic, 1>;
		
		using StorageUInt = unsigned int;
		using StorageFloat = Float;

		const int msize;

		size_t matTo1D(int row, int col) const {
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
#ifdef ONE_CONTIGUOUS_ARRAY
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
			const Math::Vector<int, 2> pixel(u * m_width, v * m_height);
			const int index = PixelTo1D(pixel[0], pixel[1]);
			return getPixelData(index);
		}

		std::string debugFolder()const
		{
			return RESULT_FOLDER + "debug/";
		}

#ifdef ONE_CONTIGUOUS_ARRAY
		// In Byte
		const unsigned int m_pixel_data_size;
		const unsigned int m_vector_ofsset;
		const unsigned int m_counter_offset;

		static_assert(sizeof(char) == 1);
		std::vector<char> m_data;
#else
		std::vector<StorageFloat> m_matrices;
		std::vector<StorageFloat> m_vectors;
		std::vector<StorageUInt> m_sampleCounts;
#endif

		//describes how many times more samples each technique has
		std::vector<unsigned int> m_sample_per_technique;

		std::mutex m_mutex;

	public:

#ifdef ONE_CONTIGUOUS_ARRAY
		DirectEstimatorImage(int N, int width, int height) :
			ImageEstimator(N, width, height),
			msize(N* (N + 1) / 2),
			m_pixel_data_size(msize * sizeof(StorageFloat) + spectrumDim() * m_numtechs * sizeof(StorageFloat) + m_numtechs * sizeof(StorageUInt)),
			m_vector_ofsset(msize * sizeof(StorageFloat)),
			m_counter_offset(msize * sizeof(StorageFloat) + spectrumDim() * m_numtechs * sizeof(StorageFloat)),
			m_data(std::vector<char>(width* height* m_pixel_data_size, (char)0)),
			m_sample_per_technique(std::vector<unsigned int>(m_numtechs, 1))
		{}

		DirectEstimatorImage(DirectEstimatorImage const& other) :
			ImageEstimator(other),
			msize(other.msize),
			m_pixel_data_size(other.m_pixel_data_size),
			m_vector_ofsset(other.m_vector_ofsset),
			m_counter_offset(other.m_counter_offset),
			m_data(other.m_data),
			m_sample_per_technique(other.m_sample_per_technique)
		{}

		DirectEstimatorImage(DirectEstimatorImage && other):
			ImageEstimator(other),
			msize(other.msize),
			m_pixel_data_size(other.m_pixel_data_size),
			m_vector_ofsset(other.m_vector_ofsset),
			m_counter_offset(other.m_counter_offset),
			m_data(std::move(other.m_data)),
			m_sample_per_technique(std::move(other.m_sample_per_technique))
		{}
#else
		DirectEstimatorImage(int N, int width, int height) :
			ImageEstimator(N, width, height),
			msize(N* (N + 1) / 2),
			m_sample_per_technique(std::vector<unsigned int>(m_numtechs, 1))
		{
			int res = width * height;
			m_matrices = std::vector<StorageFloat>(res * msize, (StorageFloat)0.0);
			m_vectors = std::vector<StorageFloat>(res * m_numtechs * spectrumDim(), (StorageFloat)0.0);
			m_sampleCounts = std::vector<StorageUInt>(res * m_numtechs, (StorageUInt)0);
		}

		DirectEstimatorImage(DirectEstimatorImage const& other) :
			ImageEstimator(other),
			msize(other.msize),
			m_matrices(other.m_matrices),
			m_vectors(other.m_vectors),
			m_sampleCounts(other.m_sampleCounts),
			m_sample_per_technique(other.m_sample_per_technique)
		{}

		DirectEstimatorImage(DirectEstimatorImage && other) :
			ImageEstimator(other),
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
			for (int i = 0; i < m_numtechs; ++i)
			{
				Float weight = balance_weights[i];
				if (weight < 0 || std::isnan(weight) || std::isinf(weight))
				{
					std::cout << weight << std::endl;
					__debugbreak();
				}
			}
		}

		virtual void addEstimate(Spectrum const& balance_estimate, const Float* balance_weights, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			PixelData data = getPixelData(u, v);
			if (thread_safe_update)
				m_mutex.lock();
			++data.sampleCount[tech_index];
			for (int i = 0; i < m_numtechs; ++i)
			{
				for (int j = 0; j <= i; ++j)
				{
					const int mat_index = matTo1D(i, j);
					Float tmp = balance_weights[i] * balance_weights[j];
					data.techMatrix[mat_index] = data.techMatrix[mat_index] + tmp;
				}
			}
			if (!balance_estimate.isBlack())
			{
				for (int k = 0; k < spectrumDim(); ++k)
				{
					StorageFloat* vector = data.contribVector + k * m_numtechs;
					for (int i = 0; i < m_numtechs; ++i)
					{
						Float tmp = balance_weights[i] * balance_estimate[k];
						vector[i] = vector[i] + tmp;
					}
				}
			}
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void addOneTechniqueEstimate(Spectrum const& balance_estimate, int tech_index, Float u, Float v, bool thread_safe_update = false) override
		{
			PixelData data = getPixelData(u, v);
			int mat_index = matTo1D(tech_index, tech_index);
			if (thread_safe_update)
				m_mutex.lock();
			++data.sampleCount[tech_index];
			data.techMatrix[mat_index] = data.techMatrix[mat_index] + 1.0;
			for (int k = 0; k < spectrumDim(); ++k)
			{
				StorageFloat* vector = data.contribVector + k * m_numtechs;
				vector[tech_index] = vector[tech_index] + balance_estimate[k];
			}
			if (thread_safe_update)
				m_mutex.unlock();
		}

		virtual void loop() override
		{}

		inline void fillMatrix(MatrixT& matrix, PixelData const& data, int iterations)const
		{
			for (int i = 0; i < m_numtechs; ++i)
			{
				for (int j = 0; j < i; ++j)
				{
					Float elem = data.techMatrix[matTo1D(i, j)];
					if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
					matrix(i, j) = elem;
					matrix(j, i) = elem;
				}
				matrix(i, i) = data.techMatrix[matTo1D(i, i)];
				// Unsampled samples
				size_t expected = m_sample_per_technique[i] * iterations;
				size_t actually = data.sampleCount[i];
				matrix(i, i) += (Float)(expected - actually);
			}
		}

#ifdef ENABLE_DEBUG

		virtual void saveColSum(int iterations)const
		{
			Image::Image<double, MAJOR> img(m_width * m_numtechs, m_height);
			int resolution = m_width * m_height;
			std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
			Parallel::ParallelFor(0, resolution, [&](int pixel)
				{
					int tid = Parallel::tid();
					Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
					PixelData data = getPixelData(pixel);
					MatrixT& matrix = matrices[tid];
					fillMatrix(matrix, data, iterations);
					for (int i = 0; i < m_numtechs; ++i)
					{
						size_t expected = m_sample_per_technique[i];
						Float sum = 0;
						for (int j = 0; j < m_numtechs; ++j)
						{
							sum += matrix(i, j);
						}
						sum /= iterations;
						Float diff = (expected - sum) / expected;
						img(indices[0] * m_numtechs + i, indices[1]) = diff;
					}
				});
			Image::ImWrite::writeEXR(debugFolder() + "OptiMIScolSum" + std::to_string(m_numtechs) + ".exr", img);
		}

		virtual void saveMatrices(int iterations)const
		{
			Image::Image<double, MAJOR> img(m_width * m_numtechs, m_height * m_numtechs);
			int resolution = m_width * m_height;
			std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
			Parallel::ParallelFor(0, resolution, [&](int pixel)
				{
					int tid = Parallel::tid();
					Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
					PixelData data = getPixelData(pixel);
					MatrixT& matrix = matrices[tid];
					fillMatrix(matrix, data, iterations);
					for (int i = 0; i < m_numtechs; ++i)
					{
						for (int j = 0; j < m_numtechs; ++j)
						{
							img(indices[0] * m_numtechs + i, indices[1] * m_numtechs + j) = matrix(i, j) / iterations;
						}
					}
				});
			Image::ImWrite::writeEXR(debugFolder() + "OptiMISMatrices" + std::to_string(m_numtechs) + ".exr", img);
		}

		virtual void saveVectors(int iterations)const
		{
			Image::Image<Spectrum, MAJOR> img(m_width * m_numtechs, m_height);
			int resolution = m_width * m_height;
			Parallel::ParallelFor(0, resolution, [&](int pixel)
				{
					int tid = Parallel::tid();
					Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
					PixelData data = getPixelData(pixel);
					for (int i = 0; i < m_numtechs; ++i)
					{
						for (int k = 0; k < spectrumDim(); ++k)
						{
							img(indices[0] * m_numtechs + i, indices[1])[k] = (data.contribVector + k * m_numtechs)[i] / iterations;
						}
					}
				});
			Image::ImWrite::writeEXR(debugFolder() + "OptiMISVectors" + std::to_string(m_numtechs) + ".exr", img);
		}

		virtual void saveAlphas(int iterations)const
		{
			Image::Image<Spectrum, MAJOR> img(m_width * m_numtechs, m_height);
			using Solver = Eigen::ColPivHouseholderQR<MatrixT>;
			VectorT MVector(m_numtechs);
			for (int i = 0; i < m_numtechs; ++i)	MVector(i) = m_sample_per_technique[i];
			std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
			std::vector<VectorT> vectors = Parallel::preAllocate(VectorT(m_numtechs));
			std::vector<Solver> solvers = Parallel::preAllocate(Solver(m_numtechs, m_numtechs));
			int resolution = m_width * m_height;
			Parallel::ParallelFor(0, resolution, [&](int pixel)
				{
					int tid = Parallel::tid();
					Math::Vector<int, 2> indices = Image::Image<double, MAJOR>::indices(pixel, m_width, m_height);
					MatrixT& matrix = matrices[tid];
					VectorT& vector = vectors[tid];
					Solver& solver = solvers[tid];

					PixelData data = getPixelData(pixel);

					bool matrix_done = false;

					for (int k = 0; k < spectrumDim(); ++k)
					{
						bool isZero = true;
						const StorageFloat* contribVector = data.contribVector + k * m_numtechs;
						for (int i = 0; i < m_numtechs; ++i)
						{
							Float elem = contribVector[i];
							if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
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
							for (int i = 0; i < m_numtechs; ++i)
							{
								vector[i] *= m_sample_per_technique[i];
								img(indices[0] * m_numtechs + i, indices[1])[k] = vector[i];
							}
						}
					}

				});
			Image::ImWrite::writeEXR(debugFolder() + "OptiMISAlphas" + std::to_string(m_numtechs) + ".exr", img);
		}

		virtual void debug(int iterations, bool col_sum, bool matrix, bool vec, bool alpha)const override
		{
#ifdef ENABLE_DEBUG
			if (col_sum)
				saveColSum(iterations);
			if (matrix)
				saveMatrices(iterations);
			if (vec)
				saveVectors(iterations);
			if (alpha)
				saveAlphas(iterations);
#else 
			std::cout<<"Warning! Image Direct Estimator Debug disabled. Nothing happens!"
#endif
		}
#endif

		virtual void solve(Image::Image<Spectrum, MAJOR>& res, int iterations) override
		{
			using Solver = Eigen::ColPivHouseholderQR<MatrixT>;
			VectorT MVector(m_numtechs);
			for (int i = 0; i < m_numtechs; ++i)	MVector(i) = m_sample_per_technique[i];
			std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
			std::vector<VectorT> vectors = Parallel::preAllocate(VectorT(m_numtechs));
			std::vector<Solver> solvers = Parallel::preAllocate(Solver(m_numtechs, m_numtechs));
			int resolution = m_width * m_height;
			Parallel::ParallelFor(0, resolution, [&](int pixel)
				{
					int tid = Parallel::tid();
					MatrixT& matrix = matrices[tid];
					VectorT& vector = vectors[tid];
					Solver& solver = solvers[tid];

					Spectrum color = 0;

					PixelData data = getPixelData(pixel);

					bool matrix_done = false;

					for (int k = 0; k < spectrumDim(); ++k)
					{
						bool isZero = true;
						const StorageFloat* contribVector = data.contribVector + k * m_numtechs;
						for (int i = 0; i < m_numtechs; ++i)
						{
							Float elem = contribVector[i];
							if (std::isnan(elem) || std::isinf(elem) || elem < 0)	elem = 0;
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
					res[pixel] += color;
				});
		}

#undef ONE_CONTIGUOUS_ARRAY
#undef ENABLE_DEBUG
	};

} // namespace MIS