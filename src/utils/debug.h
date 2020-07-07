#pragma once

namespace MIS
{

std::string debugFolder()const
{
    return "debug/";
}

virtual void saveColSum(int iterations)const
{
    Image::Image<double, ROW_MAJOR> img(m_width * m_numtechs, m_height);
    int resolution = m_width * m_height;
    std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
    Parallel::ParallelFor(0, resolution, [&](int pixel)
        {
            int tid = Parallel::tid();
            Math::Vector<int, 2> indices = Image::Image<double, ROW_MAJOR>::indices(pixel, m_width, m_height);
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
    Image::Image<double, ROW_MAJOR> img(m_width * m_numtechs, m_height * m_numtechs);
    int resolution = m_width * m_height;
    std::vector<MatrixT> matrices = Parallel::preAllocate(MatrixT(m_numtechs, m_numtechs));
    Parallel::ParallelFor(0, resolution, [&](int pixel)
        {
            int tid = Parallel::tid();
            Math::Vector<int, 2> indices = Image::Image<double, ROW_MAJOR>::indices(pixel, m_width, m_height);
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
    Image::Image<Spectrum, ROW_MAJOR> img(m_width * m_numtechs, m_height);
    int resolution = m_width * m_height;
    Parallel::ParallelFor(0, resolution, [&](int pixel)
        {
            int tid = Parallel::tid();
            Math::Vector<int, 2> indices = Image::Image<double, ROW_MAJOR>::indices(pixel, m_width, m_height);
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
    Image::Image<Spectrum, ROW_MAJOR> img(m_width * m_numtechs, m_height);
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
            Math::Vector<int, 2> indices = Image::Image<double, ROW_MAJOR>::indices(pixel, m_width, m_height);
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

}