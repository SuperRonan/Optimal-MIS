#pragma once

#include <Eigen/Dense>
#include <cstdint>

namespace MIS
{
	template <class Float>	
	struct LinearSystem
	{
		using MatrixT = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
		MatrixT tech_matrix;
		// One column per spectrum channel
		MatrixT contrib_vectors, alphas;

		LinearSystem(int N, int S):
			tech_matrix(N, N),
			contrib_vectors(N, S),
			alphas(N, S)
		{}
	};
}