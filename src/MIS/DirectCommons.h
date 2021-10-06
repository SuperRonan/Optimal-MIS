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
	};
}