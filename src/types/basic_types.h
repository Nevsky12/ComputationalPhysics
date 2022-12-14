#pragma once
#include <Eigen/Dense>
#include <vector>
namespace ComputationalPhysics::Types::BasicTypes
{

using scalar = double;
using mat = Eigen::MatrixXd;
using vec = Eigen::VectorXd;
using vecSTL = std::vector<scalar>;

using vecSeriesSTL = std::vector<vec>;

template<unsigned N, typename T>
using arrSeries = std::vector<std::array<T, N>>;

} // ComputationalPhysics::Types::BasicTypes