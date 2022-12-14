#pragma once
#include <Eigen/Dense>
namespace ComputationalPhysics::Utils::Tensors
{

Eigen::MatrixXd kronekerProduct( Eigen::MatrixXd const &m1
                               , Eigen::MatrixXd const &m2
                               ) noexcept
{
    Eigen::MatrixXd m3( m1.rows() * m2.rows()
                      , m1.cols() * m2.cols()
                      );
    for (unsigned i = 0; i < m1.cols(); ++i)
        for (unsigned j = 0; j < m1.rows(); ++j)
            m3.block(i * m2.rows(), j * m2.cols(),
                        m2.rows(),    m2.cols()) =  m1(i,j) * m2;
    return m3;
}

} // namespace ComputationalPhysics::Utils::Tensors
