#include <gtest/gtest.h>
#include <solvers/three_diagonal_solver.h>
#include <utils/special_matrix/boundary_matrix.h>
#include <iostream>
#include <fstream>

TEST(SOLVERS, BOUNDARY_PROBLEM)
{
    ComputationalPhysics::Types::BasicTypes::scalar const h = 0.1;
    ComputationalPhysics::Types::BasicTypes::scalar const lBound = 0.;
    ComputationalPhysics::Types::BasicTypes::scalar const rBound = M_PI;
    ComputationalPhysics::Types::BasicTypes::scalar const lVal = 0.;
    ComputationalPhysics::Types::BasicTypes::scalar const rVal = M_PI * M_PI;
    unsigned const N = static_cast<unsigned>((rBound - lBound) / h) + 1u;


    auto const f = +[](ComputationalPhysics::Types::BasicTypes::scalar const xn) noexcept
                                  -> ComputationalPhysics::Types::BasicTypes::scalar
    {
        return 2.
             - 6. * xn
             + 2. * std::pow(xn, 3.)
             + (xn * xn - 3.) * exp(xn) * sin(xn) * (1. + cos(xn))
             + cos(xn) * (exp(xn) + (xn * xn - 1.) + std::pow(xn, 4.) - 3. * xn * xn);
    };
    auto const an = +[](ComputationalPhysics::Types::BasicTypes::scalar const xn) noexcept
                                    -> ComputationalPhysics::Types::BasicTypes::scalar
    {
        return xn * xn - 3.;
    };
    auto const bn = +[](ComputationalPhysics::Types::BasicTypes::scalar const xn) noexcept
                                     -> ComputationalPhysics::Types::BasicTypes::scalar
    {
        return (xn * xn - 3.) * std::cos(xn);
    };


    std::pair<ComputationalPhysics::Utils::SpecialMatrix::ThreeDiagonalMatrix, ComputationalPhysics::Types::BasicTypes::vecSTL>
            A = ComputationalPhysics::Utils::SpecialMatrix::ExpandedMatrixForLinearBoundaryValueProblem3( lBound, rBound
                                                                                                        , lVal  , rVal
                                                                                                        , an, bn, f
                                                                                                        , N
                                                                                                        );

    auto const &result = ComputationalPhysics::Solvers::solveThreeDiagonal(A.first, A.second);

    ComputationalPhysics::Types::BasicTypes::scalar xi = lBound;
    ComputationalPhysics::Types::BasicTypes::vecSTL x;
    for (unsigned i = 1u; i < N + 1u; ++i)
    {
        x.emplace_back(   xi);
        xi += h;
    }
    x.emplace_back(rBound);

    std::ofstream out("boundary.csv");
    unsigned i = 0u;
    out << "X" << '\t' << "Y" << '\n';
    for (auto const stuff: result)
    {
        std::cout << "X: " << x[i] << '\t' << "/ Y: " <<  stuff << std::endl;
        out << x[i] << '\t';
        out << stuff << '\n';
        ++i;
    }
    out.close();
}
