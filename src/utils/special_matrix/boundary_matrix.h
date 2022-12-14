#pragma once
#include <functional>
#include "three_diagonal_matrix.h"
namespace ComputationalPhysics::Utils::SpecialMatrix
{

std::pair<ThreeDiagonalMatrix, Types::BasicTypes::vecSTL>
         ExpandedMatrixForLinearBoundaryValueProblem3( Types::BasicTypes::scalar const lbx
                                                     , Types::BasicTypes::scalar const rbx
                                                     , Types::BasicTypes::scalar const lby
                                                     , Types::BasicTypes::scalar const rby
                                                     , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar const)> const &a
                                                     , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar const)> const &b
                                                     , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar const)> const &f
                                                     , unsigned const N
                                                     ) noexcept
{
        auto const h = (rbx - lbx) / N;
        ThreeDiagonalMatrix TDM = ThreeDiagonalMatrix(N + 1);
        unsigned const RTDM = TDM.rows();

        TDM(0                         , 1) = static_cast<Types::BasicTypes::scalar>(1.);
        TDM(static_cast<int>(RTDM) - 1, 1) = static_cast<Types::BasicTypes::scalar>(1.);
        for(unsigned i = 1u; i < RTDM - 1u; ++i)
            TDM.fill_row(i, Types::CoreTypes::Calc<Types::CoreTypes::ColumnIdx::first >::calc(a, b, h, lbx + h * i),
                                Types::CoreTypes::Calc<Types::CoreTypes::ColumnIdx::second>::calc(a, b, h, lbx + h * i),
                                Types::CoreTypes::Calc<Types::CoreTypes::ColumnIdx::third >::calc(a, b, h, lbx + h * i)
                        );

        Types::BasicTypes::vecSTL y(N + 1);
        y.front() = lby;
        for(unsigned i = 1u; i < y.size() - 1u; ++i)
            y[i] = f(lbx + h * i);
        y.back() = rby;

        return {TDM, y};
}

} // ComputationalPhysics::Utils::SpecialMatrix