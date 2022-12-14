#pragma once
#include <iostream>
#include <types/core_types.h>
namespace ComputationalPhysics::Utils::SpecialMatrix
{

struct ThreeDiagonalMatrix
{
    explicit inline ThreeDiagonalMatrix(unsigned const size) noexcept: core_(size){};

    [[nodiscard]] static ThreeDiagonalMatrix     Zero(unsigned const size) noexcept;
    [[nodiscard]] static ThreeDiagonalMatrix Identity(unsigned const size) noexcept;
    [[nodiscard]] static ThreeDiagonalMatrix ThreeNumbers( unsigned const size
                                                         , Types::BasicTypes::scalar const a
                                                         , Types::BasicTypes::scalar const b
                                                         , Types::BasicTypes::scalar const c
                                                         ) noexcept;

    [[nodiscard]] inline       double & operator()(int const i, int const j)       noexcept {return core_[i][j];};
    [[nodiscard]] inline const double & operator()(int const i, int const j) const noexcept {return core_[i][j];};

    [[nodiscard]] inline unsigned rows() const noexcept {return core_.size();};

    inline void fill_row( unsigned const ind
                        , Types::BasicTypes::scalar const a
                        , Types::BasicTypes::scalar const b
                        , Types::BasicTypes::scalar const c
                        ) noexcept;

    void check_diagonal_domimance() const noexcept;

    Types::CoreTypes::threeDiagSeries core_;
};

ThreeDiagonalMatrix ThreeDiagonalMatrix::Zero(unsigned const size) noexcept
{
    ThreeDiagonalMatrix result(size);
    for (unsigned i = 0u; i < size; ++i)
        result.core_[i] = {0., 0., 0.};
    return result;
}
ThreeDiagonalMatrix ThreeDiagonalMatrix::Identity(unsigned const size) noexcept
{
    ThreeDiagonalMatrix result(size);
    for (unsigned i = 0u; i < size; ++i)
        result.core_[i] = {0., 1., 0.};
    return result;
}

ThreeDiagonalMatrix ThreeDiagonalMatrix::ThreeNumbers( unsigned const size
                                                     , Types::BasicTypes::scalar const a
                                                     , Types::BasicTypes::scalar const b
                                                     , Types::BasicTypes::scalar const c
                                                     ) noexcept
{
    ThreeDiagonalMatrix result(size);
    for (unsigned i = 0u; i < size; ++i)
        result.core_[i] = {a, b, c};
    return result;
}

inline void ThreeDiagonalMatrix::fill_row( unsigned const ind
                                         , Types::BasicTypes::scalar const a
                                         , Types::BasicTypes::scalar const b
                                         , Types::BasicTypes::scalar const c
                                         ) noexcept
{
    core_[ind][0] = a;
    core_[ind][1] = b;
    core_[ind][2] = c;
}

void ThreeDiagonalMatrix::check_diagonal_domimance() const noexcept
{
    bool one_strict_inequality = false;
    bool non_strict_inequality = true;
    for(unsigned i = 0u; i < this->rows(); ++i)
    {
        if (abs(core_[i][1]) <  abs(core_[i][0]) + abs(core_[i][2]))
        {
            std::cout << "The sufficient condition of diagonal dominance is not fulfilled in row: "
                      << i
                      << std::endl;
            non_strict_inequality = !non_strict_inequality;
            break;
        }
        if (abs(core_[i][1]) > abs(core_[i][0]) + abs(core_[i][2]))
            one_strict_inequality = !one_strict_inequality;
    }
    if (one_strict_inequality && non_strict_inequality)
        std::cout << "The sufficient condition of diagonal dominance is fulfilled!"
                     "Strict inequality is satisfied in at least one row! "
                  << std::endl;
}

} // ComputationalPhysics::Utils::SpecialMatrix
