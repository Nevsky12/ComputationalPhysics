#pragma once
#include <utils/special_matrix/three_diagonal_matrix.h>
namespace ComputationalPhysics::Solvers
{

Types::BasicTypes::vecSTL solveThreeDiagonal( Utils::SpecialMatrix::ThreeDiagonalMatrix const &A
                                            , Types::BasicTypes::vecSTL const &b
                                            ) noexcept
{
    A.check_diagonal_domimance();
    unsigned const N = b.size();
    Types::BasicTypes::vecSTL res(N);

    std::vector<std::array<double, 2>> pr(N - 1u);

    pr[0][0] = A(0,2) / A(0,1);
    pr[0][1] = b[0]        /  A(0, 1);

    for (unsigned i = 1u; i < N - 1u; ++i)
    {
        int j = static_cast<int>(i);
        pr[i][0] =  A(j,2)
                / ( A(j,1)
                  - A(j,0) * pr[j - 1][0]
                  );
        pr[i][1] = ( b[i]
                   - A(j, 0) * pr[j - 1][1]
                   )
                 / ( A(j,1)
                   - A(j,0) * pr[j - 1][0]
                   );

    }

    int n = static_cast<int>(N);

    res.back() = ( b.back()
                 - A(n - 1, 0) * pr[n - 2][1]
                 )
               / ( A(n - 1, 1)
                 - A(n - 1,0) * pr[n - 2][0]
                 );

    for (int i = n - 2; i >= 0; --i)
        res[i] = pr[i][1] - pr[i][0] * res[i+1];
    return res;

}

[[nodiscard]] Types::BasicTypes::vecSTL threeDiagonalSolver( Types::BasicTypes::   mat const &A
                                                           , Types::BasicTypes::vecSTL const &b
                                                           ) noexcept
{
    unsigned const size = A.rows();
    Types::BasicTypes::arrSeries<2, Types::BasicTypes::scalar> runThroughCoefficients(size);

    runThroughCoefficients[0][0] = -A(0, 1) / A(0, 0);
    runThroughCoefficients[0][1] = b[0] / A(0, 0);

    for (unsigned i = 1; i < size - 1; ++i)
    {
        runThroughCoefficients[i][0] = -A(i, i + 1)
                                       / ( runThroughCoefficients[i - 1][0] * A(i, i - 1)
                                           + A(i, i)
                                       );

        runThroughCoefficients[i][1] = ( b[i]
                                         - runThroughCoefficients[i - 1][1] * A(i, i - 1) )
                                       / ( runThroughCoefficients[i - 1][0] * A(i, i - 1)
                                           + A(i, i)
                                       );
    }

    Types::BasicTypes::vecSTL res(size);
    res.back() = (b[size - 1] - A(size - 1, 0) * runThroughCoefficients[size - 1][1]) / (A(size - 1, 0) * runThroughCoefficients[size - 2][0] + A(size - 1, 2));

    for (int i = static_cast<int>(size) - 2; i >= 0; i--)
        res[i] = runThroughCoefficients[i][0] * res[i + 1] + runThroughCoefficients[i][1];
    return res;
}

} // ComputationalPhysics::Solvers
