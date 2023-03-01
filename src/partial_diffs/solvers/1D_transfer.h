#pragma once
#include <partial_diffs/grid.h>
#include <partial_diffs/schemes/explicit/leftAngle.h>
namespace ComputationalPhysics::PartialDiffs::Solvers::Transfer1D
{

using S = Types::BasicTypes::scalar;
template<typename F>
struct Solve1DTransferLeftAngle // with periodic boundary conditions
{
    [[nodiscard]] static std::vector<std::vector<Schemes::Node<S>>> calcScheme( Schemes::Grid               <S>       &grid
                                                                              , Schemes::Explicit::LeftAngle<S> const &leftAngle
                                                                              , S const T
                                                                              , S const dt
                                                                              ) noexcept
    {
        unsigned const N = 1u + static_cast<unsigned>(T / dt);
        unsigned const Nx = grid.Nx;

        S const h = grid.h;

        for (unsigned i = 0; i < Nx; ++i)
            grid.grid[i].y = F::calc(h * i, grid.Nx * h);

        std::vector<std::vector<Schemes::Node<S>>> results;

        results.template emplace_back(grid.grid);

        for (unsigned t = 1; t < N; ++t)
        {
            for(unsigned i = 1; i < Nx; ++i)
            {
                grid.grid[i].y = leftAngle.calc(results[t - 1][i].y, results[t - 1][i - 1].y);
            }
            grid.grid[0].y = grid.grid.back().y;
            results.template emplace_back(grid.grid);
        }

        return results;
    }
};

} // namespace ComputationalPhysics::PartialDiffs::Solvers