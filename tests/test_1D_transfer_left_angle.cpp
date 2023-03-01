#include <gtest/gtest.h>
#include <partial_diffs/solvers/1D_transfer.h>
#include <partial_diffs/initial_funcs.h>
#include <array>
#include <iostream>
#include <fstream>

TEST(SCHEMES, TRANSFER_1D_LEFT_ANGLE)
{
    std::ofstream results;

    ComputationalPhysics::Types::BasicTypes::scalar const L = 20;
    ComputationalPhysics::Types::BasicTypes::scalar const T = 18;
    ComputationalPhysics::Types::BasicTypes::scalar const h  = 0.5;
    std::array<ComputationalPhysics::Types::BasicTypes::scalar , 3> arrayCFL = {0.6, 1., 1.01};

    using F = ComputationalPhysics::Schemes::XIV102;

    ComputationalPhysics::Schemes::Grid<ComputationalPhysics::Types::BasicTypes::scalar> grid(L, h);

    ComputationalPhysics::PartialDiffs::Solvers::Transfer1D::Solve1DTransferLeftAngle<F> solve1DTransferLeftAngle;

    for (auto const CFL: arrayCFL)
    {
        ComputationalPhysics::Types::BasicTypes::scalar const dt = CFL * h;
        std::string cfl = "CFL=";
        cfl += std::to_string(CFL);
        cfl += " / dt=";
        cfl += std::to_string(dt);
        results.open(cfl);
        results << "x,y,";

        std::cout << "CFL: " << CFL << std::endl;
        ComputationalPhysics::Schemes::Explicit::LeftAngle<ComputationalPhysics::Types::BasicTypes::scalar> const leftAngle(CFL);

        bool dtGot = false;
        results << std::to_string(dt) << "\n";
        auto const result = solve1DTransferLeftAngle.calcScheme(grid, leftAngle, T, dt);

        ComputationalPhysics::Types::BasicTypes::scalar t = 0.;
        for (auto const &layer: result)
        {
            std::cout << "t: " << t << std::endl;
            for (auto const &[y, x]: layer)
            {
                results << std::to_string(x) << "," << std::to_string(y) << "\n";
                std::cout << "x: " << x << " / " << "y: " << y << std::endl;
            }
            std::cout << "---------------------" << std::endl;
            t += dt;
        }
        results.close();
    }
}