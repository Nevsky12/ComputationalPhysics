#pragma once
#include <functional>
#include <types/core_types.h>
#include <utils/tensor_ops/kroneker_product.h>
namespace ComputationalPhysics::Integrators::Implicit
{

using Val = ComputationalPhysics::Types::BasicTypes::vec;
using Mat = ComputationalPhysics::Types::BasicTypes::mat;
using State = ComputationalPhysics::Types::CoreTypes::State;
using Arg   = ComputationalPhysics::Types::CoreTypes::  arg;
using Res      = ComputationalPhysics::Types::BasicTypes::vecSeriesSTL;
using ResState = ComputationalPhysics::Types::CoreTypes::stateSeries;

using Step = decltype(std::declval<Arg>() - std::declval<Arg>());

template<unsigned S>
using ButcherTable = Types::CoreTypes::ButcherTable<S>;

template<unsigned S>
[[nodiscard]] Res implicitRK( State const &init
                            , Step  const h
                            , Val   const eps
                            , ButcherTable<S> const &bTable
                            , std::function<Val(Arg const, Val const&)> const &rightPart
                            , std::function<Mat(Arg const, Val const&)> const &jacobi
                            ) noexcept;

template<unsigned S>
[[nodiscard]] ResState implicitRKFinite( State const &init
                                       , Step  const h
                                       , Arg   const broad
                                       , ButcherTable<S> const &bT
                                       , std::function<Val(Arg const, Val const&)> const &rightPart
                                       , std::function<Mat(Arg const, Val const&)> const &jacobi
                                       ) noexcept;

template<unsigned S>
[[nodiscard]] ResState implicitRKFinite( State const &init
                                        , Step const h
                                        , Arg  const broad
                                        , ButcherTable<S> const &bTable
                                        , std::function<Val(Arg const, Val const&)> const &rightPart
                                        , std::function<Mat(Arg const, Val const&)> const &jacobi
                                        ) noexcept
{
    ResState result;
    auto cInit = init;
    auto const &[b, c, A] = bTable;
    result.emplace_back(State{cInit.state, cInit.t});

    unsigned const N = static_cast<unsigned>(broad / h);
    unsigned const M = cInit.state.size();
    Types::BasicTypes::mat I = Types::BasicTypes::mat::Identity(S * M, S * M);

    for (unsigned k = 0u; k < N; ++k)
    {
        auto const &solve = [=]( Types::BasicTypes::scalar const t
                               , Types::BasicTypes::vec const &y
                               , Types::BasicTypes::mat const &Val
                               , Types::BasicTypes::mat const &J
                               ) noexcept
        {
            Types::BasicTypes::mat cVal = Val;
            Types::BasicTypes::mat const &JJ = I - h * Utils::Tensors::kronekerProduct(A, J).eval();
            auto const &LU = JJ.fullPivLu();
            auto const &newtonStep= [&]( Types::BasicTypes::scalar const t
                                       , Types::BasicTypes::   vec const &y
                                       , Types::BasicTypes::   vec const &initState
                                       ) noexcept
            {
                auto const &calcF = [&]( Types::BasicTypes::   vec const &stageDer
                                       , Types::BasicTypes::scalar const t
                                       , Types::BasicTypes::   vec const &y
                                       ) noexcept
                {
                    Types::BasicTypes::mat stageDerNext(S, M);
                    Types::BasicTypes::mat stageDerC = stageDer.reshaped(S, M);
                    Types::BasicTypes::vec stageVal(M);
                    for (unsigned i = 0u; i < S; ++i)
                    {
                        stageVal = y + h * Types::BasicTypes::vec{(A.row(i) * stageDerC).reshaped()};
                        stageDerNext.row(i) = rightPart(t + h * c(i), stageVal);
                    }
                    Types::BasicTypes::vec const &dY = Types::BasicTypes::vec{stageDerNext.reshaped()} - stageDer;
                    return dY;
                };

                Types::BasicTypes::vec const &rhs = calcF(initState, t, y);
                Types::BasicTypes::vec const &d = LU.solve(rhs);
                Types::BasicTypes::vec const &res = initState + d;
                return res;
            };

            // 3 per one newton stage
            for (unsigned i = 0u; i < 3u; ++i)
                cVal = newtonStep(t, y, cVal);
            return cVal;
        };

        Types::BasicTypes::mat const &J = jacobi(cInit.t, cInit.state);
        Types::BasicTypes::vec const &stageDer0 =  rightPart(cInit.t, cInit.state);
        Types::BasicTypes::mat  stageDer(S, M);
        for (unsigned i = 0u; i < S; ++i)
            stageDer.row(i) = stageDer0;
        Types::BasicTypes::mat const &stageVal = solve(cInit.t, cInit.state, Types::BasicTypes::vec{stageDer.reshaped()}, J);
        cInit.t += h;
        cInit.state +=  h * stageVal.reshaped(M, S) * b;
        result.emplace_back(State{cInit.state, cInit.t});
    }
    return result;
}

} // ComputationalPhysics::Integrators::Robust
