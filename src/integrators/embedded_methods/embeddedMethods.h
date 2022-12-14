#pragma once
#include <functional>
#include <iostream>
#include <types/core_types.h>
namespace ComputationalPhysics::Integrators::Embedded
{

using Val = ComputationalPhysics::Types::BasicTypes::vec;
using State = ComputationalPhysics::Types::CoreTypes::State;
using Arg   = ComputationalPhysics::Types::CoreTypes::  arg;
using Res      = ComputationalPhysics::Types::BasicTypes::vecSeriesSTL;
using ResState = ComputationalPhysics::Types::CoreTypes::stateSeries;

using Step = decltype(std::declval<Arg>() - std::declval<Arg>());

template<unsigned S>
using EmbeddedTable = Types::CoreTypes::EmbeddedTable<S>;

template<unsigned S, unsigned P>
[[nodiscard]] ResState embeddedMethod( State const &init
                                     , Step  const h
                                     , Arg   const broad
                                     , Step  const tolerance
                                     , EmbeddedTable<S> const &eT
                                     , std::function<Val(Arg const,  Val const&)> const &rightPart
                                     , std::function<Arg(Val const&, Val const&)> const &errFunc
                                     ) noexcept
{
    ResState result;

    auto cInit = init;
    auto const A  = eT.A;
    auto const b  = eT.b;
    auto const c  = eT.c;
    auto const bp = eT.bp;
    result.emplace_back(State{cInit.state, cInit.t});
    Types::BasicTypes::scalar hn = h;
    unsigned const M = cInit.state.size();

    auto const &calcRhs = [&]( Types::BasicTypes::   vec const &stageDer
                             , Types::BasicTypes::scalar const t
                             , Types::BasicTypes::   vec const &y
                             ) noexcept -> Types::BasicTypes::mat
    {
        Types::BasicTypes::mat K = Types::BasicTypes::mat::Zero(S, M);
        K.row(0) = stageDer;
        Types::BasicTypes::vec stageVal(M);
        for (unsigned i = 1u; i < S; ++i)
        {
            stageVal = y + hn * Types::BasicTypes::vec{(A.row(i) * K).reshaped()};
            K.row(i) = rightPart(t + hn * c(i), stageVal);
        }
        return K;
    };

    while (cInit.t < broad)
    {
        cInit.t += hn;
        Types::BasicTypes::vec const &k0 = rightPart(cInit.t, cInit.state);
        Types::BasicTypes::mat const &rhs = calcRhs(k0, cInit.t, cInit.state);

        Types::BasicTypes::vec const &y     = cInit.state + hn * rhs.transpose() * b;
        Types::BasicTypes::vec const &yPerm = cInit.state + hn * rhs.transpose() * bp;
        Types::BasicTypes::scalar const delta = errFunc(y, yPerm);
        Types::BasicTypes::scalar const s = std::pow(hn * tolerance / (2 * broad * delta), 1 / (P + 1u));

        if (s >= static_cast<Types::BasicTypes::scalar>(2.))
        {
            cInit.state = y;
            result.emplace_back(State{cInit.state, cInit.t});
            hn *= static_cast<Types::BasicTypes::scalar>(2.);
            if(cInit.t + hn > broad)
                hn = broad - cInit.t;
        } else if (s >= static_cast<Types::BasicTypes::scalar>(1.))
        {
            cInit.state = y;
            result.emplace_back(State{cInit.state, cInit.t});
            if(cInit.t + hn > broad)
                hn = broad - cInit.t;
        } else if (s < static_cast<Types::BasicTypes::scalar>(1.))
            hn /= static_cast<Types::BasicTypes::scalar>(2.);
    }
    return result;
}

} // ComputationalPhysics::Integrators::Embedded