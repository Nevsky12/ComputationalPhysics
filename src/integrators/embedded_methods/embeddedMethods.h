#pragma once
#include <functional>
#include <iostream>
#include <types/core_types.h>
namespace ComputationalPhysics::Integrators::Embedded
{

using Val = ComputationalPhysics::Types::BasicTypes::vec;
using Mat = ComputationalPhysics::Types::BasicTypes::mat;
using State = ComputationalPhysics::Types::CoreTypes::State;
using Arg   = ComputationalPhysics::Types::CoreTypes::  arg;
using Res      = ComputationalPhysics::Types::BasicTypes::vecSeriesSTL;
using ResState = ComputationalPhysics::Types::CoreTypes::stateSeries;

using Step = decltype(std::declval<Arg>() - std::declval<Arg>());

template<unsigned S>
using EmbeddedTable = Types::CoreTypes::EmbeddedTable<S>;

template<unsigned S>
[[nodiscard]] ResState embeddedMethod( State const &init
                                     , Step  const h
                                     , Arg   const broad
                                     , Step  const tolerance
                                     , EmbeddedTable<S> const &eT
                                     , std::function<Val(Arg const, Val const&)> const &rightPart
                                     , std::function<Arg(Arg const, Arg const )> const &errFunc
                                     ) noexcept
{
    ResState result;

    auto cInit = init;
    auto const A  = eT.A;
    auto const b  = eT.b;
    auto const c  = eT.c;
    auto const bp = eT.bp;
    result.emplace_back(State{cInit.state, cInit.t});
    Types::BasicTypes::scalar tn = cInit.t;
    Types::BasicTypes::scalar hn = h;
    unsigned const M = cInit.state.size();

    auto const &calcRhs = [&]( Types::BasicTypes::   vec const &stageDer
                                                                 , Types::BasicTypes::scalar const t
                                                                 , Types::BasicTypes::   vec const &y
                                                                 ) noexcept -> Types::BasicTypes::vec
    {
        Types::BasicTypes::mat K = Types::BasicTypes::mat::Zero(S, M);
        K.row(0) = stageDer;
        Types::BasicTypes::vec stageDerNext(M);
        Types::BasicTypes::vec stageVal(M);
        for (unsigned i = 1u; i < S; ++i)
        {
            stageVal = y + hn * Types::BasicTypes::vec{(A.row(i) * K).reshaped()};
            K.row(i) = rightPart(t + hn * c(i), stageVal);
        }
        Types::BasicTypes::vec const dY = Types::BasicTypes::vec{K.colw.reshaped()} - stageDer;
        return dY;
    };

    auto const calcStep = [=](Types::BasicTypes::vec const &err) noexcept
                                           -> Types::BasicTypes::scalar
    {
        return err.squaredNorm();
    };

    while (cInit.t < broad)
    {
        Types::BasicTypes::vec const k0 = rightPart(cInit.t, cInit.state);
        Types::BasicTypes::vec const rhs = calcRhs(k0, cInit.t, cInit.state);
        Types::BasicTypes::vec const y     = cInit.state + hn * b  * rhs;
        Types::BasicTypes::vec const yPerm = cInit.state + hn * bp * rhs;
        Types::BasicTypes::scalar const err = (y - yPerm).squaredNorm();
        Types::BasicTypes::scalar const delta = errFunc(err, tolerance);

        if (err < tolerance)
        {
            cInit.t += hn;
            cInit.state += hn * rhs;
        }
        if        (delta <= static_cast<Types::BasicTypes::scalar>(0.1))
        {
            hn *= static_cast<Types::BasicTypes::scalar>(0.1);
        } else if (delta >= static_cast<Types::BasicTypes::scalar>(4. ))
        {
            hn *= static_cast<Types::BasicTypes::scalar>(4. );
        } else {
            hn *= delta;
        }
        if (cInit.t + hn > broad)
            hn = broad - cInit.t;


        std::cout << "Y: " << std::endl << cInit.state << std::endl;
        std::cout << "X: " << tn << " / hn: " << hn << std::endl;


        cInit.t += hn;
    }
    return result;
}

} // ComputationalPhysics::Integrators::Embedded