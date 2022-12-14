#pragma once
#include "basic_types.h"
namespace ComputationalPhysics::Types::CoreTypes
{

using scalar = Types::BasicTypes::scalar;
using vec = Types::BasicTypes::vec;
using mat = Types::BasicTypes::mat;

using arg = double;

struct State
{
    vec state;
    arg t;
};
using stateSeries = std::vector<State>;

template<unsigned S>
struct ButcherTable
{
    vec b;
    vec c;
    mat A;
};
template<unsigned S>
struct EmbeddedTable: ButcherTable<S>{vec bp;};

struct SplineCoefficients
{
    scalar a;
    scalar b;
    scalar c;
    scalar d;
    scalar x;
};

using splineSeries = std::vector<SplineCoefficients>;

using threeDiagSeries = std::vector<std::array<double, 3>>;

enum ColumnIdx
{
    first = 1,
    second = 2,
    third = 3,
};

template<unsigned ColumnIdx>
struct Calc
{
    [[nodiscard]] static Types::BasicTypes::scalar calc( std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &a
                                                       , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &b
                                                       , Types::BasicTypes::scalar const h
                                                       , Types::BasicTypes::scalar const x
                                                       ) noexcept;
};

template<>
Types::BasicTypes::scalar Calc<ColumnIdx::first>::calc( std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &a
                                                        , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &b
                                                        , Types::BasicTypes::scalar const h
                                                        , Types::BasicTypes::scalar const x
                                                        ) noexcept
{
    return static_cast<Types::BasicTypes::scalar>(1.) / (h * h)
                                          - a(x) / (static_cast<Types::BasicTypes::scalar>(2.) * h);
}
template<>
Types::BasicTypes::scalar Calc<ColumnIdx::second>::calc( std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &a
                                                        , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &b
                                                        , Types::BasicTypes::scalar const h
                                                        , Types::BasicTypes::scalar const x
                                                        ) noexcept
{
    return static_cast<Types::BasicTypes::scalar>(-2.) / (h * h)
                                                       + b(x);
}
template<>
Types::BasicTypes::scalar Calc<ColumnIdx::third>::calc( std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &a
                                                      , std::function<Types::BasicTypes::scalar(Types::BasicTypes::scalar)> const &b
                                                      , Types::BasicTypes::scalar const h
                                                      , Types::BasicTypes::scalar const x
                                                      ) noexcept
{
    return static_cast<Types::BasicTypes::scalar>(1.) / (h * h)
                                                      + a(x)
                                                      / (static_cast<Types::BasicTypes::scalar>(2.) * h);
}

} // ComputationalPhysics::Types::CoreTypes
