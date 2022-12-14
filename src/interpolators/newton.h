#pragma once
#include <types/core_types.h>
namespace ComputationalPhysics::Interpolators
{

Types::BasicTypes::scalar newton( Types::BasicTypes::vecSTL const &x
                                , Types::BasicTypes::vecSTL const &y
                                , Types::BasicTypes::scalar const val
                                ) noexcept
{
    Types::BasicTypes::scalar result = y[0];
    unsigned const n = y.size();

    for(unsigned i = 1u; i < n; ++i)
    {
        Types::BasicTypes::scalar f_i = 0.;
        for (unsigned j = 0u; j <= i; ++j)
        {
            Types::BasicTypes::scalar dx = 1.;
            for (unsigned k = 0u; k <= i; ++k)
                if (k != j)
                    dx *= (x[j] - x[k]);
            f_i += y[j] / dx;
        }
        for (unsigned j = 0u; j < i; ++j)
            f_i *= (val - x[j]);
        result += f_i;
    }
    return result;
}

} // ComputationalPhysics::Interpolators