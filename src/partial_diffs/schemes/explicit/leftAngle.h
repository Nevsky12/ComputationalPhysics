#pragma once
#include <types/basic_types.h>
namespace ComputationalPhysics::Schemes::Explicit
{

using S = Types::BasicTypes::scalar;
template<typename Y>
/*
 *   du / dt + du / dx = 0
 */
struct LeftAngle
{
    LeftAngle(S const CFL) noexcept
        : CFL(        CFL)
    {}

    [[nodiscard]] Y calc(Y const &current, Y const &previous) const noexcept
    {
        return current + CFL * (current - previous);
    }

    S CFL;
};

} // namespace ComputationalPhysics::Schemes
