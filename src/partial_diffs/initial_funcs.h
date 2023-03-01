#pragma once
#include <types/basic_types.h>
namespace ComputationalPhysics::Schemes
{

using S = Types::BasicTypes::scalar;
struct XIV102
{
    [[nodiscard]] static S calc(S const x, S const L) noexcept
    {
        return std::sin(4 * M_PI * x / L);
    };
};

} // namespace ComputationalPhysics::Schemes