#pragma once
#include <types/basic_types.h>
namespace ComputationalPhysics::Utils::SpecialFuncs
{

template<typename F, typename S>
auto fmap(F const &f, std::vector<S> const &s) noexcept
{
    using V = decltype(F::calc())
}

} // namespace ComputationalPhysics::Utils::SpecialFunctions
