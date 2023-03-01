#pragma once
#include <types/basic_types.h>
namespace ComputationalPhysics::Schemes
{

using S = Types::BasicTypes::scalar;
template<typename V>
struct Node
{
    V y;
    S x;
};

using S = Types::BasicTypes::scalar;
template<typename V>
struct Grid
{

    Grid(S const L, S const h) noexcept
        : h(h)
    {
        Nx = 1u + static_cast<unsigned>(L / h);

        S xt = static_cast<S>(0);
        for (unsigned i = 0; i < Nx; ++i)
        {
            grid.template emplace_back(Node<V>{static_cast<V>(0), xt});
            xt += h;
        }
    }

    std::vector<Node<V>> grid;
    Types::BasicTypes::scalar h;
    unsigned Nx;

};

} // namespace ComputationalPhysics::Schemes