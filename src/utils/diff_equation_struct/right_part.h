#pragma once
#include "../../types/core_types.h"
namespace ComputationalPhysics::Utils::DiffEquation
{

using State = ComputationalPhysics::Types::CoreTypes::State;
using Arg   = ComputationalPhysics::Types::CoreTypes::  arg;
using Val   = ComputationalPhysics::Types::CoreTypes::  vec;

[[nodiscard]] Val rightPart( Arg   arg
                           , State state
                           ) noexcept;

} // ComputationalPhysics::Utils::DiffEquation
