#pragma once

#include <memory>

namespace ui {
namespace task {

class Functor;
using FunctorSP  = std::shared_ptr< Functor >;
using FunctorCSP = std::shared_ptr< Functor const >;


} // namespace task
} // namespace ui
