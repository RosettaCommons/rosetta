#pragma once

#include <memory>
#include <map>
#include <string>

namespace ui {
namespace network {

class Argument;
using ArgumentSP  = std::shared_ptr< Argument >;
using ArgumentCSP = std::shared_ptr< Argument const >;

using Arguments = std::map<std::string, ArgumentSP>;
using ArgumentsSP  = std::shared_ptr< Arguments >;
using ArgumentsCSP = std::shared_ptr< Arguments const >;

} // namespace network
} // namespace ui
