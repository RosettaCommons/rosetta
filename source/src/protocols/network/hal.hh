// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/hal.hh
/// @brief: HAL implementaion
///
/// @author Sergey Lyskov

#ifdef ZEROMQ

#pragma once

#include <json.hpp>

#include <functional>
#include <string>

namespace protocols {
namespace network {


using SpecificationCallBack = std::function< std::string(void) >;
using ExecutionerCallBack   = std::function< json(json const &) >;

/// start HAL listener
/// use main function `argc` and `argv` if you want auto-restart/abort functionality to work properly
void hal(SpecificationCallBack const &, ExecutionerCallBack const &, int argc, char * argv []);



} // namespace network
} // namespace protocols

#endif // ZEROMQ
