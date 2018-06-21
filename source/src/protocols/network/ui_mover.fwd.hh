// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/ui_mover.fwd.hh
/// @brief: UIMover, apply to send in-progress results to UI
///
/// @author Sergey Lyskov

#pragma once

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace network {

class UIMover;
using UIMoverSP  = utility::pointer::shared_ptr< UIMover >;
using UIMoverCSP = utility::pointer::shared_ptr< UIMover const >;

// deprecated
using UIMoverOP  = UIMoverSP;
using UIMoverCOP = UIMoverCSP;

} // namespace network
} // namespace protocols
