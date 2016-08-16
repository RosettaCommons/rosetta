// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LocalPosition.fwd.hh
/// @brief definition of the LocalPosition class
/// @author

#ifndef INCLUDED_protocols_environment_LocalPosition_fwd_hh
#define INCLUDED_protocols_environment_LocalPosition_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>
#include <utility/vector1.hh>

// Package headers

namespace core {
namespace environment {

class LocalPosition;
typedef utility::pointer::shared_ptr< LocalPosition > LocalPositionOP;
typedef utility::pointer::shared_ptr< LocalPosition const > LocalPositionCOP;

typedef utility::vector1< LocalPositionOP > LocalPositions;

} // environment
} // core

#endif //INCLUDED_protocols_moves_LocalPosition_fwd_HH
