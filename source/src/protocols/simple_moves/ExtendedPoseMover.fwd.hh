// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ExtendedPoseMover.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_ExtendedPoseMover_FWD_HH
#define INCLUDED_protocols_simple_moves_ExtendedPoseMover_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace simple_moves {

class ExtendedPoseMover;
typedef utility::pointer::shared_ptr<ExtendedPoseMover> ExtendedPoseMoverOP;
typedef utility::pointer::shared_ptr<ExtendedPoseMover const> ExtendedPoseMoverCOP;

}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_protocols_simple_moves_ExtendedPoseMover_FWD_HH
