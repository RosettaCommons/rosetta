// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/MoveMapMover.hh
/// @brief A MoveMapMover is nothing more than a mover that has a function that sets its movemap.
///
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_moves_MoveMapMover_hh
#define INCLUDED_protocols_moves_MoveMapMover_hh

// Unit Headers
#include <protocols/moves/MoveMapMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/kinematics/MoveMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility Headers

// C++ Headers

namespace protocols {
namespace moves {


class MoveMapMover : public Mover {
	typedef Mover Parent;
public:
	MoveMapMover(): Parent() {}
	MoveMapMover( std::string const& name ) : Parent( name ) {}
	MoveMapMover( MoveMapMover const & ) = default;

	~MoveMapMover() override = default;

	virtual void set_movemap( core::kinematics::MoveMapCOP ) = 0;

	/// @brief Get a movemap for the given pose
	/// (As the movemap may vary based on the Pose (e.g. MoveMapFactory), we need the Pose information)
	virtual core::kinematics::MoveMapCOP movemap( core::pose::Pose const & ) const = 0;

	virtual void initialize( core::pose::Pose& ) {}


}; // end MoveMapMover base class

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_MoveMapMover_HH
