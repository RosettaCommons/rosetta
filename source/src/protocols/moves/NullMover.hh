// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   NullMover.hh
///
/// @brief
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_moves_NullMover_hh
#define INCLUDED_protocols_moves_NullMover_hh

#include <protocols/moves/NullMover.fwd.hh>
#include <protocols/moves/MoveMapMover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {


class NullMover : public moves::MoveMapMover
{
public:
	/// @brief
	/// 	empty constructor fills values with the values
	///		read in from the commandline
	NullMover();
	virtual void apply( core::pose::Pose &  );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
  virtual void set_movemap( core::kinematics::MoveMapCOP ) {}
  virtual core::kinematics::MoveMapCOP movemap() const { return core::kinematics::MoveMapCOP( new core::kinematics::MoveMap ); }
	virtual ~NullMover();
	virtual void test_move( core::pose::Pose &  ){};

private:
};


} // moves
} // protocols

#endif

