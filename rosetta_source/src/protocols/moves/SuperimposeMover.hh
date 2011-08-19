// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SuperimposeMover.hh
/// @brief
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_moves_SuperimposeMover_hh
#define INCLUDED_protocols_moves_SuperimposeMover_hh

// Unit headers
#include <protocols/moves/SuperimposeMover.fwd.hh>
#include <protocols/moves/Mover.hh> // we need to store a pose
#include <core/pose/Pose.hh>

namespace protocols {
namespace moves {

class SuperimposeMover : public moves::Mover {

public:
	/// @brief
	/// 	empty constructor
	SuperimposeMover();

	SuperimposeMover( core::pose::Pose const & pose );

	~SuperimposeMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void set_reference_pose( core::pose::Pose const & pose );

private:

	core::pose::Pose ref_pose_;

};

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_SuperimposeMover_HH
