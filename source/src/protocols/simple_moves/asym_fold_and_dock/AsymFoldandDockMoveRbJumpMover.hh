// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockMoveRbJumpMover_hh
#define INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockMoveRbJumpMover_hh

// Unit headers
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockMoveRbJumpMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace simple_moves {
namespace asym_fold_and_dock {
///////////////////////////////////////////////////////////////////////////////
class AsymFoldandDockMoveRbJumpMover : public moves::Mover
{
public:

	// default constructor

	AsymFoldandDockMoveRbJumpMover( core::Size chain_start );

	~AsymFoldandDockMoveRbJumpMover(){}

	void apply( core::pose::Pose & pose );
	void find_new_jump_residue( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Size chain_start_;

};

}
} // asym_fold_and_dock
} // rosetta
#endif
