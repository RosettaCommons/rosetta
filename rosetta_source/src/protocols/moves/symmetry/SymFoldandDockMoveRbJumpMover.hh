// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_moves_symmetry_SymFoldandDockMoveRbJumpMover_hh

#define INCLUDED_protocols_moves_symmetry_SymFoldandDockMoveRbJumpMover_hh

// Unit headers
#include <protocols/moves/symmetry/SymFoldandDockMoveRbJumpMover.fwd.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers

namespace protocols {
namespace moves {
namespace symmetry {
///////////////////////////////////////////////////////////////////////////////
class SymFoldandDockMoveRbJumpMover : public Mover
{
public:

	// default constructor

	SymFoldandDockMoveRbJumpMover();

	~SymFoldandDockMoveRbJumpMover(){}

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

};

} // symmetry
} // moves
} // rosetta
#endif
