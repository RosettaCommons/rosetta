// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/ChangeFoldTreeMover.hh
/// @brief ChangeFoldTreeMover
/// @author

#ifndef INCLUDED_protocols_moves_ChangeFoldTreeMover_hh
#define INCLUDED_protocols_moves_ChangeFoldTreeMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/ChangeFoldTreeMover.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace moves {


class ChangeFoldTreeMover : public Mover {
public:
	// default constructor
	ChangeFoldTreeMover():
		Mover()
	{};

	ChangeFoldTreeMover(
		core::kinematics::FoldTree fold_tree
	) :
		Mover(), fold_tree_(fold_tree)
	{}
	
	void set_foldtree(core::kinematics::FoldTree fold_tree){
		fold_tree_ = fold_tree;
	}
	
	/// @brief Apply the stored fold tree to the pose
	virtual void apply( core::pose::Pose & pose ) { pose.fold_tree( fold_tree_ ); }
	virtual std::string get_name() const { return "ChangeFoldTreeMover"; }

private:
	core::kinematics::FoldTree fold_tree_;
};


} // moves
} // protocols


#endif
